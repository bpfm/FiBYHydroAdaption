#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"
#include "kernel.h"

#ifdef BG_STELLAR_EVOLUTION
#include "bg_proto.h"
#endif



/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct densdata_in
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat Hsml;
  int NodeList[NODELISTLENGTH];
}
 *DensDataIn, *DensDataGet;


static struct densdata_out
{
  MyLongDouble Rho;
  MyLongDouble Div, Rot[3];
  MyLongDouble DhsmlDensity;
  MyLongDouble Ngb;

#ifdef PRESSURE_ENTROPY_SPH
  MyLongDouble DhsmlPressure;
  MyLongDouble WeightedDensity;
#endif

#ifdef BG_METALSMOOTHING
  MyLongDouble MetalsSmoothed[BG_NELEMENTS];
  MyLongDouble MetallicitySmoothed;
/* #ifdef BG_DUST */
/*   MyLongDouble DusticitySmoothed; */
/* #endif */
#ifdef BG_SNIA_IRON
  MyLongDouble IronFromSNIaSmoothed;
#endif
#ifdef BG_MOL_NETWORK
  MyLongDouble x_H2_Smoothed;
  MyLongDouble x_H2p_Smoothed;
  MyLongDouble x_HD_Smoothed;
#endif
#endif

#if defined(BLACK_HOLES)
  MyLongDouble SmoothedEntr;
  MyLongDouble GasVel[3];
#endif

}
 *DensDataResult, *DensDataOut;


/*! \file density.c 
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and some
 *  auxiliary quantities are computed.  There is also functionality that
 *  corrects the smoothing length if needed.
 */


/*! This function computes the local density for each active SPH particle, the
 * number of neighbours in the current smoothing radius, and the divergence
 * and rotation of the velocity field.  The pressure is updated as well.  If a
 * particle with its smoothing region is fully inside the local domain, it is
 * not exported to the other processors. The function also detects particles
 * that have a number of neighbours outside the allowed tolerance range. For
 * these particles, the smoothing length is adjusted accordingly, and the
 * density() computation is called again.  Note that the smoothing length is
 * not allowed to fall below the lower bound set by MinGasHsml (this may mean
 * that one has to deal with substantially more than normal number of
 * neighbours.)
 */
void density(void)
{
  MyFloat *Left, *Right;

  int i, j, ndone, ndone_flag, npleft, dt_step, dummy, iter = 0;

  int ngrp, sendTask, recvTask, place, nexport, nimport;

  long long ntot;

  double dmax1, dmax2, fac;

  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait;

  double dt_entr, tstart, tend, t0, t1;

  double desnumngb;

#if defined(BG_METALSMOOTHING)
  int k;
#endif

  Left = (MyFloat *) mymalloc(NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc(NumPart * sizeof(MyFloat));

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(density_isactive(i))
	{
	  Left[i] = Right[i] = 0;

#ifdef BLACK_HOLES
	  P[i].SwallowID = 0;
#endif
/* #if ((defined(BLACK_HOLES) && defined(BH_THERMALFEEDBACK)) || defined(BG_SNII_THERMAL_FEEDBACK)) && defined(FLTROUNDOFFREDUCTION) */
/*           if(P[i].Type == 0) */
/*             SphP[i].i.dInjected_Energy = SphP[i].i.Injected_Energy; */
/* #endif */
	}
    }

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct densdata_in) + sizeof(struct densdata_out) +
					     sizemax(sizeof(struct densdata_in),
						     sizeof(struct densdata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));


  CPU_Step[CPU_DENSMISC] += measure_time();
  t0 = second();

  desnumngb = All.DesNumNgb;

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      i = FirstActiveParticle;	/* begin with this index */

      do
	{
	  for(j = 0; j < NTask; j++)
	    {
	      Send_count[j] = 0;
	      Exportflag[j] = -1;
	    }

	  /* do local particles and prepare export list */
	  tstart = second();
	  for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	    {
	      if(density_isactive(i))
		{
		  if(density_evaluate(i, 0, &nexport, Send_count) < 0)
		    break;
		}
	    }
	  tend = second();
	  timecomp1 += timediff(tstart, tend);

#ifdef MYSORT
	  mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#else
	  qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#endif

	  tstart = second();

	  MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

	  tend = second();
	  timewait1 += timediff(tstart, tend);

	  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	    {
	      Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
	      nimport += Recv_count[j];

	      if(j > 0)
		{
		  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
		  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
		}
	    }

	  DensDataGet = (struct densdata_in *) mymalloc(nimport * sizeof(struct densdata_in));
	  DensDataIn = (struct densdata_in *) mymalloc(nexport * sizeof(struct densdata_in));

	  /* prepare particle data for export */
	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      DensDataIn[j].Pos[0] = P[place].Pos[0];
	      DensDataIn[j].Pos[1] = P[place].Pos[1];
	      DensDataIn[j].Pos[2] = P[place].Pos[2];
	      DensDataIn[j].Hsml = PPP[place].Hsml;

	      memcpy(DensDataIn[j].NodeList,
		     DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));

#if defined(BLACK_HOLES) || defined (BG_SFR)
	      if(P[place].Type != 0)
		{
		  DensDataIn[j].Vel[0] = 0;
		  DensDataIn[j].Vel[1] = 0;
		  DensDataIn[j].Vel[2] = 0;
		}
	      else
#endif
		{
		  DensDataIn[j].Vel[0] = SphP[place].VelPred[0];
		  DensDataIn[j].Vel[1] = SphP[place].VelPred[1];
		  DensDataIn[j].Vel[2] = SphP[place].VelPred[2];
		}
	    }

	  /* exchange particle data */
	  tstart = second();
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&DensDataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A,
				   &DensDataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }
	  tend = second();
	  timecommsumm1 += timediff(tstart, tend);

	  myfree(DensDataIn);
	  DensDataResult = (struct densdata_out *) mymalloc(nimport * sizeof(struct densdata_out));
	  DensDataOut = (struct densdata_out *) mymalloc(nexport * sizeof(struct densdata_out));


	  /* now do the particles that were sent to us */

	  tstart = second();
	  for(j = 0; j < nimport; j++)
	    density_evaluate(j, 1, &dummy, &dummy);
	  tend = second();
	  timecomp2 += timediff(tstart, tend);

	  if(i < 0)
	    ndone_flag = 1;
	  else
	    ndone_flag = 0;

	  tstart = second();
	  MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	  tend = second();
	  timewait2 += timediff(tstart, tend);


	  /* get the result */
	  tstart = second();
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&DensDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B,
				   &DensDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}

	    }
	  tend = second();
	  timecommsumm2 += timediff(tstart, tend);


	  /* add the result to the local particles */
	  tstart = second();
	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      PPP[place].n.dNumNgb += DensDataOut[j].Ngb;

	      if(P[place].Type == 0)
		{
		  SphP[place].d.dDensity += DensDataOut[j].Rho;
		  SphP[place].v.dDivVel += DensDataOut[j].Div;
		  SphP[place].r.dRot[0] += DensDataOut[j].Rot[0];
		  SphP[place].r.dRot[1] += DensDataOut[j].Rot[1];
		  SphP[place].r.dRot[2] += DensDataOut[j].Rot[2];
		  SphP[place].h.dDhsmlDensityFactor += DensDataOut[j].DhsmlDensity;

#ifdef PRESSURE_ENTROPY_SPH
		  SphP[place].cky.dWeightedDensity += DensDataOut[j].WeightedDensity;
		  SphP[place].dky.dDhsmlPressure += DensDataOut[j].DhsmlPressure;
#endif

#if defined(BG_METALSMOOTHING)
		  for(k = 0; k < BG_NELEMENTS; k++)
		    SphP[place].MetalsSmoothed[k] += (MyFloat) DensDataOut[j].MetalsSmoothed[k];
		  SphP[place].MetallicitySmoothed += (MyFloat) DensDataOut[j].MetallicitySmoothed;
/* #ifdef BG_DUST */
/*                   SphP[place].DusticitySmoothed += (MyFloat) DensDataOut[j].DusticitySmoothed; */
/* #endif */
#ifdef BG_SNIA_IRON
		  SphP[place].IronFromSNIaSmoothed += (MyFloat) DensDataOut[j].IronFromSNIaSmoothed;
#endif
#ifdef BG_MOL_NETWORK
		  SphP[place].x_H2_Smoothed += (MyFloat) DensDataOut[j].x_H2_Smoothed;
		  SphP[place].x_H2p_Smoothed += (MyFloat) DensDataOut[j].x_H2p_Smoothed;
		  SphP[place].x_HD_Smoothed += (MyFloat) DensDataOut[j].x_HD_Smoothed;
#endif
#endif

		}
#ifdef BLACK_HOLES
	      if(P[place].Type == 5)
		{
		  P[place].b1.dBH_Density += DensDataOut[j].Rho;
		  P[place].b2.dBH_Entropy += DensDataOut[j].SmoothedEntr;
		  P[place].b3.dBH_SurroundingGasVel[0] += DensDataOut[j].GasVel[0];
		  P[place].b3.dBH_SurroundingGasVel[1] += DensDataOut[j].GasVel[1];
		  P[place].b3.dBH_SurroundingGasVel[2] += DensDataOut[j].GasVel[2];
		}
#endif
	    }
	  tend = second();
	  timecomp1 += timediff(tstart, tend);


	  myfree(DensDataOut);
	  myfree(DensDataResult);
	  myfree(DensDataGet);
	}
      while(ndone < NTask);

#ifdef FLTROUNDOFFREDUCTION
      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	if(density_isactive(i))
	  {
	    PPP[i].n.NumNgb = FLT(PPP[i].n.dNumNgb);

	    if(P[i].Type == 0)
	      {
		SphP[i].d.Density = FLT(SphP[i].d.dDensity);
		SphP[i].h.DhsmlDensityFactor = FLT(SphP[i].h.dDhsmlDensityFactor);
#ifdef PRESSURE_ENTROPY_SPH
		SphP[i].cky.WeightedDensity = FLT(SphP[i].cky.dWeightedDensity);
		SphP[i].dky.DhsmlPressure = FLT(SphP[i].dky.dDhsmlPressure);
#endif
		SphP[i].v.DivVel = FLT(SphP[i].v.dDivVel);
		for(j = 0; j < 3; j++)
		  SphP[i].r.Rot[j] = FLT(SphP[i].r.dRot[j]);
	      }

#ifdef BLACK_HOLES
	    if(P[i].Type == 5)
	      {
		P[i].b1.BH_Density = FLT(P[i].b1.dBH_Density);
		P[i].b2.BH_Entropy = FLT(P[i].b2.dBH_Entropy);
		for(j = 0; j < 3; j++)
		  P[i].b3.BH_SurroundingGasVel[j] = FLT(P[i].b3.dBH_SurroundingGasVel[j]);
	      }
#endif
	  }
#endif


      /* do final operations on results */
      tstart = second();
      for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
	{
	  if(density_isactive(i))
	    {
	      if(P[i].Type == 0)
		{
		  if(SphP[i].d.Density > 0)
		    {
		      SphP[i].h.DhsmlDensityFactor *= PPP[i].Hsml / (NUMDIMS * SphP[i].d.Density);
		      if(SphP[i].h.DhsmlDensityFactor > -0.9)	/* note: this would be -1 if only a single particle at zero lag is found */
			SphP[i].h.DhsmlDensityFactor = 1 / (1 + SphP[i].h.DhsmlDensityFactor);
		      else
			SphP[i].h.DhsmlDensityFactor = 1;

		      SphP[i].r.CurlVel = sqrt(SphP[i].r.Rot[0] * SphP[i].r.Rot[0] +
					       SphP[i].r.Rot[1] * SphP[i].r.Rot[1] +
					       SphP[i].r.Rot[2] * SphP[i].r.Rot[2]) / SphP[i].d.Density;

		      SphP[i].v.DivVel /= SphP[i].d.Density;

#ifdef PRESSURE_ENTROPY_SPH
                            SphP[i].dky.DhsmlPressure *= PPP[i].Hsml / (NUMDIMS * SphP[i].EntropyVarPred * SphP[i].d.Density);
                            SphP[i].cky.WeightedDensity /= SphP[i].EntropyVarPred;
#endif

#ifdef BG_METALSMOOTHING
		      /* smoothed mass of elements */
		      for(k = 0; k < BG_NELEMENTS; k++)
			SphP[i].MetalsSmoothed[k] *= (P[i].Mass / SphP[i].d.Density);

		      /* smoothed mass fraction of elements heavier that Helium */
		      SphP[i].MetallicitySmoothed /= SphP[i].d.Density;
/* #ifdef BG_DUST */
/*                       SphP[i].DusticitySmoothed /= SphP[i].d.Density; */
/* #endif */
#ifdef BG_SNIA_IRON
		      /* smoothed mass of Iron from SNIa */
		      SphP[i].IronFromSNIaSmoothed *= (P[i].Mass / SphP[i].d.Density);
#endif
#ifdef BG_MOL_NETWORK
		      SphP[i].x_H2_Smoothed /= SphP[i].d.Density;
		      SphP[i].x_H2p_Smoothed /= SphP[i].d.Density;
		      SphP[i].x_HD_Smoothed /= SphP[i].d.Density;
#endif
#endif
		    }

		  /*
		  dt_step = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0);
		  dt_entr = (All.Ti_Current - (P[i].Ti_begstep + dt_step / 2)) * All.Timebase_interval;
		  */

#ifndef PRESSURE_ENTROPY_SPH
		  SphP[i].Pressure =
		    SphP[i].EntropyPred * pow(SphP[i].d.Density, GAMMA);
#else
		  SphP[i].Pressure =
		    SphP[i].EntropyPred * pow(SphP[i].cky.WeightedDensity, GAMMA);
#endif
		}

#ifdef BLACK_HOLES
	      if(P[i].Type == 5)
		{
		  if(P[i].b1.BH_Density > 0)
		    {
		      P[i].b2.BH_Entropy /= P[i].b1.BH_Density;
		      P[i].b3.BH_SurroundingGasVel[0] /= P[i].b1.BH_Density;
		      P[i].b3.BH_SurroundingGasVel[1] /= P[i].b1.BH_Density;
		      P[i].b3.BH_SurroundingGasVel[2] /= P[i].b1.BH_Density;
		    }
		}
#endif

	      /* now check whether we had enough neighbours */
#ifdef BLACK_HOLES
	      desnumngb = All.DesNumNgb;
#endif

	      /*
              if(P[i].Type == 0)
                {
                  PPP[i].n.NumNgb = NORM_COEFF * PPP[i].Hsml * PPP[i].Hsml * SphP[i].d.Density / P[i].Mass;
#ifndef TWODIMS
                  PPP[i].n.NumNgb *= PPP[i].Hsml;
#endif
                }
	      */

	      if(PPP[i].n.NumNgb < (desnumngb - All.MaxNumNgbDeviation) ||
		 (PPP[i].n.NumNgb > (desnumngb + All.MaxNumNgbDeviation)
		  && PPP[i].Hsml > (1.01 * All.MinGasHsml)))
		{
		  /* need to redo this particle */
		  npleft++;

		  if(Left[i] > 0 && Right[i] > 0)
		    if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
		      {
			/* this one should be ok */
			npleft--;
			P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
			continue;
		      }

		  if(PPP[i].n.NumNgb < (desnumngb - All.MaxNumNgbDeviation))
		    Left[i] = DMAX(PPP[i].Hsml, Left[i]);
		  else
		    {
		      if(Right[i] != 0)
			{
			  if(PPP[i].Hsml < Right[i])
			    Right[i] = PPP[i].Hsml;
			}
		      else
			Right[i] = PPP[i].Hsml;
		    }

		  if(iter >= MAXITER - 10)
		    {
		      printf
			("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			 i, ThisTask, (int) P[i].ID, PPP[i].Hsml, Left[i], Right[i],
			 (float) PPP[i].n.NumNgb, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
		      fflush(stdout);
		    }

		  if(Right[i] > 0 && Left[i] > 0)
		    PPP[i].Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
		  else
		    {
		      if(Right[i] == 0 && Left[i] == 0)
			endrun(8188);	/* can't occur */

		      if(Right[i] == 0 && Left[i] > 0)
			{
			  if(P[i].Type == 0 && fabs(PPP[i].n.NumNgb - desnumngb) < 0.5 * desnumngb)
			    {
			      fac = 1 - (PPP[i].n.NumNgb -
					 desnumngb) / (NUMDIMS * PPP[i].n.NumNgb) *
				SphP[i].h.DhsmlDensityFactor;

			      if(fac < 1.26)
				PPP[i].Hsml *= fac;
			      else
				PPP[i].Hsml *= 1.26;
			    }
			  else
			    PPP[i].Hsml *= 1.26;
			}

		      if(Right[i] > 0 && Left[i] == 0)
			{
			  if(P[i].Type == 0 && fabs(PPP[i].n.NumNgb - desnumngb) < 0.5 * desnumngb)
			    {
			      fac = 1 - (PPP[i].n.NumNgb -
					 desnumngb) / (NUMDIMS * PPP[i].n.NumNgb) *
				SphP[i].h.DhsmlDensityFactor;

			      if(fac > 1 / 1.26)
				PPP[i].Hsml *= fac;
			      else
				PPP[i].Hsml /= 1.26;
			    }
			  else
			    PPP[i].Hsml /= 1.26;
			}
		    }

		  if(PPP[i].Hsml < All.MinGasHsml)
		    PPP[i].Hsml = All.MinGasHsml;
		}
	      else
		P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
	    }
	}
      tend = second();
      timecomp1 += timediff(tstart, tend);

      sumup_large_ints(1, &npleft, &ntot);

      if(ntot > 0)
	{
	  iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in neighbour iteration in density()\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	}
    }
  while(ntot > 0);


  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
  myfree(Right);
  myfree(Left);

  /* mark as active again */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].TimeBin < 0)
      P[i].TimeBin = -P[i].TimeBin - 1;

  /* collect some timing information */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_DENSCOMPUTE] += timecomp;
  CPU_Step[CPU_DENSWAIT] += timewait;
  CPU_Step[CPU_DENSCOMM] += timecomm;
  CPU_Step[CPU_DENSMISC] += timeall - (timecomp + timewait + timecomm);
}


/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int density_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n;
  int startnode, numngb, numngb_inbox, listindex = 0;

  double h, h2, fac, hinv, hinv3, hinv4;

  double wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  double dvx, dvy, dvz;

  MyLongDouble rho;
  MyLongDouble divv, rotv[3];
  MyLongDouble weighted_numngb;
  MyLongDouble dhsmlrho;

#ifdef PRESSURE_ENTROPY_SPH
  MyLongDouble dhsmlpressure;
  MyLongDouble weighteddensity;
  dhsmlpressure = weighteddensity = 0;
#endif

#ifdef BLACK_HOLES
  MyLongDouble gasvel[3];
#endif

  MyDouble *pos;
  MyFloat *vel;

  static MyFloat veldummy[3] = { 0, 0, 0 };

#if defined(BLACK_HOLES)
  MyLongDouble smoothentr;
  smoothentr = 0;
#endif

#ifdef BG_METALSMOOTHING
  int k;

  MyLongDouble metals_smoothed[BG_NELEMENTS], metallicity_smoothed = 0;
/* #ifdef BG_DUST */
/*   MyLongDouble dusticity_smoothed = 0; */
/* #endif */
#ifdef BG_SNIA_IRON
  MyLongDouble iron_from_snia_smoothed = 0;
#endif
#ifdef BG_MOL_NETWORK
  MyLongDouble x_h2_smoothed = 0, x_h2p_smoothed = 0;
  MyLongDouble x_hd_smoothed = 0;
#endif

  for(k = 0; k < BG_NELEMENTS; k++)
    metals_smoothed[k] = 0;
#endif


  divv = rotv[0] = rotv[1] = rotv[2] = 0;
#ifdef BLACK_HOLES
  gasvel[0] = gasvel[1] = gasvel[2] = 0;
#endif
  rho = weighted_numngb = dhsmlrho = 0;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
      if(P[target].Type == 0)
	{
	  vel = SphP[target].VelPred;
	}
      else
	vel = veldummy;
    }
  else
    {
      pos = DensDataGet[target].Pos;
      vel = DensDataGet[target].Vel;
      h = DensDataGet[target].Hsml;
    }


  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv; // / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;


  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = DensDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  numngb = 0;

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb_inbox = ngb_treefind_variable(pos, h, target, &startnode, mode, nexport, nsend_local);

	  if(numngb_inbox < 0)
	    return -1;

	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = Ngblist[n];
#ifdef BLACK_HOLES
	      if(P[j].Mass == 0)
		continue;
#endif
	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	      if(dx > boxHalf_X)
		dx -= boxSize_X;
	      if(dx < -boxHalf_X)
		dx += boxSize_X;
	      if(dy > boxHalf_Y)
		dy -= boxSize_Y;
	      if(dy < -boxHalf_Y)
		dy += boxSize_Y;
	      if(dz > boxHalf_Z)
		dz -= boxSize_Z;
	      if(dz < -boxHalf_Z)
		dz += boxSize_Z;
#endif
	      r2 = dx * dx + dy * dy + dz * dz;

	      if(r2 < h2)
		{
		  numngb++;

		  r = sqrt(r2);

		  u = r * hinv;

/* 		  if(u < 0.5) */
/* 		    { */
/* 		      wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u); */
/* 		      dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4); */
/* 		    } */
/* 		  else */
/* 		    { */
/* 		      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u); */
/* 		      dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u); */
/* 		    } */

		  wk = kernel_wk(u, hinv3);

		  dwk = kernel_dwk(u, hinv4);

		  mass_j = P[j].Mass;

		  rho += FLT(mass_j * wk);

		  weighted_numngb += FLT(NORM_COEFF * wk / hinv3);	/* 4.0/3 * PI = 4.188790204786 */

		  dhsmlrho += FLT(-mass_j * (NUMDIMS * hinv * wk + u * dwk));

#ifdef PRESSURE_ENTROPY_SPH
		  dhsmlpressure += FLT(-mass_j * SphP[j].EntropyVarPred * (NUMDIMS * hinv * wk + u * dwk));
		  weighteddensity += FLT(mass_j * SphP[j].EntropyVarPred * wk);
#endif

#ifdef BG_METALSMOOTHING
		  for(k = 0; k < BG_NELEMENTS; k++)
		    metals_smoothed[k] += SphP[j].Metals[k] * wk;

		  metallicity_smoothed += SphP[j].Metallicity * mass_j * wk;
/* #ifdef BG_DUST */
/*                   dusticity_smoothed += SphP[j].Dusticity * mass_j * wk; */
/* #endif */
#ifdef BG_SNIA_IRON
		  iron_from_snia_smoothed += SphP[j].IronFromSNIa * wk;
#endif
#ifdef BG_MOL_NETWORK
		  x_h2_smoothed += SphP[j].x_H2 * mass_j * wk;
		  x_h2p_smoothed += SphP[j].x_H2p * mass_j * wk;
		  x_hd_smoothed += SphP[j].x_HD * mass_j * wk;
#endif
#endif


#ifdef BLACK_HOLES
		  smoothentr += FLT(mass_j * wk * SphP[j].Entropy);
		  gasvel[0] += FLT(mass_j * wk * SphP[j].VelPred[0]);
		  gasvel[1] += FLT(mass_j * wk * SphP[j].VelPred[1]);
		  gasvel[2] += FLT(mass_j * wk * SphP[j].VelPred[2]);
#endif

		  if(r > 0)
		    {
		      fac = mass_j * dwk / r;

		      dvx = vel[0] - SphP[j].VelPred[0];
		      dvy = vel[1] - SphP[j].VelPred[1];
		      dvz = vel[2] - SphP[j].VelPred[2];

		      divv += FLT(-fac * (dx * dvx + dy * dvy + dz * dvz));

		      rotv[0] += FLT(fac * (dz * dvy - dy * dvz));
		      rotv[1] += FLT(fac * (dx * dvz - dz * dvx));
		      rotv[2] += FLT(fac * (dy * dvx - dx * dvy));
		    }
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = DensDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  if(mode == 0)
    {
      PPP[target].n.dNumNgb = weighted_numngb;
      if(P[target].Type == 0)
	{
	  SphP[target].d.dDensity = rho;
	  SphP[target].v.dDivVel = divv;
	  SphP[target].h.dDhsmlDensityFactor = dhsmlrho;
	  SphP[target].r.dRot[0] = rotv[0];
	  SphP[target].r.dRot[1] = rotv[1];
	  SphP[target].r.dRot[2] = rotv[2];

#ifdef PRESSURE_ENTROPY_SPH
	  SphP[target].dky.dDhsmlPressure = dhsmlpressure;
	  SphP[target].cky.dWeightedDensity = weighteddensity;
#endif

#ifdef BG_METALSMOOTHING
	  for(k = 0; k < BG_NELEMENTS; k++)
	    SphP[target].MetalsSmoothed[k] = (MyFloat) metals_smoothed[k];
	  SphP[target].MetallicitySmoothed = (MyFloat) metallicity_smoothed;
/* #ifdef BG_DUST */
/*           SphP[target].DusticitySmoothed = (MyFloat) dusticity_smoothed; */
/* #endif */
#ifdef BG_SNIA_IRON
	  SphP[target].IronFromSNIaSmoothed = (MyFloat) iron_from_snia_smoothed;
#endif
#ifdef BG_MOL_NETWORK
	  SphP[target].x_H2_Smoothed = (MyFloat) x_h2_smoothed;
	  SphP[target].x_H2p_Smoothed = (MyFloat) x_h2p_smoothed;
	  SphP[target].x_HD_Smoothed = (MyFloat) x_hd_smoothed;
#endif
#endif
	}
#ifdef BLACK_HOLES
      P[target].b1.dBH_Density = rho;
      P[target].b2.dBH_Entropy = smoothentr;
      P[target].b3.dBH_SurroundingGasVel[0] = gasvel[0];
      P[target].b3.dBH_SurroundingGasVel[1] = gasvel[1];
      P[target].b3.dBH_SurroundingGasVel[2] = gasvel[2];
#endif
    }
  else
    {
      DensDataResult[target].Rho = rho;
      DensDataResult[target].Div = divv;
      DensDataResult[target].Ngb = weighted_numngb;
      DensDataResult[target].DhsmlDensity = dhsmlrho;
      DensDataResult[target].Rot[0] = rotv[0];
      DensDataResult[target].Rot[1] = rotv[1];
      DensDataResult[target].Rot[2] = rotv[2];

#ifdef PRESSURE_ENTROPY_SPH
      DensDataResult[target].DhsmlPressure = dhsmlpressure;
      DensDataResult[target].WeightedDensity = weighteddensity;
#endif

#ifdef BG_METALSMOOTHING
      for(k = 0; k < BG_NELEMENTS; k++)
	DensDataResult[target].MetalsSmoothed[k] = metals_smoothed[k];
      DensDataResult[target].MetallicitySmoothed = metallicity_smoothed;
/* #ifdef BG_DUST */
/*       DensDataResult[target].DusticitySmoothed = dusticity_smoothed; */
/* #endif */
#ifdef BG_SNIA_IRON
      DensDataResult[target].IronFromSNIaSmoothed = iron_from_snia_smoothed;
#endif
#ifdef BG_MOL_NETWORK
      DensDataResult[target].x_H2_Smoothed = x_h2_smoothed;
      DensDataResult[target].x_H2p_Smoothed = x_h2p_smoothed;
      DensDataResult[target].x_HD_Smoothed = x_hd_smoothed;
#endif
#endif

#if defined(BLACK_HOLES)
      DensDataResult[target].SmoothedEntr = smoothentr;
      DensDataResult[target].GasVel[0] = gasvel[0];
      DensDataResult[target].GasVel[1] = gasvel[1];
      DensDataResult[target].GasVel[2] = gasvel[2];
#endif
    }

  return 0;
}





int density_isactive(int n)
{
  if(P[n].TimeBin < 0)
    return 0;

#ifdef BLACK_HOLES
  if(P[n].Type == 5)
    return 1;
#endif
#ifdef BG_SFR
  if(P[n].Type == 4)
    return 1;
#endif

  if(P[n].Type == 0)
    return 1;

  return 0;
}
