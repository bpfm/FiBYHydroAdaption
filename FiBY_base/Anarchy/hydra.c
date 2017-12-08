#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>


#include "allvars.h"
#include "proto.h"
#include "kernel.h"

/* /\* DEBUG *\/ */
/* #include "bg_proto.h" */
/* /\* DEBUG *\/ */

/*! \file hydra.c
*  \brief Computation of SPH forces and rate of entropy generation
*
*  This file contains the "second SPH loop", where the SPH forces are
*  computed, and where the rate of change of entropy due to the shock heating
*  (via artificial viscosity) is computed.
*/

static double hubble_a, atime, hubble_a2, fac_mu, fac_vsic_fix, a3inv, fac_egy;


static struct hydrodata_in
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat Hsml;
  MyFloat Mass;
  MyFloat Density;
  MyFloat Pressure;
  MyFloat F1;
  MyFloat DhsmlDensityFactor;

  MyFloat ArtViscParam;
  MyFloat ArtDiffParam;

#ifdef PRESSURE_ENTROPY_SPH
  MyFloat EntropyPred;
  MyFloat EntropyVarPred;
  MyFloat WeightedDensity;
  MyFloat DhsmlPressure;
#endif

/*   /\* DEBUG *\/ */
/*   MyIDType ID; */
/*   double Temperature; */
/*   /\* DEBUG *\/ */
  int Timestep;
  int NodeList[NODELISTLENGTH];
}
 *HydroDataIn, *HydroDataGet;


static struct hydrodata_out
{
  MyLongDouble Acc[3];
  MyLongDouble DtEntropy;
/* #ifdef ALTERNATIVE_VISCOUS_TIMESTEP */
/*   MyFloat MinViscousDt; */
/* #else */
  MyFloat MaxSignalVel;
/* #endif */
}
 *HydroDataResult, *HydroDataOut;



/*! This function is the driver routine for the calculation of hydrodynamical
*  force and rate of change of entropy due to shock heating for all active
*  particles .
*/
void hydro_force(void)
{
  int i, j, k, ngrp, ndone, ndone_flag, dummy, normal_flag;

  int sendTask, recvTask, nexport, nimport, nimport_max, place;

  double soundspeed_i;

  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait, tstart, tend, t0, t1;


  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time)
#ifdef DARKENERGY
	+ DarkEnergy_a(All.Time);
#else
	+ All.OmegaLambda;
#endif
      hubble_a = All.Hubble * sqrt(hubble_a);
      hubble_a2 = All.Time * All.Time * hubble_a;

      fac_mu = pow(All.Time, 3 * (GAMMA - 1) / 2) / All.Time;

      fac_egy = pow(All.Time, 3 * (GAMMA - 1));

      fac_vsic_fix = hubble_a * pow(All.Time, 3 * GAMMA_MINUS1);

      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;
    }
  else
    hubble_a = hubble_a2 = atime = fac_mu = fac_vsic_fix = a3inv = fac_egy = 1.0;


  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct hydrodata_in) +
					     sizeof(struct hydrodata_out) +
					     sizemax(sizeof(struct hydrodata_in),
						     sizeof(struct hydrodata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));


  CPU_Step[CPU_HYDMISC] += measure_time();
  t0 = second();

  i = FirstActiveParticle;	/* first particle for this task */

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
	if(P[i].Type == 0)
	  {
	    if(hydro_evaluate(i, 0, &nexport, Send_count) < 0)
	      break;
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

      MPI_Allreduce(&nimport, &nimport_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

      if(nimport_max > 4 * All.BunchSize)
	{
	  /* we have a big imbalance between exports and imports on a certain
	     processor.  to prevent big imbalance in temporary memory
	     allocation, we process the pair-wise exchanges on the fly. But
	     note that this will be slower because it doesn't allow a
	     reduction of work-load imbalance losses by accumulating the total
	     work. */

	  if(ThisTask == 0)
	    {
	      printf("use special (slower) exchange in hydra.c to avoid large memory imbalance\n"
		     "we have BunchSize=%d and nimport_max/BunchSize=%g\n",
		     All.BunchSize, nimport_max / ((double) All.BunchSize));
	      fflush(stdout);
	    }
	  normal_flag = 0;
	}
      else
	normal_flag = 1;

      if(normal_flag)
	HydroDataGet =
	  (struct hydrodata_in *) mymalloc_msg(nimport * sizeof(struct hydrodata_in), "hydrodata_in");
      else
	HydroDataOut =
	  (struct hydrodata_out *) mymalloc_msg(nexport * sizeof(struct hydrodata_out), "hydrodata_out");

      HydroDataIn =
	(struct hydrodata_in *) mymalloc_msg(nexport * sizeof(struct hydrodata_in), "hydrodata_in");

      /* prepare particle data for export */

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    {
	      HydroDataIn[j].Pos[k] = P[place].Pos[k];
	      HydroDataIn[j].Vel[k] = SphP[place].VelPred[k];
	    }
	  HydroDataIn[j].Hsml = PPP[place].Hsml;
	  HydroDataIn[j].Mass = P[place].Mass;
	  HydroDataIn[j].DhsmlDensityFactor = SphP[place].h.DhsmlDensityFactor;
	  HydroDataIn[j].Density = SphP[place].d.Density;
	  HydroDataIn[j].Pressure = SphP[place].Pressure;
	  HydroDataIn[j].Timestep = (P[place].TimeBin ? (1 << P[place].TimeBin) : 0);
	  
	  //HydroDataIn[j].ArtViscParam = SphP[place].ArtViscParam;
	  //HydroDataIn[j].ArtDiffParam = SphP[place].ArtDiffParam;

#ifdef PRESSURE_ENTROPY_SPH
	  HydroDataIn[j].EntropyPred = SphP[place].EntropyPred;
	  HydroDataIn[j].EntropyVarPred = SphP[place].EntropyVarPred;
	  HydroDataIn[j].WeightedDensity = SphP[place].cky.WeightedDensity;
	  HydroDataIn[j].DhsmlPressure = SphP[place].dky.DhsmlPressure;
#endif

/* 	  /\* DEBUG *\/ */
/* 	  HydroDataIn[j].ID = P[place].ID; */
/* 	  HydroDataIn[j].Temperature = bg_get_temperature(place); */
/* 	  /\* DEBUG *\/ */

	  /* calculation of F1 */
#ifndef ALTVISCOSITY
	  soundspeed_i = sqrt(GAMMA * SphP[place].Pressure / SphP[place].d.Density);
	  HydroDataIn[j].F1 = fabs(SphP[place].v.DivVel) /
	    (fabs(SphP[place].v.DivVel) + SphP[place].r.CurlVel +
	     0.0001 * soundspeed_i / PPP[place].Hsml / fac_mu);
#else
	  HydroDataIn[j].F1 = SphP[place].v.DivVel;
#endif

	  memcpy(HydroDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));

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
		  if(!normal_flag)
		    {
		      HydroDataGet =
			(struct hydrodata_in *) mymalloc_msg(Recv_count[recvTask] *
							     sizeof(struct hydrodata_in), "hyd_get");
		      HydroDataGet -= Recv_offset[recvTask];
		      HydroDataResult =
			(struct hydrodata_out *) mymalloc_msg(Recv_count[recvTask] *
							      sizeof(struct hydrodata_out), "hyd_out");
		      HydroDataResult -= Recv_offset[recvTask];
		    }

		  /* get the particles */
		  MPI_Sendrecv(&HydroDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A,
			       &HydroDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		  if(!normal_flag)
		    {
		      double t0, t1;

		      t0 = second();

		      for(j = Recv_offset[recvTask]; j < Recv_offset[recvTask] + Recv_count[recvTask]; j++)
			hydro_evaluate(j, 1, &dummy, &dummy);

		      t1 = second();
		      timecomp2 += timediff(t0, t1);
		      timecommsumm1 -= timediff(t0, t1);

		      MPI_Sendrecv(&HydroDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B,
				   &HydroDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		      myfree(HydroDataResult + Recv_offset[recvTask]);
		      myfree(HydroDataGet + Recv_offset[recvTask]);
		    }
		}
	    }
	}
      tend = second();
      timecommsumm1 += timediff(tstart, tend);

      myfree(HydroDataIn);

      if(normal_flag)
	{
	  HydroDataResult =
	    (struct hydrodata_out *) mymalloc_msg(nimport * sizeof(struct hydrodata_out), "hyd_result");
	  HydroDataOut =
	    (struct hydrodata_out *) mymalloc_msg(nexport * sizeof(struct hydrodata_out), "hyd_out");

	  /* now do the particles that were sent to us */

	  tstart = second();
	  for(j = 0; j < nimport; j++)
	    {
	      hydro_evaluate(j, 1, &dummy, &dummy);
	    }
	  tend = second();
	  timecomp2 += timediff(tstart, tend);
	}

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      tstart = second();
      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      tend = second();
      timewait2 += timediff(tstart, tend);

      if(normal_flag)
	{
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
		      MPI_Sendrecv(&HydroDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B,
				   &HydroDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }
	  tend = second();
	  timecommsumm2 += timediff(tstart, tend);
	}


      /* add the result to the local particles */
      tstart = second();
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    {
	      SphP[place].a.dHydroAccel[k] += HydroDataOut[j].Acc[k];
	    }

	  SphP[place].e.dDtEntropy += HydroDataOut[j].DtEntropy;

/* #ifdef ALTERNATIVE_VISCOUS_TIMESTEP */
/* 	  if(SphP[place].MinViscousDt > HydroDataOut[j].MinViscousDt) */
/* 	    SphP[place].MinViscousDt = HydroDataOut[j].MinViscousDt; */
/* #else */
	  if(SphP[place].MaxSignalVel < HydroDataOut[j].MaxSignalVel)
	    SphP[place].MaxSignalVel = HydroDataOut[j].MaxSignalVel;
/* #endif */

	}
      tend = second();
      timecomp1 += timediff(tstart, tend);

      if(normal_flag)
	{
	  myfree(HydroDataOut);
	  myfree(HydroDataResult);
	  myfree(HydroDataGet);
	}
      else
	myfree(HydroDataOut);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);

  myfree(Ngblist);


  /* do final operations on results */


#ifdef FLTROUNDOFFREDUCTION
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
	SphP[i].e.DtEntropy = FLT(SphP[i].e.dDtEntropy);

	for(j = 0; j < 3; j++)
	  SphP[i].a.HydroAccel[j] = FLT(SphP[i].a.dHydroAccel[j]);
      }
#endif



  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
	/* Translate energy change rate into entropy change rate */
#ifndef PRESSURE_ENTROPY_SPH
	SphP[i].e.DtEntropy *= GAMMA_MINUS1 / (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1));
#else
	SphP[i].e.DtEntropy *= GAMMA_MINUS1 / (hubble_a2 * pow(SphP[i].cky.WeightedDensity, GAMMA_MINUS1));
#endif

#ifdef SPH_BND_PARTICLES
	if(P[i].ID == 0)
	  {
	    SphP[i].e.DtEntropy = 0;
	    for(k = 0; k < 3; k++)
	      SphP[i].a.HydroAccel[k] = 0;
	  }
#endif
      }

  /* collect some timing information */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_HYDCOMPUTE] += timecomp;
  CPU_Step[CPU_HYDWAIT] += timewait;
  CPU_Step[CPU_HYDCOMM] += timecomm;
  CPU_Step[CPU_HYDMISC] += timeall - (timecomp + timewait + timecomm);
}




/*! This function is the 'core' of the SPH force computation. A target
*  particle is specified which may either be local, or reside in the
*  communication buffer.
*/
int hydro_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int startnode, numngb, listindex = 0;

  int j, k, n, timestep;

  MyDouble *pos;

  MyFloat *vel;

  MyFloat mass, h_i, dhsmlDensityFactor, rho, pressure, f1, f2;

  MyLongDouble acc[3], dtEntropy;

/* #ifdef ALTERNATIVE_VISCOUS_TIMESTEP */
/*   MyFloat minViscousDt; */
/* #else */
  MyFloat maxSignalVel;
/* #endif */
  double dx, dy, dz, dvx, dvy, dvz;

  double h_i2, hinv, hinv4;

  double p_over_rho2_i, p_over_rho2_j, soundspeed_i, soundspeed_j;

  double hfc, dwk_i, vdotr, vdotr2, visc, mu_ij, rho_ij, vsig;

  double h_j, dwk_j;

  double r, r2, u;

  double hfc_visc;

  double dmin1, dmin2;

  double BulkVisc_ij;

  int imax1, imax2;

#ifdef ALTVISCOSITY
  double mu_i, mu_j;
#endif

#ifndef NOVISCOSITYLIMITER
  double dt;
#endif

#ifdef CONVENTIONAL_VISCOSITY
  double c_ij, h_ij;
#endif

#ifdef PRESSURE_ENTROPY_SPH
  double entropyPred, entropyVarPred, weightedDensity, dhsmlPressure, f_i, f_j;
#endif

  double u_i, u_j;

/*   /\* DEBUG *\/ */
/*   MyIDType id; */
/*   double temperature; */
/*   /\* DEBUG *\/ */

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h_i = PPP[target].Hsml;
      mass = P[target].Mass;
      dhsmlDensityFactor = SphP[target].h.DhsmlDensityFactor;
      rho = SphP[target].d.Density;
      pressure = SphP[target].Pressure;
      timestep = TISTEP(P[target].TimeBin);

#ifdef PRESSURE_ENTROPY_SPH
      entropyPred = SphP[target].EntropyPred;
      entropyVarPred = SphP[target].EntropyVarPred;
      weightedDensity = SphP[target].cky.WeightedDensity;
      dhsmlPressure = SphP[target].dky.DhsmlPressure;
#endif

#ifndef PRESSURE_ENTROPY_SPH
      p_over_rho2_i = pressure / (rho * rho);
      soundspeed_i = sqrt(GAMMA * pressure / rho);
#else
      p_over_rho2_i = pressure / (SphP[target].cky.WeightedDensity * SphP[target].cky.WeightedDensity);
      soundspeed_i = sqrt(GAMMA * pressure / SphP[target].cky.WeightedDensity);
      u_i = p_over_rho2_i * SphP[target].cky.WeightedDensity / GAMMA_MINUS1;
#endif

/*       /\* DEBUG *\/ */
/*       id = P[target].ID; */
/*       temperature = bg_get_temperature(target); */
/*       /\* DEBUG *\/ */

#ifndef ALTVISCOSITY
      f1 = fabs(SphP[target].v.DivVel) /
	(fabs(SphP[target].v.DivVel) + SphP[target].r.CurlVel +
	 0.0001 * soundspeed_i / PPP[target].Hsml / fac_mu);
#else
      f1 = SphP[target].v.DivVel;
#endif
    }
  else
    {
      pos = HydroDataGet[target].Pos;
      vel = HydroDataGet[target].Vel;
      h_i = HydroDataGet[target].Hsml;
      mass = HydroDataGet[target].Mass;
      dhsmlDensityFactor = HydroDataGet[target].DhsmlDensityFactor;
      rho = HydroDataGet[target].Density;
      pressure = HydroDataGet[target].Pressure;
      timestep = HydroDataGet[target].Timestep;

#ifdef PRESSURE_ENTROPY_SPH
      entropyPred = HydroDataGet[target].EntropyPred;
      entropyVarPred = HydroDataGet[target].EntropyVarPred;
      weightedDensity = HydroDataGet[target].WeightedDensity;
      dhsmlPressure = HydroDataGet[target].DhsmlPressure;
#endif

#ifndef PRESSURE_ENTROPY_SPH
      p_over_rho2_i = pressure / (rho * rho);
      soundspeed_i = sqrt(GAMMA * pressure / rho);
#else
      p_over_rho2_i = pressure / (weightedDensity * weightedDensity); 
      soundspeed_i = sqrt(GAMMA * pressure / weightedDensity);
      u_i = p_over_rho2_i * weightedDensity / GAMMA_MINUS1;
#endif

      f1 = HydroDataGet[target].F1;

/*       /\* DEBUG *\/ */
/*       id = HydroDataGet[target].ID; */
/*       temperature = HydroDataGet[target].Temperature; */
/*       /\* DEBUG *\/ */
    }



  /* initialize variables before SPH loop is started */

#ifndef PRESSURE_ENTROPY_SPH
  p_over_rho2_i *= dhsmlDensityFactor;
#endif

  h_i2 = h_i * h_i;



  acc[0] = acc[1] = acc[2] = dtEntropy = 0;

/* #ifdef ALTERNATIVE_VISCOUS_TIMESTEP */
/*   minViscousDt = 1.0e32; */
/* #else */
  maxSignalVel = 0;
/* #endif */


  /* Now start the actual SPH computation for this particle */

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = HydroDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb = ngb_treefind_pairs(pos, h_i, target, &startnode, mode, nexport, nsend_local);

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
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
	      h_j = PPP[j].Hsml;
	      if(r2 < h_i2 || r2 < h_j * h_j)
		{
		  r = sqrt(r2);
		  if(r > 0)
		    {
#ifndef PRESSURE_ENTROPY_SPH
		      p_over_rho2_j = SphP[j].Pressure / (SphP[j].d.Density * SphP[j].d.Density);
		      soundspeed_j = sqrt(GAMMA * p_over_rho2_j * SphP[j].d.Density);
#else
		      p_over_rho2_j = SphP[j].Pressure / (SphP[j].cky.WeightedDensity * SphP[j].cky.WeightedDensity);
		      soundspeed_j = sqrt(GAMMA * p_over_rho2_j * SphP[j].cky.WeightedDensity);
		      u_j = p_over_rho2_j * SphP[j].cky.WeightedDensity / GAMMA_MINUS1;
#endif
		      dvx = vel[0] - SphP[j].VelPred[0];
		      dvy = vel[1] - SphP[j].VelPred[1];
		      dvz = vel[2] - SphP[j].VelPred[2];
		      vdotr = dx * dvx + dy * dvy + dz * dvz;

		      if(All.ComovingIntegrationOn)
			vdotr2 = vdotr + hubble_a2 * r2;
		      else
			vdotr2 = vdotr;

		      if(r2 < h_i2)
			{
			  hinv = 1.0 / h_i;
#ifndef  TWODIMS
			  hinv4 = hinv * hinv * hinv * hinv;
#else
			  hinv4 = hinv * hinv * hinv; // / boxSize_Z;
#endif
			  u = r * hinv;

/* 			  if(u < 0.5) */
/* 			    dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4); */
/* 			  else */
/* 			    dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u); */

			  dwk_i = kernel_dwk(u, hinv4);
			}
		      else
			{
			  dwk_i = 0;
			}

		      if(r2 < h_j * h_j)
			{
			  hinv = 1.0 / h_j;
#ifndef  TWODIMS
			  hinv4 = hinv * hinv * hinv * hinv;
#else
			  hinv4 = hinv * hinv * hinv; // / boxSize_Z;
#endif
			  u = r * hinv;

/* 			  if(u < 0.5) */
/* 			    dwk_j = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4); */
/* 			  else */
/* 			    dwk_j = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u); */

			  dwk_j = kernel_dwk(u, hinv4);
			}
		      else
			{
			  dwk_j = 0;
			}

		      vsig = soundspeed_i + soundspeed_j;

/* #ifndef ALTERNATIVE_VISCOUS_TIMESTEP */
		      if(vsig > maxSignalVel)
			maxSignalVel = vsig;
/* #endif */
		      if(vdotr2 < 0)	/* ... artificial viscosity */
			{
#ifndef ALTVISCOSITY
#ifndef CONVENTIONAL_VISCOSITY
			  mu_ij = fac_mu * vdotr2 / r;	/* note: this is negative! */
#else
			  c_ij = 0.5 * (soundspeed_i + soundspeed_j);
			  h_ij = 0.5 * (h_i + h_j);
			  mu_ij = fac_mu * h_ij * vdotr2 / (r2 + 0.0001 * h_ij * h_ij);
#endif
			  vsig -= 3 * mu_ij;

/* #ifndef ALTERNATIVE_VISCOUS_TIMESTEP */
			  if(vsig > maxSignalVel)
			    maxSignalVel = vsig;
/* #endif */

			  rho_ij = 0.5 * (rho + SphP[j].d.Density);
			  f2 =
			    fabs(SphP[j].v.DivVel) / (fabs(SphP[j].v.DivVel) + SphP[j].r.CurlVel +
						      0.0001 * soundspeed_j / fac_mu / PPP[j].Hsml);
#ifdef NO_SHEAR_VISCOSITY_LIMITER
			  f1 = f2 = 1;
#endif

			  BulkVisc_ij = All.ArtBulkViscConst;

#ifndef CONVENTIONAL_VISCOSITY
			  visc = 0.25 * BulkVisc_ij * vsig * (-mu_ij) / rho_ij * (f1 + f2);
			  //visc = 0.5 * BulkVisc_ij * vsig * (-mu_ij) / rho_ij;
#else
			  visc =
			    (-BulkVisc_ij * mu_ij * c_ij + 2 * BulkVisc_ij * mu_ij * mu_ij) /
			    rho_ij * (f1 + f2) * 0.5;
#endif

#else /* start of ALTVISCOSITY block */
			  if(f1 < 0)
			    mu_i = h_i * fabs(f1);	/* f1 hold here the velocity divergence of particle i */
			  else
			    mu_i = 0;
			  if(SphP[j].u.s.a4.DivVel < 0)
			    mu_j = h_j * fabs(SphP[j].u.s.a4.DivVel);
			  else
			    mu_j = 0;
			  visc = All.ArtBulkViscConst * ((soundspeed_i + mu_i) * mu_i / rho +
							 (soundspeed_j + mu_j) * mu_j / SphP[j].d.Density);
#endif /* end of ALTVISCOSITY block */


			  /* .... end artificial viscosity evaluation */
			  /* now make sure that viscous acceleration is not too large */
/* #ifdef ALTERNATIVE_VISCOUS_TIMESTEP */
/* 			  if(visc > 0) */
/* 			    { */
/* 			      dt = fac_vsic_fix * vdotr2 / */
/* 				(0.5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * visc); */

/* 			      dt /= hubble_a; */

/* 			      if(dt < minViscousDt) */
/* 				minViscousDt = dt; */
/* 			    } */
/* #endif */

#ifndef NOVISCOSITYLIMITER
			  dt = 2 * IMAX(timestep, TISTEP(P[j].TimeBin)) * All.Timebase_interval;
			  if(dt > 0 && (dwk_i + dwk_j) < 0)
			    {
#ifdef BLACK_HOLES
			      if((mass + P[j].Mass) > 0)
#endif
				visc = DMIN(visc, 0.5 * fac_vsic_fix * vdotr2 /
					    (0.5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * dt));
			    }
#endif
			}
		      else
			{
			  visc = 0;
			}

#ifndef PRESSURE_ENTROPY_SPH
		      p_over_rho2_j *= SphP[j].h.DhsmlDensityFactor;
#endif

		      hfc_visc = 0.5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;

		      //#ifndef RPSPH
#ifdef PRESSURE_ENTROPY_SPH
		      f_i = dhsmlPressure * dhsmlDensityFactor;
                      f_j = SphP[j].dky.DhsmlPressure * SphP[j].h.DhsmlDensityFactor;

                      hfc = hfc_visc + P[j].Mass *
                        ((SphP[j].EntropyVarPred / entropyVarPred - f_i) * dwk_i * p_over_rho2_i +
                         (entropyVarPred / SphP[j].EntropyVarPred - f_j) * dwk_j * p_over_rho2_j) / r;
#else
		      /* Formulation derived from the Lagrangian */
		      hfc = hfc_visc + P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;
#endif
		      //#else
		      /* hfc = hfc_visc + P[j].Mass/SphP[j].Density * (SphP[j].Pressure-pressure)/SphP[j].Density*(dwk_j+dwk_i)/r/2; */

		      //		      hfc = hfc_visc + P[j].Mass * (SphP[j].Pressure - pressure) * (dwk_j + dwk_i) /
		      //			(2 * r * SphP[j].d.Density * SphP[j].d.Density);
		      //#endif

#ifndef NOACCEL
		      acc[0] += FLT(-hfc * dx);
		      acc[1] += FLT(-hfc * dy);
		      acc[2] += FLT(-hfc * dz);
#endif

		      dtEntropy += FLT(0.5 * hfc_visc * vdotr2);

/* 		      /\* DEBUG *\/ */
/* 		      if(temperature > 5e7) */
/* 			{ */
/* 			  printf("[HYDRA %2d] PID=%d, NID=%d, hfc=%g, hfc_visc=%g, vdotr2=%g, r=%g, p_over_rho2_i=%g, p_over_rho2_j=%g\n", */
/* 				 ThisTask, id, P[j].ID, hfc, hfc_visc, vdotr2, r, p_over_rho2_i, p_over_rho2_j); */
/* 			  printf("[HYDRA %2d] PID=%d, NID=%d, BulkVisc_ij=%g, vsig=%g, mu_ij=%g, rho_ij=%g, f1=%g, f2=%g\n", */
/* 				 ThisTask, id, P[j].ID, BulkVisc_ij, vsig, mu_ij, rho_ij, f1, f2); */
/* 			  printf("[HYDRA %2d] PID=%d, NID=%d, p_i=%g, p_j=%g, c_i=%g, c_j=%g, dwk_i=%g, dwk_j=%g\n", */
/* 				 ThisTask, id, P[j].ID, pressure, SphP[j].Pressure, soundspeed_i, soundspeed_j, dwk_i, dwk_j); */
/* 			  printf("[HYDRA %2d] PID=%d, NID=%d, h_i=%g, h_j=%g, rho_i=%g, rho_j=%g\n", */
/* 				 ThisTask, id, P[j].ID, h_i, h_j, rho, SphP[j].d.Density); */
/* 			  printf("[HYDRA %2d] PID=%d, NID=%d, pos_i=(%g,%g,%g), pos_j=(%g,%g,%g)\n", */
/* 				 ThisTask, id, P[j].ID, pos[0], pos[1], pos[2], P[j].Pos[0], P[j].Pos[1], P[j].Pos[2]); */
/* 			  printf("[HYDRA %2d] PID=%d, NID=%d, vel_i=(%g,%g,%g), vel_j=(%g,%g,%g)\n", */
/* 				 ThisTask, id, P[j].ID, vel[0], vel[1], vel[2], P[j].Vel[0], P[j].Vel[1], P[j].Vel[2]); */
/* 			} */
/* 		      /\* DEBUG *\/ */
		    }
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = HydroDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }


  /* Now collect the result at the right place */
  if(mode == 0)
    {
      for(k = 0; k < 3; k++)
	SphP[target].a.dHydroAccel[k] = acc[k];
      SphP[target].e.dDtEntropy = dtEntropy;
/* #ifdef ALTERNATIVE_VISCOUS_TIMESTEP */
/*       SphP[target].MinViscousDt = minViscousDt; */
/* #else */
      SphP[target].MaxSignalVel = maxSignalVel;
/* #endif */
    }
  else
    {
      for(k = 0; k < 3; k++)
	HydroDataResult[target].Acc[k] = acc[k];
      HydroDataResult[target].DtEntropy = dtEntropy;
/* #ifdef ALTERNATIVE_VISCOUS_TIMESTEP */
/*       HydroDataResult[target].MinViscousDt = minViscousDt; */
/* #else */
      HydroDataResult[target].MaxSignalVel = maxSignalVel;
/* #endif */
    }

  return 0;
}
