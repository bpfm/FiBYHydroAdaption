#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

/*! \file gravtree.c 
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for
 *  all active local particles, and particles are exported to other
 *  processors if needed, where they can receive additional force
 *  contributions. If the TreePM algorithm is enabled, the force computed
 *  will only be the short-range part.
 */



/*! This function computes the gravitational forces for all active
 *  particles.  If needed, a new tree is constructed, otherwise the
 *  dynamically updated tree is used.  Particles are only exported to other
 *  processors when really needed, thereby allowing a good use of the
 *  communication buffer.
 */
void gravity_tree(void)
{
  long long n_exported = 0, n_nodesinlist = 0;

  int i, j, iter = 0, ewald_iter, ewald_max, ret;

  int maxnumnodes, nexport, nimport, nodesinlist;

  double tstart, tend, t0, t1;

  double timeall = 0, timetree1 = 0, timetree2 = 0;

  double timetree, timewait, timecomm;

  double timecommsumm1 = 0, timecommsumm2 = 0, timewait1 = 0, timewait2 = 0;

  double ewaldcount, costtotal, sum_costtotal, ewaldtot;

  double maxt, sumt, maxt1, sumt1, maxt2, sumt2, sumcommall, sumwaitall;

  double plb, plb_max;

#ifndef NOGRAVITY
  int ndone, ndone_flag, ngrp;

  int k, place, dummy;

  int sendTask, recvTask;

  double ax, ay, az;

  MPI_Status status;
#endif

  CPU_Step[CPU_MISC] += measure_time();

  /* set new softening lengths */
  if(All.ComovingIntegrationOn)
    set_softenings();

  /* contruct tree if needed */

  if(TreeReconstructFlag)
    {
      if(ThisTask == 0)
	printf("Tree construction.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

      CPU_Step[CPU_MISC] += measure_time();

      force_treebuild(NumPart, NULL);

      TreeReconstructFlag = 0;

      if(ThisTask == 0)
	printf("Tree construction done.\n");
    }

  CPU_Step[CPU_TREEMISC] += measure_time();
  t0 = second();

#ifndef NOGRAVITY
  if(ThisTask == 0)
    printf("Begin tree force.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct gravdata_in) + sizeof(struct gravdata_out) +
					     sizemax(sizeof(struct gravdata_in),
						     sizeof(struct gravdata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));

  if(ThisTask == 0)
    printf("All.BunchSize=%d\n", All.BunchSize);

  costtotal = ewaldcount = 0;

#if defined(PERIODIC) && !defined(PMGRID)
  ewald_max = 1;
#else
  ewald_max = 0;
#endif

  for(ewald_iter = 0; ewald_iter <= ewald_max; ewald_iter++)
    {

      i = FirstActiveParticle;	/* beginn with this index */

      do
	{
	  iter++;

	  for(j = 0; j < NTask; j++)
	    {
	      Send_count[j] = 0;
	      Exportflag[j] = -1;
	    }

	  /* do local particles and prepare export list */
	  tstart = second();
	  for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	    {
#if !defined(PMGRID)
#if defined(PERIODIC)
	      if(ewald_iter)
		{
		  ret = force_treeevaluate_ewald_correction(i, 0, &nexport, Send_count);
		  if(ret >= 0)
		    ewaldcount += ret;
		  else
		    break;	/* export buffer has filled up */
		}
	      else
#endif
		{
		  ret = force_treeevaluate(i, 0, &nexport, Send_count);
		  if(ret >= 0)
		    {
		      costtotal += ret;
		    }
		  else
		    break;	/* export buffer has filled up */
		}
#else
	      ret = force_treeevaluate_shortrange(i, 0, &nexport, Send_count);
	      if(ret >= 0)
		{
		  costtotal += ret;
		}
	      else
		break;		/* export buffer has filled up */
#endif
	    }
	  tend = second();
	  timetree1 += timediff(tstart, tend);

	  n_exported += nexport;

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

	  GravDataGet = (struct gravdata_in *) mymalloc(nimport * sizeof(struct gravdata_in));
	  GravDataIn = (struct gravdata_in *) mymalloc(nexport * sizeof(struct gravdata_in));

	  /* prepare particle data for export */

	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      for(k = 0; k < 3; k++)
		GravDataIn[j].Pos[k] = P[place].Pos[k];

#ifdef UNEQUALSOFTENINGS
	      GravDataIn[j].Type = P[place].Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	      if(P[place].Type == 0)
		GravDataIn[j].Soft = SphP[place].Hsml;
#endif
#endif
	      GravDataIn[j].OldAcc = P[place].OldAcc;

	      memcpy(GravDataIn[j].NodeList,
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
		      /* get the particles */
		      MPI_Sendrecv(&GravDataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A,
				   &GravDataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);
		    }
		}
	    }
	  tend = second();
	  timecommsumm1 += timediff(tstart, tend);


	  myfree(GravDataIn);
	  GravDataResult = (struct gravdata_out *) mymalloc(nimport * sizeof(struct gravdata_out));
	  GravDataOut = (struct gravdata_out *) mymalloc(nexport * sizeof(struct gravdata_out));


	  /* now do the particles that were sent to us */
	  tstart = second();
	  for(j = 0; j < nimport; j++)
	    {
#if !defined(PMGRID)
#if defined(PERIODIC)
	      if(ewald_iter)
		ewaldcount += force_treeevaluate_ewald_correction(j, 1, &dummy, &dummy);
	      else
#endif
		{
		  costtotal += force_treeevaluate(j, 1, &nodesinlist, &dummy);
		  n_nodesinlist += nodesinlist;
		}
#else
	      costtotal += force_treeevaluate_shortrange(j, 1, &nodesinlist, &dummy);
#endif
	    }
	  tend = second();
	  timetree2 += timediff(tstart, tend);

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
		      MPI_Sendrecv(&GravDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct gravdata_out),
				   MPI_BYTE, recvTask, TAG_GRAV_B,
				   &GravDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct gravdata_out),
				   MPI_BYTE, recvTask, TAG_GRAV_B, MPI_COMM_WORLD, &status);
		    }
		}

	    }
	  tend = second();
	  timecommsumm2 += timediff(tstart, tend);


	  /* add the results to the local particles */
	  tstart = second();
	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      for(k = 0; k < 3; k++)
		P[place].g.dGravAccel[k] += GravDataOut[j].Acc[k];

	      P[place].GravCost += GravDataOut[j].Ninteractions;
#ifdef EVALPOTENTIAL
	      P[place].p.dPotential += GravDataOut[j].Potential;
#endif
	    }
	  tend = second();
	  timetree1 += timediff(tstart, tend);

	  myfree(GravDataOut);
	  myfree(GravDataResult);
	  myfree(GravDataGet);
	}
      while(ndone < NTask);
    }

  myfree(DataNodeList);
  myfree(DataIndexTable);



#ifdef FLTROUNDOFFREDUCTION
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#ifdef EVALPOTENTIAL
      P[i].p.Potential = FLT(P[i].p.dPotential);
#endif
      for(j = 0; j < 3; j++)
	P[i].g.GravAccel[j] = FLT(P[i].g.dGravAccel[j]);
    }
#endif

  /* now add things for comoving integration */

#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn)
    {
      double fac = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;

      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	for(j = 0; j < 3; j++)
	  P[i].g.GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#ifdef PMGRID
      ax = P[i].g.GravAccel[0] + P[i].GravPM[0] / All.G;
      ay = P[i].g.GravAccel[1] + P[i].GravPM[1] / All.G;
      az = P[i].g.GravAccel[2] + P[i].GravPM[2] / All.G;
#else
      ax = P[i].g.GravAccel[0];
      ay = P[i].g.GravAccel[1];
      az = P[i].g.GravAccel[2];
#endif
      P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);
    }


  if(All.TypeOfOpeningCriterion == 1)
    All.ErrTolTheta = 0;	/* This will switch to the relative opening criterion for the following force computations */

  /*  muliply by G */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      for(j = 0; j < 3; j++)
	P[i].g.GravAccel[j] *= All.G;

#ifdef EVALPOTENTIAL
      /* remove self-potential */
      P[i].p.Potential += P[i].Mass / All.SofteningTable[P[i].Type];

      if(All.ComovingIntegrationOn)
	if(All.PeriodicBoundariesOn)
	  P[i].p.Potential -= 2.8372975 * pow(P[i].Mass, 2.0 / 3) *
	    pow(All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G), 1.0 / 3);

      P[i].p.Potential *= All.G;

#ifdef PMGRID
      P[i].p.Potential += P[i].PM_Potential;	/* add in long-range potential */
#endif

      if(All.ComovingIntegrationOn)
	{
#ifndef PERIODIC
	  double fac, r2;

	  fac = -0.5 * All.Omega0 * All.Hubble * All.Hubble;

	  for(k = 0, r2 = 0; k < 3; k++)
	    r2 += P[i].Pos[k] * P[i].Pos[k];

	  P[i].p.Potential += fac * r2;
#endif
	}
      else
	{
	  double fac, r2;

	  fac = -0.5 * All.OmegaLambda * All.Hubble * All.Hubble;

	  if(fac != 0)
	    {
	      for(k = 0, r2 = 0; k < 3; k++)
		r2 += P[i].Pos[k] * P[i].Pos[k];

	      P[i].p.Potential += fac * r2;
	    }
	}
#endif
    }

  /* Finally, the following factor allows a computation of a cosmological simulation 
     with vacuum energy in physical coordinates */
#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn == 0)
    {
      double fac = All.OmegaLambda * All.Hubble * All.Hubble;

      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	for(j = 0; j < 3; j++)
	  P[i].g.GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif


  if(ThisTask == 0)
    printf("tree is done.\n");

#else /* gravity is switched off */

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    for(j = 0; j < 3; j++)
      P[i].g.GravAccel[j] = 0;

#endif /* end of NOGRAVITY */


#ifdef STATICNFW
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      r = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]);
      m = enclosed_mass(r);
#ifdef NFW_DARKFRACTION
      m *= NFW_DARKFRACTION;
#endif
      if(r > 0)
	{
	  for(k = 0; k < 3; k++)
	    P[i].g.GravAccel[k] += -All.G * m * P[i].Pos[k] / (r * r * r);
	}
    }
#endif



#ifdef STATICHQ
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      r = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]);

      a = pow(All.G * HQ_M200 / (100 * All.Hubble * All.Hubble), 1.0 / 3) / HQ_C *
	sqrt(2 * (log(1 + HQ_C) - HQ_C / (1 + HQ_C)));

      m = HQ_M200 * pow(r / (r + a), 2);

#ifdef HQ_DARKFRACTION
      m *= HQ_DARKFRACTION;
#endif
      if(r > 0)
	{
	  for(k = 0; k < 3; k++)
	    P[i].g.GravAccel[k] += -All.G * m * P[i].Pos[k] / (r * r * r);
	}
    }
#endif


  /* Now the force computation is finished */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  /*  gather some diagnostic information */

  timetree = timetree1 + timetree2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  MPI_Reduce(&timetree, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timetree, &maxt, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timetree1, &sumt1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timetree1, &maxt1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timetree2, &sumt2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timetree2, &maxt2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timewait, &sumwaitall, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecomm, &sumcommall, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&costtotal, &sum_costtotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&ewaldcount, &ewaldtot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  sumup_longs(1, &n_exported, &n_exported);
  sumup_longs(1, &n_nodesinlist, &n_nodesinlist);

  All.TotNumOfForces += GlobNumForceUpdate;

  plb = (NumPart / ((double) All.TotNumPart)) * NTask;
  MPI_Reduce(&plb, &plb_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Numnodestree, &maxnumnodes, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  /*
  if(ThisTask == 0)
    printf("timeall=%g, timetree=%g, timewait=%g, timecomm=%g, timetree1=%g, timetree2=%g, timecommsumm1=%g, timecommsumm2=%g, timewait1=%g, timewait2=%g\n",
	   timeall, timetree, timewait, timecomm, timetree1, timetree2, timecommsumm1, timecommsumm2, timewait1, timewait2);
  */

  CPU_Step[CPU_TREEMISC] += timeall - (timetree + timewait + timecomm);
  CPU_Step[CPU_TREEWALK1] += timetree1;
  CPU_Step[CPU_TREEWALK2] += timetree2;
  CPU_Step[CPU_TREESEND] += timecommsumm1;
  CPU_Step[CPU_TREERECV] += timecommsumm2;
  CPU_Step[CPU_TREEWAIT1] += timewait1;
  CPU_Step[CPU_TREEWAIT2] += timewait2;


  if(ThisTask == 0)
    {
      fprintf(FdTimings, "Step= %d  t= %g  dt= %g \n", All.NumCurrentTiStep, All.Time, All.TimeStep);
      fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g (%g) iter= %d\n",
	      (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000),
	      (int) (All.TotNumOfForces / 1000000000), (int) (All.TotNumOfForces % 1000000000),
	      n_exported / ((double) GlobNumForceUpdate), n_nodesinlist / ((double) n_exported + 1.0e-10),
	      iter);
      /* note: on Linux, the 8-byte integer could be printed with the format identifier "%qd", but doesn't work on AIX */

      fprintf(FdTimings, "work-load balance: %g (%g %g) rel1to2=%g   max=%g avg=%g\n",
	      maxt / (sumt / NTask), maxt1 / (sumt1 / NTask), maxt2 / (sumt2 / NTask),
	      sumt1 / (sumt1 + sumt2), maxt, sumt / NTask);
      fprintf(FdTimings, "particle-load balance: %g\n", plb_max);
      fprintf(FdTimings, "max. nodes: %d, filled: %g\n", maxnumnodes,
	      maxnumnodes / (All.TreeAllocFactor * All.MaxPart + NTopnodes));
      fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g (%g)\n", GlobNumForceUpdate / (sumt + 1.0e-20),
	      GlobNumForceUpdate / (maxt * NTask), ((double) (sum_costtotal)) / GlobNumForceUpdate,
	      ((double) ewaldtot) / GlobNumForceUpdate);
      fprintf(FdTimings, "\n");

      fflush(FdTimings);

/* #ifdef FORCE_DOMAIN_DECOMPOSITION */
/*       PartsPerSecond = GlobNumForceUpdate / (sumt + 1.0e-20); */
/*       printf("Actual part/sec=%f\n", PartsPerSecond); */
/*       fflush(stdout); */
/* #endif */
    }

/* #ifdef FORCE_DOMAIN_DECOMPOSITION */
/*   MPI_Bcast(&PartsPerSecond, sizeof(PartsPerSecond), MPI_BYTE, 0, MPI_COMM_WORLD); */
/* #endif */

  CPU_Step[CPU_TREEMISC] += measure_time();
}


/*! This function sets the (comoving) softening length of all particle
 *  types in the table All.SofteningTable[...].  We check that the physical
 *  softening length is bounded by the Softening-MaxPhys values.
 */
void set_softenings(void)
{
  int i;

  if(All.ComovingIntegrationOn)
    {
      if(All.SofteningGas * All.Time > All.SofteningGasMaxPhys)
	All.SofteningTable[0] = All.SofteningGasMaxPhys / All.Time;
      else
	All.SofteningTable[0] = All.SofteningGas;

      if(All.SofteningHalo * All.Time > All.SofteningHaloMaxPhys)
	All.SofteningTable[1] = All.SofteningHaloMaxPhys / All.Time;
      else
	All.SofteningTable[1] = All.SofteningHalo;

      if(All.SofteningDisk * All.Time > All.SofteningDiskMaxPhys)
	All.SofteningTable[2] = All.SofteningDiskMaxPhys / All.Time;
      else
	All.SofteningTable[2] = All.SofteningDisk;

      if(All.SofteningBulge * All.Time > All.SofteningBulgeMaxPhys)
	All.SofteningTable[3] = All.SofteningBulgeMaxPhys / All.Time;
      else
	All.SofteningTable[3] = All.SofteningBulge;

      if(All.SofteningStars * All.Time > All.SofteningStarsMaxPhys)
	All.SofteningTable[4] = All.SofteningStarsMaxPhys / All.Time;
      else
	All.SofteningTable[4] = All.SofteningStars;

      if(All.SofteningBndry * All.Time > All.SofteningBndryMaxPhys)
	All.SofteningTable[5] = All.SofteningBndryMaxPhys / All.Time;
      else
	All.SofteningTable[5] = All.SofteningBndry;
    }
  else
    {
      All.SofteningTable[0] = All.SofteningGas;
      All.SofteningTable[1] = All.SofteningHalo;
      All.SofteningTable[2] = All.SofteningDisk;
      All.SofteningTable[3] = All.SofteningBulge;
      All.SofteningTable[4] = All.SofteningStars;
      All.SofteningTable[5] = All.SofteningBndry;
    }

  for(i = 0; i < 6; i++)
    All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];

  All.MinGasHsml = All.MinGasHsmlFractional * All.ForceSoftening[0];
}


/*! This function is used as a comparison kernel in a sort routine. It is
 *  used to group particles in the communication buffer that are going to
 *  be sent to the same CPU.
    */
int data_index_compare(const void *a, const void *b)
{
  if(((struct data_index *) a)->Task < (((struct data_index *) b)->Task))
    return -1;

  if(((struct data_index *) a)->Task > (((struct data_index *) b)->Task))
    return +1;

  if(((struct data_index *) a)->Index < (((struct data_index *) b)->Index))
    return -1;

  if(((struct data_index *) a)->Index > (((struct data_index *) b)->Index))
    return +1;

  return 0;
}

static void msort_dataindex_with_tmp(struct data_index *b, size_t n, struct data_index *t)
{
  struct data_index *tmp;

  struct data_index *b1, *b2;

  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = b;
  b2 = b + n1;

  msort_dataindex_with_tmp(b1, n1, t);
  msort_dataindex_with_tmp(b2, n2, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(b1->Task < b2->Task || (b1->Task == b2->Task && b1->Index <= b2->Index))
	{
	  --n1;
	  *tmp++ = *b1++;
	}
      else
	{
	  --n2;
	  *tmp++ = *b2++;
	}
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * sizeof(struct data_index));

  memcpy(b, t, (n - n2) * sizeof(struct data_index));
}

void mysort_dataindex(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *))
{
  const size_t size = n * s;

  struct data_index *tmp = (struct data_index *) mymalloc(size);

  msort_dataindex_with_tmp((struct data_index *) b, n, tmp);

  myfree(tmp);
}
