#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

/*! \file predict.c
 *  \brief drift particles by a small time interval
 *
 *  This function contains code to implement a drift operation on all the
 *  particles, which represents one part of the leapfrog integration scheme.
 */


void reconstruct_timebins(void)
{
  int i, n, prev, bin;

  for(bin = 0; bin < TIMEBINS; bin++)
    {
      TimeBinCount[bin] = 0;
      TimeBinCountSph[bin] = 0;
      FirstInTimeBin[bin] = -1;
      LastInTimeBin[bin] = -1;
#ifdef BG_SFR
      TimeBinSfr[bin] = 0;
#ifdef BG_POPIII
      TimeBinPOPIIISfr[bin] = 0;
#endif
#endif
#ifdef BLACK_HOLES
      TimeBin_BH_mass[bin] = 0;
      TimeBin_BH_dynamicalmass[bin] = 0;
      TimeBin_BH_Mdot[bin] = 0;
      TimeBin_BH_Medd[bin] = 0;
#endif
    }

  for(i = 0; i < NumPart; i++)
    {
      bin = P[i].TimeBin;

      if(TimeBinCount[bin] > 0)
	{
	  PrevInTimeBin[i] = LastInTimeBin[bin];
	  NextInTimeBin[i] = -1;
	  NextInTimeBin[LastInTimeBin[bin]] = i;
	  LastInTimeBin[bin] = i;
	}
      else
	{
	  FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
	  PrevInTimeBin[i] = NextInTimeBin[i] = -1;
	}
      TimeBinCount[bin]++;
      if(P[i].Type == 0)
	TimeBinCountSph[bin]++;

#ifdef BG_SFR
      if(P[i].Type == 0)
	{
	  TimeBinSfr[bin] += SphP[i].Sfr;
#ifdef BG_POPIII
#ifdef BG_METALSMOOTHING
	  if(SphP[i].MetallicitySmoothed < All.POPIII_MetallicityThreshold)
	    TimeBinPOPIIISfr[bin] += SphP[i].Sfr;
#else
	  if(SphP[i].Metallicity < All.POPIII_MetallicityThreshold)
	    TimeBinPOPIIISfr[bin] += SphP[i].Sfr;
#endif
#endif
	}
#endif
#if BLACK_HOLES
      if(P[i].Type == 5)
	{
	  TimeBin_BH_mass[bin] += P[i].BH_Mass;
	  TimeBin_BH_dynamicalmass[bin] += P[i].Mass;
	  TimeBin_BH_Mdot[bin] += P[i].BH_Mdot;
	  TimeBin_BH_Medd[bin] += P[i].BH_Mdot / P[i].BH_Mass;
	}
#endif
    }

  FirstActiveParticle = -1;

  for(n = 0, prev = -1; n < TIMEBINS; n++)
    {
      if(TimeBinActive[n])
	for(i = FirstInTimeBin[n]; i >= 0; i = NextInTimeBin[i])
	  {
	    if(prev == -1)
	      FirstActiveParticle = i;

	    if(prev >= 0)
	      NextActiveParticle[prev] = i;

	    prev = i;
	  }
    }

  if(prev >= 0)
    NextActiveParticle[prev] = -1;

  long long sum1, sum2;
  int *temp;

  /* MPI_Allreduce(&NumForceUpdate, &sum1, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); */

  temp = (int *) mymalloc(NTask * sizeof(int));

  MPI_Allgather(&NumForceUpdate, 1, MPI_INT, temp, 1, MPI_INT, MPI_COMM_WORLD);

  sum1 = 0;
  for(i = 0; i < NTask; i++)
    sum1 += temp[i];

  myfree(temp);

  for(i = FirstActiveParticle, NumForceUpdate = 0; i >= 0; i = NextActiveParticle[i])
    {
      NumForceUpdate++;
      if(i >= NumPart)
	{
	  printf("Bummer i=%d\n", i);
	  endrun(12);
	}
    }

  /* MPI_Allreduce(&NumForceUpdate, &sum2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); */

  temp = (int *) mymalloc(NTask * sizeof(int));

  MPI_Allgather(&NumForceUpdate, 1, MPI_INT, temp, 1, MPI_INT, MPI_COMM_WORLD);

  sum2 = 0;
  for(i = 0; i < NTask; i++)
    sum2 += temp[i];

  myfree(temp);

  if(ThisTask == 0)
    {
      printf("sum1=%d%09d sum2=%d%09d\n",
	     (int) (sum1 / 1000000000), (int) (sum1 % 1000000000),
	     (int) (sum2 / 1000000000), (int) (sum2 % 1000000000));
    }

  if(sum1 != sum2 && All.NumCurrentTiStep > 0)
    endrun(121);
}



void drift_particle(int i, int time1)
{
  int j, time0, dt_step;

  double dt_drift, dt_gravkick, dt_hydrokick, dt_entr;

  time0 = P[i].Ti_current;

  if(time1 < time0)
    {
      printf("i=%d time0=%d time1=%d\n", i, time0, time1);
      endrun(12);
    }

  if(time1 == time0)
    return;

  if(All.ComovingIntegrationOn)
    {
      dt_drift = get_drift_factor(time0, time1);
      dt_gravkick = get_gravkick_factor(time0, time1);
      dt_hydrokick = get_hydrokick_factor(time0, time1);
    }
  else
    {
      dt_drift = dt_gravkick = dt_hydrokick = (time1 - time0) * All.Timebase_interval;
    }


  for(j = 0; j < 3; j++)
    P[i].Pos[j] += P[i].Vel[j] * dt_drift;

  if(P[i].Type == 0)
    {
#ifdef PMGRID
      for(j = 0; j < 3; j++)
	SphP[i].VelPred[j] +=
	  (P[i].g.GravAccel[j] + P[i].GravPM[j]) * dt_gravkick + SphP[i].a.HydroAccel[j] * dt_hydrokick;
#else
      for(j = 0; j < 3; j++)
	SphP[i].VelPred[j] += P[i].g.GravAccel[j] * dt_gravkick + SphP[i].a.HydroAccel[j] * dt_hydrokick;
#endif

      SphP[i].d.Density *= exp(-SphP[i].v.DivVel * dt_drift);
      PPP[i].Hsml *= exp(NUMDIMS_INV * SphP[i].v.DivVel * dt_drift);

	/*
      if(fact > 1.2)
	fact = 1.2;
      PPP[i].Hsml *= fact;
	*/

      if(PPP[i].Hsml < All.MinGasHsml)
	PPP[i].Hsml = All.MinGasHsml;

      dt_step = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0);
      dt_entr = (time1 - (P[i].Ti_begstep + dt_step / 2)) * All.Timebase_interval;

      SphP[i].EntropyPred = SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr;

#ifdef PRESSURE_ENTROPY_SPH
      SphP[i].cky.WeightedDensity *= exp(-SphP[i].v.DivVel * dt_drift);
      SphP[i].EntropyVarPred = pow(SphP[i].EntropyPred, GAMMA_INV);
#endif

#ifndef PRESSURE_ENTROPY_SPH
      SphP[i].Pressure = SphP[i].EntropyPred * pow(SphP[i].d.Density, GAMMA);
#else
      SphP[i].Pressure = SphP[i].EntropyPred * pow(SphP[i].cky.WeightedDensity, GAMMA);
#endif

      /*
      if(SphP[i].Pressure <= 0)
	{
	  printf("[PREDICT %2d] PID=%d, pressure=%g, entropy=%g, dtentropy=%g, rho=%g\n",
		 ThisTask, P[i].ID, SphP[i].Pressure, SphP[i].Entropy, SphP[i].e.DtEntropy, SphP[i].d.Density);
	  printf("[PREDICT %2d] PID=%d, dt_step=%d, time1=%d, dt_entr=%g, ti_begstep=%d, h=%g\n",
		 ThisTask, P[i].ID, dt_step, time1, dt_entr, P[i].Ti_begstep, PPP[i].Hsml);
	}
      */
    }

  P[i].Ti_current = time1;
}



void move_particles(int time1)
{
  int i;

  if(ThisTask == 0)
    printf("MOVE\n");

  for(i = 0; i < NumPart; i++)
    drift_particle(i, time1);
}



/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef PERIODIC
void do_box_wrapping(void)
{
  int i, j;

  double boxsize[3];

  for(j = 0; j < 3; j++)
    boxsize[j] = All.BoxSize;

#ifdef LONG_X
  boxsize[0] *= LONG_X;
#endif
#ifdef LONG_Y
  boxsize[1] *= LONG_Y;
#endif
#ifdef LONG_Z
  boxsize[2] *= LONG_Z;
#endif

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      {
	while(P[i].Pos[j] < 0)
	  P[i].Pos[j] += boxsize[j];

	while(P[i].Pos[j] >= boxsize[j])
	  P[i].Pos[j] -= boxsize[j];
      }
}
#endif
