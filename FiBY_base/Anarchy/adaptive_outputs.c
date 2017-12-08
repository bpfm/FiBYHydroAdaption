#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include "allvars.h"
#include "proto.h"

#ifdef BG_COOLING
#include "bg_cooling.h"
#endif

#ifdef BG_SFR
#include "bg_vars.h"
#endif

#if defined(BG_SFR) || defined(BG_COOLING)
#include "bg_proto.h"
#endif

/*! \file io.c
 *  \brief Output of an adaptive snapshot file to disk.
 */

#ifdef ADAPTIVE_OUTPUT

static int n_type[6];

static long long ntot_type_all[6];

/*! This function writes a snapshot of the particle distribution to one or
 * several files using Gadget's default file format.  If
 * NumFilesPerSnapshot>1, the snapshot is distributed into several files,
 * which are written simultaneously. Each file contains data from a group of
 * processors of size roughly NTask/NumFilesPerSnapshot.
 */

void write_adaptive_output()
{
  char buf[500];

  int i, j, *temp, n, filenr, gr, ngroups, masterTask, lastTask,nactive;

  size_t bytes;

  double t00, t10;

  if(ThisTask == 0)
    {
	printf("\nwriting adaptive output file... \n");
    }
  t00 = second();

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
      P[i].StepsSinceLastOutput++;

#ifdef MEMDEBUG
  size_t nbefore;

  nbefore = check_for_largest_memory_block();
  printf("Task=%d can alloacte maximal %g MB before I/O of snapshot (ceiling is %g MB)\n", ThisTask,
	 nbefore / (1024.0 * 1024.0), (nbefore + AllocatedBytes) / (1024.0 * 1024.0));
#endif

#if defined(BG_SFR) || defined(BLACK_HOLES)
  rearrange_particle_sequence();
  /* ensures that new tree will be constructed */
  All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TreeDomainUpdateFrequency * All.TotNumPart);
#endif

  if(!(CommBuffer = mymalloc(bytes = All.BufferSize * 1024 * 1024)))
    {
      printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(2);
    }

  /* determine global and local particle numbers */
  for(n = 0; n < 6; n++)
    n_type[n] = 0;

  for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
    if (P[n].StepsSinceLastOutput > All.AO_SubSamplingFactor)
      n_type[P[n].Type]++;


  /* because ntot_type_all[] is of type `long long', we cannot do a simple
   * MPI_Allreduce() to sum the total particle numbers 
   */
  temp = (int *) mymalloc(NTask * 6 * sizeof(int));
  MPI_Allgather(n_type, 6, MPI_INT, temp, 6, MPI_INT, MPI_COMM_WORLD);
  for(i = 0; i < 6; i++)
    {
      ntot_type_all[i] = 0;
      for(j = 0; j < NTask; j++)
	ntot_type_all[i] += temp[j * 6 + i];
    }
  myfree(temp);

  nactive = 0;
  for(n = 0; n < 6; n++)
    nactive += ntot_type_all[n];

  if(nactive > 0) {

    /* assign processors to output files */
    distribute_file(All.NumFilesPerSnapshot, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);
    
    if(All.NumFilesPerSnapshot > 1)
      {
	sprintf(buf, "%s/snapshot_%03d/adaptive/adapt_%03d.%d", All.OutputDir, All.SnapshotFileCount-1, All.SnapshotFileCount-1, filenr);
      }
    else
      {
	sprintf(buf, "%s/snapshot_%03d/adaptive/adapt_%03d", All.OutputDir, All.SnapshotFileCount-1, All.SnapshotFileCount-1);
      }
    
    ngroups = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
    if((All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel))
      ngroups++;
    
    for(gr = 0; gr < ngroups; gr++)
      {
	if((filenr / All.NumFilesWrittenInParallel) == gr)	/* ok, it's this processor's turn */
	  write_adaptive_file(buf, masterTask, lastTask);
	MPI_Barrier(MPI_COMM_WORLD);
      }
  }

  myfree(CommBuffer);
  if(ThisTask == 0)
    printf("done with adaptive snapshot.\n");

#ifdef USE_HDF5_FIX
  H5close();
  hdf5_memory_cleanup();
#endif

#ifdef MEMDEBUG
  size_t nafter;

  nafter = check_for_largest_memory_block();
  if(nafter < nbefore)
    printf("Task=%d can allocate %g MB after snapshot, LESS than before (%g)\n", ThisTask,
	   nafter / (1024.0 * 1024.0), nbefore / (1024.0 * 1024.0));
#endif

  for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
      if (P[n].StepsSinceLastOutput > All.AO_SubSamplingFactor)
		  P[n].StepsSinceLastOutput = 0;

  t10 = second();
  
  if(ThisTask == 0.0)
	 {
		printf(" total adaptive write took %g sec)\n", timediff(t00, t10));
		fflush(stdout);
	 }

}



/*! This function fills the write buffer with particle data. New output blocks can in
 *  principle be added here.
 */
void fill_adaptive_write_buffer(enum iofields blocknr, int *startindex, int pc, int type)
{
  int n, k, pindex, dt_step;

  float *fp;

  double dmax1, dmax2;

#ifdef LONGIDS
  long long *ip;
#else
  int *ip;
#endif

#ifdef PERIODIC
  MyFloat boxSize;
#endif
#ifdef PMGRID
  double dt_gravkick_pm = 0;
#endif
  double dt_gravkick, dt_hydrokick, a3inv = 1, fac1, fac2;

#ifdef BG_COOLING
  double temperature;
#endif

  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      fac1 = 1 / (All.Time * All.Time);
      fac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
    }
  else
    a3inv = fac1 = fac2 = 1;

#ifdef PMGRID
  if(All.ComovingIntegrationOn)
    dt_gravkick_pm =
      get_gravkick_factor(All.PM_Ti_begstep, All.Ti_Current) -
      get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
  else
    dt_gravkick_pm = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;
#endif

  fp = (float *) CommBuffer;
  ip = (int *) CommBuffer;

  pindex = *startindex;

  switch (blocknr)
    {
    case IO_POS:		/* positions */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].Pos[k];
#ifdef PERIODIC
		boxSize = All.BoxSize;
#ifdef LONG_X
		if(k == 0)
		  boxSize = All.BoxSize * LONG_X;
#endif
#ifdef LONG_Y
		if(k == 1)
		  boxSize = All.BoxSize * LONG_Y;
#endif
#ifdef LONG_Z
		if(k == 2)
		  boxSize = All.BoxSize * LONG_Z;
#endif
		while(fp[k] < 0)
		  fp[k] += boxSize;
		while(fp[k] >= boxSize)
		  fp[k] -= boxSize;
#endif
	      }
	    n++;
	    fp += 3;
	  }
      break;

    case IO_VEL:		/* velocities */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    dt_step = TISTEP(P[pindex].TimeBin);

	    if(All.ComovingIntegrationOn)
	      {
		dt_gravkick =
		  get_gravkick_factor(P[pindex].Ti_begstep,
				      All.Ti_Current) -
		  get_gravkick_factor(P[pindex].Ti_begstep, P[pindex].Ti_begstep + dt_step / 2);
		dt_hydrokick =
		  get_hydrokick_factor(P[pindex].Ti_begstep,
				       All.Ti_Current) -
		  get_hydrokick_factor(P[pindex].Ti_begstep, P[pindex].Ti_begstep + dt_step / 2);
	      }
	    else
	      dt_gravkick = dt_hydrokick =
		(All.Ti_Current - (P[pindex].Ti_begstep + dt_step / 2)) * All.Timebase_interval;

	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].Vel[k] + P[pindex].g.GravAccel[k] * dt_gravkick;
		if(P[pindex].Type == 0)
		  fp[k] += SphP[pindex].a.HydroAccel[k] * dt_hydrokick;
	      }
#ifdef PMGRID
	    for(k = 0; k < 3; k++)
	      fp[k] += P[pindex].GravPM[k] * dt_gravkick_pm;
#endif
	    for(k = 0; k < 3; k++)
	      fp[k] *= sqrt(a3inv);

	    n++;
	    fp += 3;
	  }
      break;

    case IO_ID:		/* particle ID */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *ip++ = P[pindex].ID;
	    n++;
	  }
      break;

    case IO_MASS:		/* particle mass */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = P[pindex].Mass;
	    n++;
	  }
      break;

    case IO_U:			/* internal energy */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ =
	      DMAX(All.MinEgySpec,
		   SphP[pindex].Entropy / GAMMA_MINUS1 * pow(SphP[pindex].d.Density * a3inv, GAMMA_MINUS1));
	    n++;
	  }
      break;

    case IO_RHO:		/* density */
      for(n = 0; n < pc; pindex++)
#ifndef BG_SFR
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].d.Density;
	    n++;
	  }
#else
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    if(type == 4)
	      *fp++ = StarP[P[pindex].StarID].GasDensity;
	    else
	      *fp++ = SphP[pindex].d.Density;
	    n++;
	  }
#endif
      break;

    case IO_NE:		/* electron abundance */
      break;

    case IO_NH:		/* neutral hydrogen fraction */
      break;

    case IO_ELECT:
    case IO_HI:
    case IO_HII:
    case IO_HeI:
    case IO_HeII:
    case IO_HeIII:
    case IO_H2I:
    case IO_H2II:
    case IO_HM:
      break;

    case IO_HSML:		/* SPH smoothing length */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = PPP[pindex].Hsml;
	    n++;
	  }
      break;

    case IO_SFR:		/* star formation rate */
#ifdef BG_SFR
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].Sfr;
	    n++;
	  }
#endif
      break;

    case IO_AGE:		/* stellar formation time */
#ifdef BG_SFR
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = StarP[P[pindex].StarID].StarBirthTime;
	    n++;
	  }
#endif
      break;

    case IO_Z:			/* gas and star metallicity */
      break;

    case IO_POT:		/* gravitational potential */
#ifdef OUTPUTPOTENTIAL
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = P[pindex].p.Potential;
	    n++;
	  }
#endif
      break;

    case IO_ACCEL:		/* acceleration */
#ifdef OUTPUTACCELERATION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    for(k = 0; k < 3; k++)
	      fp[k] = fac1 * P[pindex].g.GravAccel[k];
#ifdef PMGRID
	    for(k = 0; k < 3; k++)
	      fp[k] += fac1 * P[pindex].GravPM[k];
#endif
	    if(P[pindex].Type == 0)
	      for(k = 0; k < 3; k++)
		fp[k] += fac2 * SphP[pindex].HydroAccel[k];
	    fp += 3;
	    n++;
	  }
#endif
      break;

    case IO_DTENTR:		/* rate of change of entropy */
#ifdef OUTPUTCHANGEOFENTROPY
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].e.DtEntropy;
	    n++;
	  }
#endif
      break;

    case IO_TSTP:		/* timestep  */
#ifdef OUTPUTTIMESTEP

      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = (P[pindex].Ti_endstep - P[pindex].Ti_begstep) * All.Timebase_interval;
	    n++;
	  }
#endif
      break;

    case IO_BFLD:		/* magnetic field  */
#ifdef MAGNETIC
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = SphP[pindex].B[k];
	    n++;
	  }
#endif
      break;

    case IO_DBDT:		/* rate of change of magnetic field  */
#ifdef DBOUTPUT
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = SphP[pindex].DtB[k];
	    n++;
	  }
#endif
      break;

    case IO_DIVB:		/* divergence of magnetic field  */
#ifdef TRACEDIVB
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].divB;
	    n++;
	  }
#endif
      break;

    case IO_ABVC:		/* artificial viscosity of particle  */
#ifdef TIME_DEP_ART_VISC
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].alpha;
	    n++;
	  }
#endif
      break;

    case IO_AMDC:		/* artificial viscosity of particle  */
#ifdef TIME_DEP_MAGN_DISP
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].Balpha;
	    n++;
	  }
#endif
      break;

    case IO_PHI:		/* divBcleaning fuction of particle  */
#ifdef DIVBCLEANING_DEDNER
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].PhiPred;
	    n++;
	  }
#endif
      break;

    case IO_COOLRATE:		/* current cooling rate of particle  */
      break;

    case IO_CONDRATE:		/* current heating/cooling due to thermal conduction  */
      break;

    case IO_BSMTH:		/* smoothed magnetic field */
      break;

    case IO_DENN:		/* density normalization factor */
      break;

    case IO_EGYPROM:
      break;

    case IO_EGYCOLD:
      break;

    case IO_CR_C0:
      break;

    case IO_CR_Q0:
      break;

    case IO_CR_P0:
      break;

    case IO_CR_E0:
      break;

    case IO_CR_n0:
      break;

    case IO_CR_ThermalizationTime:
      break;

    case IO_CR_DissipationTime:
      break;

    case IO_BHMASS:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = P[pindex].BH_Mass;
	    n++;
	  }
#endif
      break;

    case IO_BHMDOT:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = P[pindex].BH_Mdot;
	    n++;
	  }
#endif
      break;

    case IO_BH_BIRTH_TIME:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = P[pindex].BH_BirthTime;
	    n++;
	  }
#endif
      break;

    case IO_MACH:
#ifdef MACHNUM
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].Shock_MachNumber;
	    n++;
	  }
#endif
      break;

    case IO_DTENERGY:
#ifdef MACHSTATISTIC
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].Shock_DtEnergy;
	    n++;
	  }
#endif
      break;

    case IO_PRESHOCK_DENSITY:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].PreShock_PhysicalDensity;
	    n++;
	  }
#endif
      break;

    case IO_PRESHOCK_ENERGY:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].PreShock_PhysicalEnergy;
	    n++;
	  }
#endif
      break;

    case IO_PRESHOCK_XCR:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].PreShock_XCR;
	    n++;
	  }
#endif
      break;

    case IO_DENSITY_JUMP:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].Shock_DensityJump;
	    n++;
	  }
#endif
      break;

    case IO_ENERGY_JUMP:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].Shock_EnergyJump;
	    n++;
	  }
#endif
      break;

    case IO_CRINJECT:
      break;

    case IO_Zs:
    case IO_iMass:
    case IO_CLDX:
      break;

    case IO_BG_TEMP:
#ifdef BG_COOLING
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    temperature = bg_get_temperature(pindex);
	    *fp++ = (float) temperature;
	    n++;
	  }
#endif
      break;

    case IO_BG_ONEOS:
#ifdef BG_SFR
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = SphP[pindex].OnEOS;
	    n++;
	  }
#endif
      break;

    case IO_BG_SNII_KINETIC_FEEDBACK:
#ifdef BG_SNII_KINETIC_FEEDBACK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    if(type == 4)
	      *fp++ = StarP[P[pindex].StarID].WindFlag;
	    else
	      *fp++ = SphP[pindex].WindFlag;
	    n++;
	  }
#endif
      break;

      /* BG metals */
    case IO_BG_METALS:
#ifdef BG_STELLAR_EVOLUTION
      for(n = 0; n < pc; pindex++)
	{
	  if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	    {
	      if(type == 4)
		for(k = 0; k < BG_NELEMENTS; k++)
		  *fp++ = StarP[P[pindex].StarID].Metals[k] / P[pindex].Mass;
	      else
		for(k = 0; k < BG_NELEMENTS; k++)
		  *fp++ = SphP[pindex].Metals[k] / P[pindex].Mass;
	      n++;
	    }
	}
#endif
      break;

      /* BG metals */
    case IO_BG_METALS_SMOOTHED:
#ifdef BG_METALSMOOTHING
      for(n = 0; n < pc; pindex++)
	{
	  if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	    {
	      if(type == 4)
		for(k = 0; k < BG_NELEMENTS; k++)
		  *fp++ = StarP[P[pindex].StarID].MetalsSmoothed[k] / StarP[P[pindex].StarID].InitialMass;
	      else
		for(k = 0; k < BG_NELEMENTS; k++)
		  *fp++ = SphP[pindex].MetalsSmoothed[k] / P[pindex].Mass;
	      n++;
	    }
	}
#endif
      break;

    case IO_BG_INITIAL_MASS:
#ifdef BG_STELLAR_EVOLUTION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    *fp++ = StarP[P[pindex].StarID].InitialMass;
	    n++;
	  }
#endif
      break;

    case IO_BG_METALLICITY:
#ifdef BG_STELLAR_EVOLUTION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    if(type == 4)
	      *fp++ = StarP[P[pindex].StarID].Metallicity;
	    else
	      *fp++ = SphP[pindex].Metallicity;
	    n++;
	  }
#endif
      break;

    case IO_BG_METALLICITY_SMOOTHED:
#ifdef BG_METALSMOOTHING
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    if(type == 4)
	      *fp++ = StarP[P[pindex].StarID].MetallicitySmoothed;
	    else
	      *fp++ = SphP[pindex].MetallicitySmoothed;
	    n++;
	  }
#endif
      break;

    case IO_BG_IRON_FROM_SNIA:
#ifdef BG_SNIA_IRON
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    if(type == 4)
	      *fp++ = StarP[P[pindex].StarID].IronFromSNIa / StarP[P[pindex].StarID].InitialMass;
	    else
	      *fp++ = SphP[pindex].IronFromSNIa / P[pindex].Mass;
	    n++;
	  }
#endif
      break;

    case IO_BG_IRON_FROM_SNIA_SMOOTHED:
#if defined(BG_SNIA_IRON) && defined(BG_METALSMOOTHING)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    if(type == 4)
	      *fp++ = StarP[P[pindex].StarID].IronFromSNIaSmoothed / StarP[P[pindex].StarID].InitialMass;
	    else
	      *fp++ = SphP[pindex].IronFromSNIaSmoothed / P[pindex].Mass;
	    n++;
	  }
#endif
      break;

    case IO_BG_MAX_ENTROPY:
#ifdef BG_EXTRA_ARRAYS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    if(type == 4)
	      *fp++ = StarP[P[pindex].StarID].MaximumEntropy;
	    else
	      *fp++ = SphP[pindex].MaximumEntropy;
	    n++;
	  }
#endif
      break;

    case IO_BG_MAX_TEMPERATURE:
#ifdef BG_EXTRA_ARRAYS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    if(type == 4)
	      *fp++ = StarP[P[pindex].StarID].MaximumTemperature;
	    else
	      *fp++ = SphP[pindex].MaximumTemperature;
	    n++;
	  }
#endif
      break;

    case IO_BG_TIME_MAX_ENTROPY:
#ifdef BG_EXTRA_ARRAYS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    if(type == 4)
	      *fp++ = StarP[P[pindex].StarID].TimeMaximumEntropy;
	    else
	      *fp++ = SphP[pindex].TimeMaximumEntropy;
	    n++;
	  }
#endif
      break;

    case IO_BG_TIME_MAX_TEMPERATURE:
#ifdef BG_EXTRA_ARRAYS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    if(type == 4)
	      *fp++ = StarP[P[pindex].StarID].TimeMaximumTemperature;
	    else
	      *fp++ = SphP[pindex].TimeMaximumTemperature;
	    n++;
	  }
#endif
      break;

    case IO_BG_METALLICITY_WEIGHTED_REDSHIFT:
#ifdef BG_Z_WEIGHTED_REDSHIFT
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type && P[pindex].StepsSinceLastOutput > All.AO_SubSamplingFactor)
	  {
	    if(type == 4)
	      *fp++ = StarP[P[pindex].StarID].MetallicityWeightedRedshift;
	    else
	      *fp++ = SphP[pindex].MetallicityWeightedRedshift;
	    n++;
	  }
#endif
      break;

    case IO_BG_METALLICITY_WEIGHTED_POTENTIAL:
      break;

    case IO_BG_STELLAR_AGE:
      break;

#ifdef SUBFIND
      //The SubFind additional units - Alan
    case IO_SUBFIND_MASS:
    case IO_SUBFIND_MASSTYPE:
    case IO_SUBFIND_POS:
    case IO_SUBFIND_CMPOS:
    case IO_SUBFIND_CMVEL:
    case IO_SUBFIND_VELDISP:
    case IO_SUBFIND_STELLARVELDISP:
    case IO_SUBFIND_STELLARVELDISPHALFPROJ:
    case IO_SUBFIND_HALFMASS:
    case IO_SUBFIND_HALFMASSPROJ:
    case IO_SUBFIND_VMAX:
    case IO_SUBFIND_VMAXRAD:
    case IO_SUBFIND_SPIN:
    case IO_SUBFIND_PARTPOS:
    case IO_SUBFIND_PARTVEL:
    case IO_SUBFIND_PARTMASS:
      break;
      //The FOF additional units - Alan
    case IO_M_MEAN:
    case IO_R_MEAN:
    case IO_M_CRIT:
    case IO_R_CRIT:
    case IO_M_TOPHAT:
    case IO_R_TOPHAT:
      break;
#endif

    case IO_LASTENTRY:
      endrun(213);
      break;
    }

  *startindex = pindex;
}




/*! This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array.
 */
int get_adaptive_particles_in_block(enum iofields blocknr, int *typelist)
{
  int i, nall, ntot_withmasses, ngas, nstars;

  nall = 0;
  ntot_withmasses = 0;

  for(i = 0; i < 6; i++)
    {
      typelist[i] = 0;

      if(ntot_type_all[i] > 0)
	{
	  nall += ntot_type_all[i];
	  typelist[i] = 1;
	}

      if(All.MassTable[i] == 0)
	ntot_withmasses += ntot_type_all[i];
    }

  ngas = ntot_type_all[0];
  nstars = ntot_type_all[4];


  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
    case IO_TSTP:
    case IO_ID:
    case IO_POT:
    case IO_AO:
      return nall;
      break;

    case IO_MASS:
      for(i = 0; i < 6; i++)
	{
	  typelist[i] = 0;
	  if(All.MassTable[i] == 0 && header.npart[i] > 0)
	    typelist[i] = 1;
	}
      return ntot_withmasses;
      break;

    case IO_U:
#ifndef BG_SFR
    case IO_RHO:
#endif
    case IO_NE:
    case IO_NH:
    case IO_ELECT:
    case IO_HI:
    case IO_HII:
    case IO_HeI:
    case IO_HeII:
    case IO_HeIII:
    case IO_H2I:
    case IO_H2II:
    case IO_HM:
#ifndef BG_SFR
    case IO_HSML:
#endif
    case IO_SFR:
    case IO_DTENTR:
    case IO_BSMTH:
    case IO_BFLD:
    case IO_DBDT:
    case IO_DIVB:
    case IO_ABVC:
    case IO_AMDC:
    case IO_PHI:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_DENN:
    case IO_CR_C0:
    case IO_CR_Q0:
    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
    case IO_MACH:
    case IO_DTENERGY:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
      for(i = 1; i < 6; i++)
	typelist[i] = 0;
      return ngas;
      break;

    case IO_AGE:
      for(i = 0; i < 6; i++)
	if(i != 4)
	  typelist[i] = 0;
      return nstars;

    case IO_Z:
    case IO_EGYPROM:
    case IO_EGYCOLD:
#ifdef BG_SFR
    case IO_HSML:
    case IO_RHO:
#endif
    case IO_BG_METALS:
    case IO_BG_METALS_SMOOTHED:
    case IO_BG_METALLICITY:
    case IO_BG_METALLICITY_SMOOTHED:
    case IO_BG_METALLICITY_WEIGHTED_REDSHIFT:
    case IO_BG_METALLICITY_WEIGHTED_POTENTIAL:
    case IO_BG_IRON_FROM_SNIA:
    case IO_BG_IRON_FROM_SNIA_SMOOTHED:
    case IO_BG_MAX_ENTROPY:
    case IO_BG_MAX_TEMPERATURE:
    case IO_BG_TIME_MAX_ENTROPY:
    case IO_BG_TIME_MAX_TEMPERATURE:
    case IO_BG_SNII_KINETIC_FEEDBACK:
      for(i = 0; i < 6; i++)
	if(i != 0 && i != 4)
	  typelist[i] = 0;
      return ngas + nstars;
      break;

    case IO_BG_INITIAL_MASS:
      for(i = 0; i < 6; i++)
	if(i != 4)
	  typelist[i] = 0;
      return nstars;
      break;

    case IO_BG_TEMP:
    case IO_BG_ONEOS:
      for(i = 0; i < 6; i++)
	if(i != 0)
	  typelist[i] = 0;
      return ngas;
      break;

    case IO_BHMASS:
    case IO_BHMDOT:
    case IO_BH_BIRTH_TIME:
      for(i = 0; i < 6; i++)
	if(i != 5)
	  typelist[i] = 0;
      return header.npart[5];
      break;
    case IO_Zs:
      for(i = 0; i < 6; i++)
	if(i != 0 && i != 4)
	  typelist[i] = 0;
      return ngas + nstars;
      break;
    case IO_iMass:
    case IO_CLDX:
      for(i = 0; i < 6; i++)
	if(i != 4)
	  typelist[i] = 0;
      return nstars;
      break;

    case IO_BG_STELLAR_AGE:
      return 0;
      break;

#ifdef SUBFIND
      //The SubFind additional units - Alan
    case IO_SUBFIND_MASS:
    case IO_SUBFIND_MASSTYPE:
    case IO_SUBFIND_POS:
    case IO_SUBFIND_CMPOS:
    case IO_SUBFIND_HALFMASS:
    case IO_SUBFIND_HALFMASSPROJ:
    case IO_SUBFIND_VMAXRAD:
    case IO_SUBFIND_SPIN:
    case IO_SUBFIND_VMAX:
    case IO_SUBFIND_CMVEL:
    case IO_SUBFIND_VELDISP:
    case IO_SUBFIND_STELLARVELDISP:
    case IO_SUBFIND_STELLARVELDISPHALFPROJ:
    case IO_SUBFIND_PARTPOS:
    case IO_SUBFIND_PARTVEL:
    case IO_SUBFIND_PARTMASS:
      return 0;
      break;

      //The FOF additional units - Alan
    case IO_M_MEAN:
    case IO_M_CRIT:
    case IO_M_TOPHAT:
    case IO_R_MEAN:
    case IO_R_CRIT:
    case IO_R_TOPHAT:
      return 0;
      break;
#endif

    case IO_LASTENTRY:
      endrun(216);
      break;
    }

  endrun(212);
  return 0;
}



void write_adaptive_file(char *fname, int writeTask, int lastTask)
{
  int type, bytes_per_blockelement, npart, nextblock, typelist[6];

  int n_for_this_task, ntask, n, p, pc, offset = 0, task;

  int blockmaxlen, ntot_type[6], nn[6];

  enum iofields blocknr;

  char label[8];

  int bnr;

  int blksize;

  MPI_Request req;

  MPI_Status status;

  FILE *fd = 0;

  int dummy, otherTask;

  char buf[500], gbuf[500];

  double t0, t1;

  int nsteps_thisfile;

#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0,
    hdf5_dataspace_memory, element_grp[6], selement_grp[6];
#ifdef BG_STELLAR_EVOLUTION
  int k;

  hid_t element_dset[BG_NELEMENTS], hdf5_element_dataspace_in_file[BG_NELEMENTS];

  hid_t selement_dset[BG_NELEMENTS], hdf5_selement_dataspace_in_file[BG_NELEMENTS];

  char mbuf[500];
#endif
  hid_t hdf5_paramgrp = 0, hdf5_unitsgrp = 0, hdf5_constgrp = 0, hdf5_timegrp = 0, this_grp;

  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0, hdf5_attribute = 0, hdf5_basegrp=0;

  herr_t hdf5_status;

  hsize_t dims[2], count[2], start[2];
  hsize_t adim[1] = { 6 };

  int rank = 0, pcsum = 0;
#endif

#define SKIP  {my_fwrite(&blksize,sizeof(int),1,fd);}

  /* determine particle numbers of each type in file */

  if(ThisTask == writeTask)
    {
      for(n = 0; n < 6; n++)
	ntot_type[n] = n_type[n];

      for(task = writeTask + 1; task <= lastTask; task++)
	{
	  MPI_Recv(&nn[0], 6, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
	  for(n = 0; n < 6; n++)
	    ntot_type[n] += nn[n];
	}

      for(task = writeTask + 1; task <= lastTask; task++)
	MPI_Send(&ntot_type[0], 6, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Send(&n_type[0], 6, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
      MPI_Recv(&ntot_type[0], 6, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
    }



  /* fill file header */

  for(n = 0; n < 6; n++)
    {
      header.npart[n] = ntot_type[n];
      header.npartTotal[n] = (unsigned int) ntot_type_all[n];
      header.npartTotalHighWord[n] = (unsigned int) (ntot_type_all[n] >> 32);
    }

  for(n = 0; n < 6; n++)
    header.mass[n] = All.MassTable[n];

  header.time = All.Time;

  if(All.ComovingIntegrationOn)
    header.redshift = 1.0 / All.Time - 1;
  else
    header.redshift = 0;

  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;
  header.num_files = All.NumFilesPerSnapshot;
  header.BoxSize = All.BoxSize;
  header.Omega0 = All.Omega0;
  header.OmegaLambda = All.OmegaLambda;
  header.HubbleParam = All.HubbleParam;

  if (ThisTask == writeTask) {
    if (All.NeedFileRefresh == 1) {

      All.NeedFileRefresh = 0;
      /* open file and write header */
      sprintf(buf, "%s.hdf5", fname);
      hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      printf("Creating file : %s\n",buf);
      hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);
      hdf5_paramgrp = H5Gcreate(hdf5_file, "/Parameters", 0);
      hdf5_unitsgrp = H5Gcreate(hdf5_file, "/Units", 0);
      hdf5_constgrp = H5Gcreate(hdf5_file, "/Constants", 0);
      
      write_header_attributes_in_hdf5(hdf5_headergrp);
      write_parameters_attributes_in_hdf5(hdf5_paramgrp);
      write_units_attributes_in_hdf5(hdf5_unitsgrp);
      write_constants_attributes_in_hdf5(hdf5_constgrp);

      H5Gclose(hdf5_headergrp);
      H5Gclose(hdf5_paramgrp);
      H5Gclose(hdf5_unitsgrp);
      H5Gclose(hdf5_constgrp);      

		//Set the number of TimeSteps in this file to 1
		/* NSteps_ThisFile */
		hdf5_dataspace_memory = H5Screate(H5S_SCALAR);
		hdf5_attribute = H5Acreate(hdf5_file, "NSteps_ThisFile", H5T_NATIVE_INT, hdf5_dataspace_memory, H5P_DEFAULT);
		nsteps_thisfile = 1;
		H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &nsteps_thisfile);
		H5Aclose(hdf5_attribute);
		H5Sclose(hdf5_dataspace_memory);

		//Set start and end times for this adaptive file
		hdf5_dataspace_memory = H5Screate(H5S_SCALAR);
		hdf5_attribute = H5Acreate(hdf5_file, "TimeFirst", H5T_NATIVE_DOUBLE, hdf5_dataspace_memory, H5P_DEFAULT);
		H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
		H5Aclose(hdf5_attribute);
		H5Sclose(hdf5_dataspace_memory);


		hdf5_dataspace_memory = H5Screate(H5S_SCALAR);
		hdf5_attribute = H5Acreate(hdf5_file, "TimeLast", H5T_NATIVE_DOUBLE, hdf5_dataspace_memory, H5P_DEFAULT);
		H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
		H5Aclose(hdf5_attribute);
		H5Sclose(hdf5_dataspace_memory);

		//Create the baser group for timestep storage
		hdf5_basegrp = H5Gcreate(hdf5_file, "/Steps", 0);

    } else { //DoNotNeedFileRefresh
      sprintf(buf, "%s.hdf5", fname);
      hdf5_file = H5Fopen(buf, H5F_ACC_RDWR, H5P_DEFAULT);

		//Increment the number of TimeSteps in this file
		hdf5_attribute = H5Aopen_name(hdf5_file, "NSteps_ThisFile");
		H5Aread(hdf5_attribute, H5T_NATIVE_INT, &nsteps_thisfile);
		nsteps_thisfile++;
		H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &nsteps_thisfile);
		H5Aclose(hdf5_attribute);

		//Update the time
		hdf5_attribute = H5Aopen_name(hdf5_file, "TimeLast");
		H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
		H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
		H5Aclose(hdf5_attribute);

		//Open the base group for timestep storage
		hdf5_basegrp = H5Gopen(hdf5_file,"Steps");

    }
  }

  ntask = lastTask - writeTask + 1;

  if(ThisTask == writeTask) {
    printf("Starting variables output\n");

    sprintf(buf,"TimeStep%05d",nsteps_thisfile-1);


    hdf5_timegrp = H5Gcreate(hdf5_basegrp, buf, 0);

    /* NumPart */
    hdf5_dataspace_memory = H5Screate(H5S_SIMPLE);
    H5Sset_extent_simple(hdf5_dataspace_memory, 1, adim, NULL);
    hdf5_attribute = H5Acreate(hdf5_timegrp, "NumPart", H5T_NATIVE_INT, hdf5_dataspace_memory, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npart);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace_memory);

    /* ExpansionFactor */
    hdf5_dataspace_memory = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(hdf5_timegrp, "ExpansionFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace_memory, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace_memory);

    /* Time_GYR */
    double time_in_gyr;
    time_in_gyr = bg_get_elapsed_time(0, header.time, 1);
    hdf5_dataspace_memory = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(hdf5_timegrp, "Time_GYR", H5T_NATIVE_DOUBLE, hdf5_dataspace_memory, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &time_in_gyr);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace_memory);

    /* Redshift */
    hdf5_dataspace_memory = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(hdf5_timegrp, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace_memory, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace_memory);

    for(type = 0; type < 6; type++)
      {
	element_grp[type] = -1;
	selement_grp[type] = -1;
	hdf5_grp[type] = -1;
	if(ntot_type_all[type] > 0)
	  {
	    sprintf(buf, "PartType%d", type);
	    hdf5_grp[type] = H5Gcreate(hdf5_timegrp, buf, 0);
#ifdef BG_STELLAR_EVOLUTION
	    /* group for metallicities */
	    if(type == 0 || type == 4)
	      element_grp[type] = H5Gcreate(hdf5_grp[type], "ElementAbundance", 0);
	    
	    /* group for smoothed metallicities */
	    if(type == 0 || type == 4)
	      selement_grp[type] = H5Gcreate(hdf5_grp[type], "SmoothedElementAbundance", 0);
#endif
	  }
      }

  }

  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum iofields) bnr;

      if(blocknr == IO_LASTENTRY)
	break;

      if(blockpresent(blocknr))
	{
	  t0 = second();

	  if(ThisTask == writeTask)
	    {
	      printf(" writing variable ");
	      get_dataset_name(blocknr, buf);
	      get_group_name(blocknr, gbuf);
	      printf(" %s, blocknr = %d   ", buf, blocknr);
	      fflush(stdout);
	    }

	  bytes_per_blockelement = get_bytes_per_blockelement(blocknr);

	  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / bytes_per_blockelement;

	  npart = get_adaptive_particles_in_block(blocknr, &typelist[0]);
	  if(ThisTask == writeTask)
	    {
	      printf("Npart = %d   ",npart);
	    }



	  if(npart > 0)
	    {
	      for(type = 0; type < 6; type++)
		{
		  if(typelist[type])
 		    {
#ifdef HAVE_HDF5
		      if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
			{
			  switch (get_datatype_in_block(blocknr))
			    {
			    case 0:
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
			      break;
			    case 1:
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
			      break;
			    case 2:
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
			      break;
			    }

			  dims[0] = header.npart[type];
			  dims[1] = get_values_per_blockelement(blocknr);
			  if(dims[1] == 1)
			    rank = 1;
			  else
			    rank = 2;

			  this_grp = hdf5_grp[type];

#ifdef BG_STELLAR_EVOLUTION
			  /* group name */
			  get_group_name(blocknr, gbuf);

			  if(strcmp(gbuf, "ElementAbundance") == 0)
			    this_grp = element_grp[type];

			  if(strcmp(gbuf, "SmoothedElementAbundance") == 0)
			    this_grp = selement_grp[type];

			  /* dataset name */
			  get_dataset_name(blocknr, buf);

			  /* create dataset handle */
			  hdf5_dataset = -1;

			  if(strcmp(buf, "ElementAbundance") == 0)
			    {
			      /* create data sets for each metal */
			      rank = 1;
			      dims[0] = header.npart[type];
			      dims[1] = 1;
			      for(k = 0; k < BG_NELEMENTS; k++)
				{
				  hdf5_element_dataspace_in_file[k] = H5Screate_simple(rank, dims, NULL);
				  strcpy(mbuf, ElementNames[k]);
				  element_dset[k] =
				    H5Dcreate(this_grp, mbuf, hdf5_datatype,
					      hdf5_element_dataspace_in_file[k], H5P_DEFAULT);
				  write_attributes_in_hdf5(blocknr, element_dset[k]);
				}
			    }
			  else if(strcmp(buf, "SmoothedElementAbundance") == 0)
			    {
			      /* create data sets for each smoothed metal */
			      rank = 1;
			      dims[0] = header.npart[type];
			      dims[1] = 1;
			      for(k = 0; k < BG_NELEMENTS; k++)
				{
				  hdf5_selement_dataspace_in_file[k] = H5Screate_simple(rank, dims, NULL);
				  strcpy(mbuf, ElementNames[k]);
				  selement_dset[k] =
				    H5Dcreate(this_grp, mbuf, hdf5_datatype,
					      hdf5_selement_dataspace_in_file[k], H5P_DEFAULT);
				  write_attributes_in_hdf5(blocknr, selement_dset[k]);
				}
			    }
			  else
			    {
#endif /* BG_STELLAR_EVOLUTION */
			      hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
			      hdf5_dataset =
				H5Dcreate(this_grp, buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);

			      if(All.ComovingIntegrationOn)
				write_attributes_in_hdf5(blocknr, hdf5_dataset);
#ifdef BG_STELLAR_EVOLUTION
			    }
#endif /* BG_STELLAR_EVOLUTION */

			  pcsum = 0;
			}
#endif /* HAVE_HDF5 */

		      for(task = writeTask, offset = 0; task <= lastTask; task++)
			{
			  if(task == ThisTask)
			    {
			      n_for_this_task = n_type[type];

			      for(p = writeTask; p <= lastTask; p++)
				if(p != ThisTask)
				  MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK, MPI_COMM_WORLD);
			    }
			  else
			    MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK, MPI_COMM_WORLD,
				     &status);

			  while(n_for_this_task > 0)
			    {
			      pc = n_for_this_task;

			      if(pc > blockmaxlen)
				pc = blockmaxlen;

			      if(ThisTask == task)
				fill_adaptive_write_buffer(blocknr, &offset, pc, type);

			      if(task != writeTask)
				{
				  if(ThisTask == writeTask)
				    MPI_Irecv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task,
					      TAG_PDATA, MPI_COMM_WORLD, &req);

				  if(ThisTask == writeTask || ThisTask == task)
				    {
				      /* exchange a dummy message to make sure that the Isend
				         is placed after the Irecv (to address IBM Bluegene troubles)
				       */
				      if(ThisTask == writeTask)
					otherTask = task;
				      else
					otherTask = writeTask;

				      MPI_Sendrecv(&dummy, 1, MPI_INT, otherTask, TAG_GRAV_A,
						   &dummy, 1, MPI_INT, otherTask, TAG_GRAV_A,
						   MPI_COMM_WORLD, &status);
				    }

				  if(ThisTask == task)
				    MPI_Isend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask,
					      TAG_PDATA, MPI_COMM_WORLD, &req);

				  if(ThisTask == writeTask || ThisTask == task)
				    MPI_Wait(&req, MPI_STATUS_IGNORE);
				}


			      if(ThisTask == writeTask)
				{
				  if(All.SnapFormat == 3)
				    {
#ifdef HAVE_HDF5
#ifdef BG_STELLAR_EVOLUTION
				      if(strcmp(buf, "ElementAbundance") == 0)
					{
					  /* write each metal separately */
					  for(k = 0; k < BG_NELEMENTS; k++)
					    {
					      /* hyperslab in file */
					      rank = 1;
					      start[0] = pcsum;
					      start[1] = 0;
					      count[0] = pc;
					      count[1] = 1;
					      H5Sselect_hyperslab(hdf5_element_dataspace_in_file[k],
								  H5S_SELECT_SET, start, NULL, count, NULL);

					      int n;

					      float *elem, *data;

					      elem = (float *) mymalloc(pc * sizeof(float));
					      data = (float *) CommBuffer;

					      for(n = 0; n < pc; n++)
						elem[n] = data[BG_NELEMENTS * n + k];

					      rank = 1;
					      dims[0] = pc;
					      hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);

					      hdf5_status =
						H5Dwrite(element_dset[k], hdf5_datatype,
							 hdf5_dataspace_memory,
							 hdf5_element_dataspace_in_file[k], H5P_DEFAULT,
							 elem);

					      myfree(elem);

					      H5Sclose(hdf5_dataspace_memory);
					    }
					  pcsum += pc;
					}
				      else if(strcmp(buf, "SmoothedElementAbundance") == 0)
					{
					  /* write each smoothed metal separately */
					  for(k = 0; k < BG_NELEMENTS; k++)
					    {
					      /* hyperslab in file */
					      rank = 1;
					      start[0] = pcsum;
					      start[1] = 0;
					      count[0] = pc;
					      count[1] = 1;
					      H5Sselect_hyperslab(hdf5_selement_dataspace_in_file[k],
								  H5S_SELECT_SET, start, NULL, count, NULL);

					      int n;

					      float *elem, *data;

					      elem = (float *) mymalloc(pc * sizeof(float));
					      data = (float *) CommBuffer;

					      for(n = 0; n < pc; n++)
						elem[n] = data[BG_NELEMENTS * n + k];

					      rank = 1;
					      dims[0] = pc;
					      hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);

					      hdf5_status =
						H5Dwrite(selement_dset[k], hdf5_datatype,
							 hdf5_dataspace_memory,
							 hdf5_selement_dataspace_in_file[k], H5P_DEFAULT,
							 elem);

					      myfree(elem);

					      H5Sclose(hdf5_dataspace_memory);
					    }
					  pcsum += pc;
					}
				      else
					{
#endif /* BG_STELLAR_EVOLUTION */
					  start[0] = pcsum;
					  start[1] = 0;

					  count[0] = pc;
					  count[1] = get_values_per_blockelement(blocknr);
					  pcsum += pc;

					  H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET,
							      start, NULL, count, NULL);

					  dims[0] = pc;
					  dims[1] = get_values_per_blockelement(blocknr);
					  hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);

					  hdf5_status =
					    H5Dwrite(hdf5_dataset, hdf5_datatype,
						     hdf5_dataspace_memory,
						     hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);

					  H5Sclose(hdf5_dataspace_memory);
#ifdef BG_STELLAR_EVOLUTION
					}
#endif /* BG_STELLAR_EVOLUTION */
#endif /* HAVE_HDF5 */
				    }
				  else
				    {
				      my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
				    }
				}

			      n_for_this_task -= pc;
			    }
			}

#ifdef HAVE_HDF5
		      if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
			{
			  if(All.SnapFormat == 3)
			    {
#ifdef BG_STELLAR_EVOLUTION
			      if(strcmp(buf, "ElementAbundance") == 0)
				{
				  for(k = 0; k < BG_NELEMENTS; k++)
				    {
				      H5Dclose(element_dset[k]);
				      H5Sclose(hdf5_element_dataspace_in_file[k]);
				    }
				}
			      else if(strcmp(buf, "SmoothedElementAbundance") == 0)
				{
				  for(k = 0; k < BG_NELEMENTS; k++)
				    {
				      H5Dclose(selement_dset[k]);
				      H5Sclose(hdf5_selement_dataspace_in_file[k]);
				    }
				}
			      else
				{
#endif /* BG_STELLAR_EVOLUTION */
				  H5Dclose(hdf5_dataset);
				  H5Sclose(hdf5_dataspace_in_file);
#ifdef BG_STELLAR_EVOLUTION
				}
#endif /* BG_STELLAR_EVOLUTION */
			      H5Tclose(hdf5_datatype);
			    }
			}
#endif /* HAVE_HDF5 */
		    }
		}
	    }

	  t1 = second();

	  if(ThisTask == writeTask)
	    {
	      printf("  (took %g sec)\n", timediff(t0, t1));
	      fflush(stdout);
	    }
	}

    }

  if(ThisTask == writeTask)
    {
      for(type = 5; type >= 0; type--)
	if(ntot_type_all[type] > 0)
	  {
#ifdef BG_STELLAR_EVOLUTION
	    if(element_grp[type] > 0)
	      H5Gclose(element_grp[type]);
	    
	    if(selement_grp[type] > 0)
	      H5Gclose(selement_grp[type]);
#endif
	    H5Gclose(hdf5_grp[type]);
	  }
      
      H5Gclose(hdf5_timegrp);
      H5Gclose(hdf5_basegrp);
      H5Fclose(hdf5_file);
    }
  
  if(ThisTask == writeTask)
    printf("Done variables output\n\n");

}

#endif
