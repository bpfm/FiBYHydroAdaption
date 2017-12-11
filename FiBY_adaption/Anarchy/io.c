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
 *  \brief Output of a snapshot file to disk.
 */

static int n_type[6];

static long long ntot_type_all[6];

/*! This function writes a snapshot of the particle distribution to one or
 * several files using Gadget's default file format.  If
 * NumFilesPerSnapshot>1, the snapshot is distributed into several files,
 * which are written simultaneously. Each file contains data from a group of
 * processors of size roughly NTask/NumFilesPerSnapshot.
 */
void savepositions(int num, int mode)
{
  char buf[500];

  int i, j, *temp, n, filenr, gr, ngroups, masterTask, lastTask;

  size_t bytes;


#ifdef BG_MOL_NETWORK
#if defined(LW_BACKGROUND) || defined(LW_LOCAL)
      lw(1);
#endif
#endif


  if(ThisTask == 0)
    {
      if(mode > 0)
	printf("\nwriting restart snapshot file... \n");
      else
	printf("\nwriting snapshot file... \n");
    }

#ifdef MEMDEBUG
  size_t nbefore;

  nbefore = check_for_largest_memory_block();
  printf("Task=%d can alloacte maximal %g MB before I/O of snapshot (ceiling is %g MB)\n", ThisTask,
	 nbefore / (1024.0 * 1024.0), (nbefore + AllocatedBytes) / (1024.0 * 1024.0));
#endif

  if(ThisTask == 0)
    {
      if(mode > 0)
	sprintf(buf, "%s/restart_%03d", All.OutputDir, num);
      else
	sprintf(buf, "%s/snapshot_%03d", All.OutputDir, num);

      mkdir(buf, 02755);

#ifdef ADAPTIVE_OUTPUT
      sprintf(buf, "%s/snapshot_%03d/adaptive", All.OutputDir, num);
      mkdir(buf, 02755);
      All.NeedFileRefresh = 1;
#endif
    }

  CPU_Step[CPU_MISC] += measure_time();

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

  for(n = 0; n < NumPart; n++)
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


  /* assign processors to output files */
  distribute_file(All.NumFilesPerSnapshot, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

  if(All.NumFilesPerSnapshot > 1)
    {
      if(mode > 0)
	sprintf(buf, "%s/restart_%03d/%s_%03d.%d", All.OutputDir, num, All.SnapshotFileBase, num, filenr);
      else
	sprintf(buf, "%s/snapshot_%03d/%s_%03d.%d", All.OutputDir, num, All.SnapshotFileBase, num, filenr);
    }
  else
    {
      if(mode > 0)
	sprintf(buf, "%s/restart_%03d/%s_%03d", All.OutputDir, num, All.SnapshotFileBase, num);
      else
	sprintf(buf, "%s/snapshot_%03d/%s_%03d", All.OutputDir, num, All.SnapshotFileBase, num);
    }

  ngroups = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
  if((All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel))
    ngroups++;

  for(gr = 0; gr < ngroups; gr++)
    {
      if((filenr / All.NumFilesWrittenInParallel) == gr)	/* ok, it's this processor's turn */
	write_file(buf, masterTask, lastTask);
      MPI_Barrier(MPI_COMM_WORLD);
    }

  myfree(CommBuffer);

  if(ThisTask == 0)
    printf("done with snapshot.\n");

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

  CPU_Step[CPU_SNAPSHOT] += measure_time();

#ifdef FOF
  if(mode == 0)
    {
      if(ThisTask == 0)
	printf("\ncomputing group catalogue...\n");

      fof_fof(num);

      if(ThisTask == 0)
	printf("done with group catalogue.\n");
    }

  CPU_Step[CPU_FOF] += measure_time();
#endif
}



/*! This function fills the write buffer with particle data. New output blocks can in
 *  principle be added here.
 */
void fill_write_buffer(enum iofields blocknr, int *startindex, int pc, int type)
{
  int n, k, pindex, dt_step;

  float *fp;

  double dmax1, dmax2;
  MyIDType *ip;

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
  ip = (MyIDType *) CommBuffer;

  pindex = *startindex;

  switch (blocknr)
    {
    case IO_POS:		/* positions */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
	  {
	    *ip++ = P[pindex].ID;
	    n++;
	  }
      break;

    case IO_MASS:		/* particle mass */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].Mass;
	    n++;
	  }
      break;

    case IO_AO:		/* particle mass */
#ifdef ADAPTIVE_OUTPUT
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].StepsSinceLastOutput;
	    n++;
	  }
#endif
      break;


    case IO_U:			/* internal energy */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].d.Density;
	    n++;
	  }
#else
	if(P[pindex].Type == type)
	  {
	    if(type == 4)
	      *fp++ = StarP[P[pindex].StarID].GasDensity;
	    else
	      *fp++ = SphP[pindex].d.Density;
	    n++;
	  }
#endif
      break;

    case IO_WRHO:		/* weighted density */
#ifdef PRESSURE_ENTROPY_SPH
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].cky.WeightedDensity;
	    n++;
	  }
#endif
      break;

    case  IO_BG_e:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].n_e;
	    n++;
	  }
#endif
	break;

    case  IO_BG_HII:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].x_Hp;
	    n++;
	  }
#endif
	break;

    case  IO_BG_Hminus:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].x_Hm;
	    n++;
	  }
#endif
	break;

    case  IO_BG_H2I:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].x_H2;
	    n++;
	  }
#endif
	break;

    case  IO_BG_H2II:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].x_H2p;
	    n++;
	  }
#endif
	break;

    case  IO_BG_HeII:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].x_Hep;
	    n++;
	  }
#endif
	break;

    case  IO_BG_HeIII:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].x_Hepp;
	    n++;
	  }
#endif
	break;

    case  IO_BG_DII:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].x_Dp;
	    n++;
	  }
#endif
	break;

    case  IO_BG_HD:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].x_HD;
	    n++;
	  }
#endif
	break;

/*     case  IO_BG_DUST: */
/* #ifdef BG_DUST */
/*       for(n = 0; n < pc; pindex++) */
/*         if(P[pindex].Type == type) */
/*           { */
/*             *fp++ = SphP[pindex].Dusticity; */
/*             n++; */
/*           } */
/* #endif */
/*       break; */

/*     case  IO_BG_DUST_SMOOTHED: */
/* #ifdef BG_METALSMOOTHING */
/* #ifdef BG_DUST */
/*       for(n = 0; n < pc; pindex++) */
/* 	if(P[pindex].Type == type) */
/*           { */
/*             *fp++ = SphP[pindex].DusticitySmoothed; */
/*             n++; */
/*           } */
/* #endif */
/* #endif */
/*       break; */

/*     case IO_BG_DUST_SIZE: */
/* #ifdef BG_DUST_DESTRUCTION_SUBLIMATION */
/*       for(n = 0; n < pc; pindex++) */
/* 	if(P[pindex].Type == type) */
/*           { */
/*             *fp++ = SphP[pindex].DustGrainSize; */
/*             n++; */
/*           } */
/* #endif */
/*       break; */

case IO_HSML:		/* SPH smoothing length */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = PPP[pindex].Hsml;
	    n++;
	  }
      break;

    case IO_SFR:		/* star formation rate */
#ifdef BG_SFR
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Sfr;
	    n++;
	  }
#endif
      break;

    case IO_AGE:		/* stellar formation time */
#ifdef BG_SFR
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].p.Potential;
	    n++;
	  }
#endif
      break;

    case IO_ACCEL:		/* acceleration */
#ifdef OUTPUTACCELERATION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].e.DtEntropy;
	    n++;
	  }
#endif
      break;

    case IO_TSTP:		/* timestep  */
#ifdef OUTPUTTIMESTEP
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = (P[pindex].TimeBin ? (1 << P[pindex].TimeBin) : 0) * All.Timebase_interval;
	    n++;
	  }
#endif
      break;

    case IO_BFLD:		/* magnetic field  */
#ifdef MAGNETIC
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].divB;
	    n++;
	  }
#endif
      break;

    case IO_ABVC:		/* artificial viscosity of particle  */
#ifdef TIME_DEP_ART_VISC
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].alpha;
	    n++;
	  }
#endif
      break;

    case IO_AMDC:		/* artificial viscosity of particle  */
#ifdef TIME_DEP_MAGN_DISP
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Balpha;
	    n++;
	  }
#endif
      break;

    case IO_PHI:		/* divBcleaning fuction of particle  */
#ifdef DIVBCLEANING_DEDNER
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].BH_Mass;
	    n++;
	  }
#endif
      break;

    case IO_BHMDOT:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].BH_Mdot;
	    n++;
	  }
#endif
      break;

    case IO_BH_BIRTH_TIME:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].BH_BirthTime;
	    n++;
	  }
#endif
      break;

    case IO_BH_ENERGY:
#if defined(BLACK_HOLES) && defined(BH_THERMALFEEDBACK)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].BH_Energy;
	    n++;
	  }
#endif
      break;

    case IO_MACH:
#ifdef MACHNUM
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_MachNumber;
	    n++;
	  }
#endif
      break;

    case IO_DTENERGY:
#ifdef MACHSTATISTIC
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_DtEnergy;
	    n++;
	  }
#endif
      break;

    case IO_PRESHOCK_DENSITY:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].PreShock_PhysicalDensity;
	    n++;
	  }
#endif
      break;

    case IO_PRESHOCK_ENERGY:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].PreShock_PhysicalEnergy;
	    n++;
	  }
#endif
      break;

    case IO_PRESHOCK_XCR:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].PreShock_XCR;
	    n++;
	  }
#endif
      break;

    case IO_DENSITY_JUMP:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_DensityJump;
	    n++;
	  }
#endif
      break;

    case IO_ENERGY_JUMP:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].OnEOS;
	    n++;
	  }
#endif
      break;

    case IO_BG_SNII_KINETIC_FEEDBACK:
#ifdef BG_SNII_KINETIC_FEEDBACK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
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
	  if(P[pindex].Type == type)
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
	  if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
	  {
	    *fp++ = StarP[P[pindex].StarID].InitialMass;
	    n++;
	  }
#endif
      break;

    case IO_BG_METALLICITY:
#ifdef BG_STELLAR_EVOLUTION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
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
	if(P[pindex].Type == type)
	  {
	    if(type == 4)
	      *fp++ = StarP[P[pindex].StarID].MetallicityWeightedRedshift;
	    else
	      *fp++ = SphP[pindex].MetallicityWeightedRedshift;
	    n++;
	  }
#endif
      break;

/*     case IO_LW_popII: */
/* #ifdef LW_LOCAL */
/*       for(n = 0; n < pc; pindex++) */
/* 	if(P[pindex].Type == type) */
/* 	  { */
/* 	    *fp++ = SphP[pindex].LWRadiation_popII; */
/* 	    n++; */
/* 	  } */
/* #endif */
/*       break; */

/*     case IO_LW_popIII: */
/* #ifdef LW_LOCAL */
/*       for(n = 0; n < pc; pindex++) */
/* 	if(P[pindex].Type == type) */
/* 	  { */
/* 	    *fp++ = SphP[pindex].LWRadiation_popIII; */
/* 	    n++; */
/* 	  } */
/* #endif */
/*       break; */

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




/*! This function tells the size of one data entry in each of the blocks
 *  defined for the output file.
 */
int get_bytes_per_blockelement(enum iofields blocknr)
{
  int bytes_per_blockelement = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_BSMTH:
    case IO_ACCEL:
    case IO_BFLD:
    case IO_DBDT:
      bytes_per_blockelement = 3 * sizeof(float);
      break;

    case IO_ID:
      bytes_per_blockelement = sizeof(MyIDType);
      break;
    case IO_AO:
      bytes_per_blockelement = sizeof(int);
    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_WRHO:
    case IO_HSML:
    case IO_SFR:
    case IO_AGE:
    case IO_POT:
    case IO_DTENTR:
    case IO_TSTP:
    case IO_DIVB:
    case IO_ABVC:
    case IO_AMDC:
    case IO_PHI:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_DENN:
    case IO_EGYPROM:
    case IO_EGYCOLD:
    case IO_CR_C0:
    case IO_CR_Q0:
    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
    case IO_BHMASS:
    case IO_BH_ENERGY:
    case IO_BHMDOT:
    case IO_BH_BIRTH_TIME:
    case IO_MACH:
    case IO_DTENERGY:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
    case IO_iMass:
    case IO_CLDX:
      bytes_per_blockelement = sizeof(float);
      break;

    case IO_Z:
      bytes_per_blockelement = sizeof(float);
      break;

    case IO_Zs:
      bytes_per_blockelement = 0;
      break;

    case IO_BG_TEMP:
#ifdef BG_COOLING
      bytes_per_blockelement = sizeof(float);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_BG_ONEOS:
#ifdef BG_SFR
      bytes_per_blockelement = sizeof(float);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_BG_SNII_KINETIC_FEEDBACK:
#ifdef BG_SNII_KINETIC_FEEDBACK
      bytes_per_blockelement = sizeof(float);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_BG_METALS:
#ifdef BG_STELLAR_EVOLUTION
      bytes_per_blockelement = (BG_NELEMENTS) * sizeof(float);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_BG_METALS_SMOOTHED:
#ifdef BG_METALSMOOTHING
      bytes_per_blockelement = (BG_NELEMENTS) * sizeof(float);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_BG_INITIAL_MASS:
    case IO_BG_METALLICITY:
#ifdef BG_STELLAR_EVOLUTION
      bytes_per_blockelement = sizeof(float);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_BG_METALLICITY_SMOOTHED:
#ifdef BG_METALSMOOTHING
      bytes_per_blockelement = sizeof(float);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_BG_IRON_FROM_SNIA:
#ifdef BG_SNIA_IRON
      bytes_per_blockelement = sizeof(float);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_BG_IRON_FROM_SNIA_SMOOTHED:
#if defined(BG_SNIA_IRON) && defined(BG_METALSMOOTHING)
      bytes_per_blockelement = sizeof(float);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_BG_MAX_ENTROPY:
    case IO_BG_MAX_TEMPERATURE:
    case IO_BG_TIME_MAX_ENTROPY:
    case IO_BG_TIME_MAX_TEMPERATURE:
#ifdef BG_EXTRA_ARRAYS
      bytes_per_blockelement = sizeof(float);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_BG_METALLICITY_WEIGHTED_REDSHIFT:
#ifdef BG_Z_WEIGHTED_REDSHIFT
      bytes_per_blockelement = sizeof(float);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_BG_METALLICITY_WEIGHTED_POTENTIAL:
      bytes_per_blockelement = 0;
      break;

    case IO_BG_STELLAR_AGE:
      bytes_per_blockelement = 0;
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
      bytes_per_blockelement = 0;
      break;
      //The FOF additional units - Alan
    case IO_M_MEAN:
    case IO_R_MEAN:
    case IO_M_CRIT:
    case IO_R_CRIT:
    case IO_M_TOPHAT:
    case IO_R_TOPHAT:
      bytes_per_blockelement = 0;
      break;
#endif

    case IO_BG_e:
    case IO_BG_HII:
    case IO_BG_Hminus:
    case IO_BG_H2I:
    case IO_BG_H2II:
    case IO_BG_HeII:
    case IO_BG_HeIII:
    case IO_BG_DII:
    case IO_BG_HD:
#ifdef BG_MOL_NETWORK
      bytes_per_blockelement = sizeof(float);
#else
      bytes_per_blockelement = 0;
#endif
      break;

/*     case IO_BG_DUST: */
/* #ifdef BG_DUST */
/*       bytes_per_blockelement = sizeof(float); */
/* #else */
/*       bytes_per_blockelement = 0; */
/* #endif */
/*       break; */

/*     case IO_BG_DUST_SMOOTHED: */
/* #ifdef BG_DUST */
/* #ifdef BG_METALSMOOTHING */
/*       bytes_per_blockelement = sizeof(float); */
/* #endif */
/* #else */
/*       bytes_per_blockelement = 0; */
/* #endif */
/*       break; */

/*     case IO_BG_DUST_SIZE: */
/* #ifdef BG_DUST_DESTRUCTION_SUBLIMATION */
/*       bytes_per_blockelement = sizeof(float); */
/* #else */
/*       bytes_per_blockelement = 0; */
/* #endif */
/*       break; */

/*     case IO_LW_popII: */
/*     case IO_LW_popIII: */
/*       bytes_per_blockelement = sizeof(float); */
/*       break; */

    case IO_LASTENTRY:
      endrun(214);
      break;
    }

  return bytes_per_blockelement;
}

int get_datatype_in_block(enum iofields blocknr)
{
  int typekey;

  switch (blocknr)
    {
    case IO_ID:
#ifndef LONGIDS
      typekey = 0;		/* native uint */
#else
      typekey = 2;		/* native uint64 */
#endif
      break;
    default:
      typekey = 1;		/* native float */
      break;
    }

  return typekey;
}



int get_values_per_blockelement(enum iofields blocknr)
{
  int values = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_BSMTH:
    case IO_ACCEL:
    case IO_BFLD:
    case IO_DBDT:
      values = 3;
      break;

    case IO_ID:
    case IO_MASS:
    case IO_AO:
    case IO_U:
    case IO_RHO:
    case IO_WRHO:
    case IO_HSML:
    case IO_SFR:
    case IO_AGE:
    case IO_POT:
    case IO_DTENTR:
    case IO_TSTP:
    case IO_DIVB:
    case IO_ABVC:
    case IO_AMDC:
    case IO_PHI:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_DENN:
    case IO_EGYPROM:
    case IO_EGYCOLD:
    case IO_CR_C0:
    case IO_CR_Q0:
    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
    case IO_BHMASS:
    case IO_BH_ENERGY:
    case IO_BHMDOT:
    case IO_BH_BIRTH_TIME:
    case IO_MACH:
    case IO_DTENERGY:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
    case IO_iMass:
    case IO_BG_INITIAL_MASS:
    case IO_BG_METALLICITY:
    case IO_BG_METALLICITY_SMOOTHED:
    case IO_BG_IRON_FROM_SNIA:
    case IO_BG_IRON_FROM_SNIA_SMOOTHED:
    case IO_BG_MAX_ENTROPY:
    case IO_BG_MAX_TEMPERATURE:
    case IO_BG_TIME_MAX_ENTROPY:
    case IO_BG_TIME_MAX_TEMPERATURE:
    case IO_BG_METALLICITY_WEIGHTED_REDSHIFT:
    case IO_BG_METALLICITY_WEIGHTED_POTENTIAL:
    case IO_CLDX:
    case IO_BG_e:
    case IO_BG_HII:
    case IO_BG_Hminus:
    case IO_BG_H2I:
    case IO_BG_H2II:
    case IO_BG_HeII:
    case IO_BG_HeIII:
    case IO_BG_DII:
    case IO_BG_HD:
/*     case IO_BG_DUST: */
/*     case IO_BG_DUST_SMOOTHED: */
/*     case IO_BG_DUST_SIZE: */
      values = 1;
      break;

    case IO_Z:
      values = 1;
      break;

    case IO_Zs:
      values = 0;
      break;

    case IO_BG_TEMP:
#ifdef BG_COOLING
      values = 1;
#else
      values = 0;
#endif
      break;

    case IO_BG_ONEOS:
#ifdef BG_SFR
      values = 1;
#else
      values = 0;
#endif
      break;

    case IO_BG_SNII_KINETIC_FEEDBACK:
#ifdef BG_SNII_KINETIC_FEEDBACK
      values = 1;
#else
      values = 0;
#endif
      break;

    case IO_BG_METALS:
#ifdef BG_STELLAR_EVOLUTION
      values = BG_NELEMENTS;
#else
      values = 0;
#endif

    case IO_BG_METALS_SMOOTHED:
#ifdef BG_METALSMOOTHING
      values = BG_NELEMENTS;
#else
      values = 0;
#endif
      break;

/*     case IO_LW_popII: */
/*     case IO_LW_popIII: */
/*       values = 1; */
/*       break; */

    case IO_BG_STELLAR_AGE:
      values = 0;
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
      values = 0;
      break;
      //The FOF additional units - Alan
    case IO_M_MEAN:
    case IO_R_MEAN:
    case IO_M_CRIT:
    case IO_R_CRIT:
    case IO_M_TOPHAT:
    case IO_R_TOPHAT:
      values = 0;
      break;
#endif

    case IO_LASTENTRY:
      endrun(215);
      break;
    }
  return values;
}




/*! This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array.
 */
int get_particles_in_block(enum iofields blocknr, int *typelist)
{
  int i, nall, ntot_withmasses, ngas, nstars;

  nall = 0;
  ntot_withmasses = 0;

  for(i = 0; i < 6; i++)
    {
      typelist[i] = 0;

      if(header.npart[i] > 0)
	{
	  nall += header.npart[i];
	  typelist[i] = 1;
	}

      if(All.MassTable[i] == 0)
	ntot_withmasses += header.npart[i];
    }

  ngas = header.npart[0];
  nstars = header.npart[4];


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
    case IO_WRHO:
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
    case IO_BG_e:
    case IO_BG_HII:
    case IO_BG_Hminus:
    case IO_BG_H2I:
    case IO_BG_H2II:
    case IO_BG_HeII:
    case IO_BG_HeIII:
    case IO_BG_DII:
    case IO_BG_HD:
/*     case IO_BG_DUST: */
/*     case IO_BG_DUST_SMOOTHED: */
/*     case IO_BG_DUST_SIZE: */
      for(i = 0; i < 6; i++)
	if(i != 0)
	  typelist[i] = 0;
      return ngas;
      break;

    case IO_BHMASS:
    case IO_BH_ENERGY:
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

/*     case IO_LW_popII: */
/*     case IO_LW_popIII: */
/*       for(i = 0; i < 6; i++) */
/* 	if(i != 0) */
/* 	  typelist[i] = 0; */
/*       return ngas; */
/*       break; */

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



/*! This function tells whether a block in the output file is present or not.
 */
int blockpresent(enum iofields blocknr)
{
  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ID:
    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_WRHO:
    case IO_HSML:
      return 1;			/* always present */


    case IO_SFR:
    case IO_AGE:
    case IO_Z:
      if(All.StarformationOn == 0)
	return 0;
      else
	{
#ifdef BG_SFR
	  if(blocknr == IO_SFR)
	    return 1;
#endif
#ifdef BG_SFR
	  if(blocknr == IO_AGE)
	    return 1;
#endif
	}
      return 0;
      break;


    case IO_POT:
#ifdef OUTPUTPOTENTIAL
      return 1;
#else
      return 0;
#endif

    case IO_AO:
#ifdef ADAPTIVE_OUTPUT
      return 1;
#else
      return 0;
#endif


    case IO_iMass:
    case IO_CLDX:
      return 0;
      break;


    case IO_ACCEL:
#ifdef OUTPUTACCELERATION
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DTENTR:
#ifdef OUTPUTCHANGEOFENTROPY
      return 1;
#else
      return 0;
#endif
      break;


    case IO_TSTP:
#ifdef IO_TSTP
      return 1;
#else
      return 0;
#endif


    case IO_BFLD:
#ifdef MAGNETIC
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DBDT:
#ifdef DBOUTPUT
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DIVB:
#ifdef TRACEDIVB
      return 1;
#else
      return 0;
#endif
      break;


    case IO_ABVC:
#ifdef TIME_DEP_ART_VISC
      return 1;
#else
      return 0;
#endif
      break;


    case IO_AMDC:
#ifdef TIME_DEP_MAGN_DISP
      return 1;
#else
      return 0;
#endif
      break;


    case IO_PHI:
#ifdef DIVBCLEANING_DEDNER
      return 1;
#else
      return 0;
#endif
      break;


    case IO_COOLRATE:
      return 0;
      break;


    case IO_CONDRATE:
      return 0;
      break;


    case IO_BSMTH:
      return 0;
      break;


    case IO_DENN:
      return 0;
      break;


    case IO_EGYPROM:
    case IO_EGYCOLD:
      return 0;
      break;


    case IO_BHMASS:
    case IO_BHMDOT:
    case IO_BH_BIRTH_TIME:
#ifdef BLACK_HOLES
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BH_ENERGY:
#if defined(BLACK_HOLES) && defined (BH_THERMALFEEDBACK)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_MACH:
#ifdef MACHNUM
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DTENERGY:
#ifdef MACHSTATISTIC
      return 1;
#else
      return 0;
#endif
      break;


    case IO_CR_C0:
    case IO_CR_Q0:
    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
      return 0;
      break;


    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      return 1;
#else
      return 0;
#endif
      break;


    case IO_CRINJECT:
#ifdef CR_OUTPUT_INJECTION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BG_TEMP:
#ifdef BG_COOLING
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BG_ONEOS:
#ifdef BG_SFR
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BG_SNII_KINETIC_FEEDBACK:
#ifdef BG_SNII_KINETIC_FEEDBACK
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BG_INITIAL_MASS:
    case IO_BG_METALLICITY:
#if defined(BG_STELLAR_EVOLUTION)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BG_METALS:
#ifdef BG_STELLAR_EVOLUTION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BG_IRON_FROM_SNIA:
#ifdef BG_SNIA_IRON
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BG_IRON_FROM_SNIA_SMOOTHED:
#if defined(BG_SNIA_IRON) && defined(BG_METALSMOOTHING)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BG_MAX_ENTROPY:
    case IO_BG_MAX_TEMPERATURE:
    case IO_BG_TIME_MAX_ENTROPY:
    case IO_BG_TIME_MAX_TEMPERATURE:
#ifdef BG_EXTRA_ARRAYS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BG_METALS_SMOOTHED:
    case IO_BG_METALLICITY_SMOOTHED:
#ifdef BG_METALSMOOTHING
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BG_METALLICITY_WEIGHTED_REDSHIFT:
#ifdef BG_Z_WEIGHTED_REDSHIFT
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BG_METALLICITY_WEIGHTED_POTENTIAL:
      return 0;
      break;

    case IO_Zs:
      return 0;
      break;

      /* this is used only for units and a / h scaling */
    case IO_BG_STELLAR_AGE:
      return 0;
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
      return 0;
      break;
      //The FOF additional units - Alan
    case IO_M_MEAN:
    case IO_R_MEAN:
    case IO_M_CRIT:
    case IO_R_CRIT:
    case IO_M_TOPHAT:
    case IO_R_TOPHAT:
      return 0;
      break;
#endif


    case IO_BG_e:
    case IO_BG_HII:
    case IO_BG_Hminus:
    case IO_BG_H2I:
    case IO_BG_H2II:
    case IO_BG_HeII:
    case IO_BG_HeIII:
    case IO_BG_DII:
    case IO_BG_HD:
#ifdef BG_MOL_NETWORK
      return 1;
#else
      return 0;
#endif
      break;

/*     case IO_BG_DUST: */
/* #ifdef BG_DUST */
/*       return 1; */
/* #else */
/*       return 0; */
/* #endif */
/*       break; */

/*     case IO_BG_DUST_SMOOTHED: */
/* #ifdef BG_DUST */
/* #ifdef BG_METALSMOOTHING */
/*       return 1; */
/* #endif */
/* #else */
/*       return 0; */
/* #endif */
/*       break; */

/*     case IO_BG_DUST_SIZE: */
/* #ifdef BG_DUST_DESTRUCTION_SUBLIMATION */
/*       return 1; */
/* #else */
/*       return 0; */
/* #endif */
/*       break; */

/*     case IO_LW_popII: */
/*     case IO_LW_popIII: */
/* #ifdef LW_LOCAL */
/*       return 1; */
/* #endif */
/*       break; */

    case IO_LASTENTRY:
      return 0;			/* will not occur */
      break;
    }

  return 0;			/* default: not present */
}




/*! This function associates a short 4-character block name with each block number.
 *  This is stored in front of each block for snapshot FileFormat=2.
 */
void get_Tab_IO_Label(enum iofields blocknr, char *label)
{
  switch (blocknr)
    {
    case IO_POS:
      strncpy(label, "POS ", 4);
      break;
    case IO_VEL:
      strncpy(label, "VEL ", 4);
      break;
    case IO_ID:
      strncpy(label, "ID  ", 4);
      break;
    case IO_MASS:
      strncpy(label, "MASS", 4);
      break;
    case IO_AO:
      strncpy(label, "AO  ", 4);
      break;
    case IO_U:
      strncpy(label, "U   ", 4);
      break;
    case IO_RHO:
      strncpy(label, "RHO ", 4);
      break;
    case IO_WRHO:
      strncpy(label, "WRHO", 4);
      break;
    case IO_HSML:
      strncpy(label, "HSML", 4);
      break;
    case IO_SFR:
      strncpy(label, "SFR ", 4);
      break;
    case IO_AGE:
      strncpy(label, "AGE ", 4);
      break;
    case IO_Z:
      strncpy(label, "Z   ", 4);
      break;
    case IO_POT:
      strncpy(label, "POT ", 4);
      break;
    case IO_ACCEL:
      strncpy(label, "ACCE", 4);
      break;
    case IO_DTENTR:
      strncpy(label, "ENDT", 4);
      break;
    case IO_TSTP:
      strncpy(label, "TSTP", 4);
      break;
    case IO_BFLD:
      strncpy(label, "BFLD", 4);
      break;
    case IO_DBDT:
      strncpy(label, "DBDT", 4);
      break;
    case IO_DIVB:
      strncpy(label, "DIVB", 4);
      break;
    case IO_ABVC:
      strncpy(label, "ABVC", 4);
      break;
    case IO_AMDC:
      strncpy(label, "AMDC", 4);
      break;
    case IO_PHI:
      strncpy(label, "PHI ", 4);
      break;
    case IO_COOLRATE:
      strncpy(label, "COOR", 4);
      break;
    case IO_CONDRATE:
      strncpy(label, "CONR", 4);
      break;
    case IO_BSMTH:
      strncpy(label, "BFSM", 4);
      break;
    case IO_DENN:
      strncpy(label, "DENN", 4);
      break;
    case IO_EGYPROM:
      strncpy(label, "EGYP", 4);
      break;
    case IO_EGYCOLD:
      strncpy(label, "EGYC", 4);
      break;
    case IO_CR_C0:
      strncpy(label, "CRC0", 4);
      break;
    case IO_CR_Q0:
      strncpy(label, "CRQ0", 4);
      break;
    case IO_CR_P0:
      strncpy(label, "CRP0", 4);
      break;
    case IO_CR_E0:
      strncpy(label, "CRE0", 4);
      break;
    case IO_CR_n0:
      strncpy(label, "CRn0", 4);
      break;
    case IO_CR_ThermalizationTime:
      strncpy(label, "CRco", 4);
      break;
    case IO_CR_DissipationTime:
      strncpy(label, "CRdi", 4);
      break;
    case IO_BHMASS:
      strncpy(label, "BHMA", 4);
      break;
    case IO_BH_ENERGY:
      strncpy(label, "BHEN", 4);
      break;
    case IO_BHMDOT:
      strncpy(label, "BHMD", 4);
      break;
    case IO_BH_BIRTH_TIME:
      strncpy(label, "BHBT", 4);
      break;
    case IO_MACH:
      strncpy(label, "MACH", 4);
      break;
    case IO_DTENERGY:
      strncpy(label, "DTEG", 4);
      break;
    case IO_PRESHOCK_DENSITY:
      strncpy(label, "PSDE", 4);
      break;
    case IO_PRESHOCK_ENERGY:
      strncpy(label, "PSEN", 4);
      break;
    case IO_PRESHOCK_XCR:
      strncpy(label, "PSXC", 4);
      break;
    case IO_DENSITY_JUMP:
      strncpy(label, "DJMP", 4);
      break;
    case IO_ENERGY_JUMP:
      strncpy(label, "EJMP", 4);
      break;
    case IO_CRINJECT:
      strncpy(label, "CRDE", 4);
      break;
    case IO_Zs:
      strncpy(label, "Zs  ", 4);
      break;
    case IO_iMass:
      strncpy(label, "iM  ", 4);
      break;
    case IO_CLDX:
      strncpy(label, "CLDX", 4);
      break;
    case IO_BG_TEMP:
      strncpy(label, "BGT ", 4);
      break;
    case IO_BG_ONEOS:
      strncpy(label, "BEOS", 4);
      break;
    case IO_BG_SNII_KINETIC_FEEDBACK:
      strncpy(label, "BWIN", 4);
      break;
    case IO_BG_METALS:
      strncpy(label, "BGMe", 4);
      break;
    case IO_BG_METALS_SMOOTHED:
      strncpy(label, "BGMS", 4);
      break;
    case IO_BG_INITIAL_MASS:
      strncpy(label, "BGIM", 4);
      break;
    case IO_BG_METALLICITY:
      strncpy(label, "BGZ ", 4);
      break;
    case IO_BG_METALLICITY_SMOOTHED:
      strncpy(label, "BGZS", 4);
      break;
    case IO_BG_IRON_FROM_SNIA:
      strncpy(label, "BGIR", 4);
      break;
    case IO_BG_IRON_FROM_SNIA_SMOOTHED:
      strncpy(label, "BIRS", 4);
      break;
    case IO_BG_MAX_ENTROPY:
      strncpy(label, "BMXE", 4);
      break;
    case IO_BG_MAX_TEMPERATURE:
      strncpy(label, "BMXT", 4);
      break;
    case IO_BG_TIME_MAX_ENTROPY:
      strncpy(label, "TMXE", 4);
      break;
    case IO_BG_TIME_MAX_TEMPERATURE:
      strncpy(label, "TMXT", 4);
      break;
    case IO_BG_METALLICITY_WEIGHTED_REDSHIFT:
      strncpy(label, "MWRZ", 4);
      break;
    case IO_BG_METALLICITY_WEIGHTED_POTENTIAL:
      strncpy(label, "MWPT", 4);
      break;
    case IO_BG_STELLAR_AGE:
      strncpy(label, "BGST", 4);
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
      strncpy(label, "nope", 4);
      break;
      //The FOF additional units - Alan
    case IO_M_MEAN:
    case IO_R_MEAN:
    case IO_M_CRIT:
    case IO_R_CRIT:
    case IO_M_TOPHAT:
    case IO_R_TOPHAT:
      strncpy(label, "nope", 4);
      break;
#endif

    case IO_BG_e:
      strncpy(label, "elec", 4);
      break;
    case IO_BG_HII:
      strncpy(label, "HII ", 4);
      break;
    case IO_BG_Hminus:
      strncpy(label, "Hmin", 4);
      break;
    case IO_BG_H2I:
      strncpy(label, "H2I ", 4);
      break;
    case IO_BG_H2II:
      strncpy(label, "H2II", 4);
      break;
    case IO_BG_HeII:
      strncpy(label, "HeII", 4);
      break;
    case IO_BG_HeIII:
      strncpy(label, "He3I", 4);
      break;
    case IO_BG_DII:
      strncpy(label, "DII ", 4);
      break;
    case IO_BG_HD:
      strncpy(label, "HD  ", 4);
      break;
/*     case IO_BG_DUST: */
/*       strncpy(label, "DUST", 4); */
/*       break; */
/*     case IO_BG_DUST_SMOOTHED: */
/*       strncpy(label, "DUSM", 4); */
/*       break; */
/*     case IO_BG_DUST_SIZE: */
/*       strncpy(label, "DUSS", 4); */
/*       break; */
/*     case IO_LW_popII: */
/*       strncpy(label, "LWp2", 4); */
/*       break; */
/*     case IO_LW_popIII: */
/*       strncpy(label, "LWp3", 4); */
/*       break; */

    case IO_LASTENTRY:
      endrun(217);
      break;
    }
}

float get_h_scaling(enum iofields blocknr)
{
  switch (blocknr)
    {
    case IO_POS:
      return -1;
      break;

    case IO_VEL:
      return 0;
      break;

    case IO_ACCEL:
      return 1;
      break;

    case IO_RHO:
    case IO_WRHO:
      return 2;
      break;

    case IO_HSML:
      return -1;
      break;

    case IO_AO:
      return 0;
      break;

    case IO_MASS:
    case IO_BHMASS:
    case IO_BG_INITIAL_MASS:
      return -1;
      break;

    case IO_BH_ENERGY:
    case IO_BHMDOT:
    case IO_BH_BIRTH_TIME:
      return 0;
      break;

    case IO_U:
      return 0;
      break;

    case IO_BG_MAX_ENTROPY:
      return 2 - 2 * GAMMA;
      break;

    case IO_BG_MAX_TEMPERATURE:
      return 0;
      break;

    case IO_POT:
    case IO_BG_METALLICITY_WEIGHTED_POTENTIAL:
      return 0;
      break;

    case IO_SFR:
      return 0;
      break;

    case IO_BG_TEMP:
      return 0;
      break;

    case IO_BG_STELLAR_AGE:
    case IO_BG_TIME_MAX_ENTROPY:
    case IO_BG_TIME_MAX_TEMPERATURE:
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
      return -1;
      break;

    case IO_SUBFIND_PARTPOS:
    case IO_SUBFIND_PARTVEL:
    case IO_SUBFIND_PARTMASS:
      return -1;
      break;

    case IO_SUBFIND_SPIN:
      return -2;
      break;

    case IO_SUBFIND_VMAX:
      return 0;
      break;

    case IO_SUBFIND_CMVEL:
    case IO_SUBFIND_VELDISP:
    case IO_SUBFIND_STELLARVELDISP:
    case IO_SUBFIND_STELLARVELDISPHALFPROJ:
      return -1;
      break;

      //The FOF additional units - Alan
    case IO_M_MEAN:
    case IO_M_CRIT:
    case IO_M_TOPHAT:
      return -1;
      break;

    case IO_R_MEAN:
    case IO_R_CRIT:
    case IO_R_TOPHAT:
      return -1;
      break;
#endif

    default:
      return 0;
    }
}

double get_units(enum iofields blocknr)
{
  switch (blocknr)
    {
    case IO_POS:
    case IO_HSML:
      return All.UnitLength_in_cm;
      break;

    case IO_VEL:
      return All.UnitVelocity_in_cm_per_s;
      break;

    case IO_ACCEL:
      return All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm;
      break;

    case IO_RHO:
    case IO_WRHO:
      return All.UnitDensity_in_cgs;
      break;

    case IO_AO:
      return 1.0;
      break;


    case IO_MASS:
    case IO_BHMASS:
    case IO_BG_INITIAL_MASS:
      return All.UnitMass_in_g;
      break;

    case IO_U:
      return All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s;
      break;

    case IO_BH_BIRTH_TIME:
      return 1.0;
      break;

    case IO_BHMDOT:
      return All.UnitMass_in_g/All.UnitTime_in_s;
      break;

    case IO_BH_ENERGY:
      return All.UnitEnergy_in_cgs;
      break;
      
    case IO_BG_STELLAR_AGE:
      return SEC_PER_YEAR * 1e9;
      break;

    case IO_BG_MAX_ENTROPY:
      return All.UnitPressure_in_cgs * pow(All.UnitDensity_in_cgs, -GAMMA);
      break;

    case IO_POT:
    case IO_BG_METALLICITY_WEIGHTED_POTENTIAL:
      return All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s;
      break;

    case IO_SFR:
      return SOLAR_MASS / SEC_PER_YEAR;
      break;

#ifdef SUBFIND
      //The SubFind additional units - Alan
    case IO_SUBFIND_MASS:
    case IO_SUBFIND_MASSTYPE:
      return All.UnitMass_in_g;
      break;

    case IO_SUBFIND_POS:
    case IO_SUBFIND_CMPOS:
    case IO_SUBFIND_HALFMASS:
    case IO_SUBFIND_HALFMASSPROJ:
    case IO_SUBFIND_VMAXRAD:
      return All.UnitLength_in_cm;
      break;

    case IO_SUBFIND_SPIN:
      return All.UnitLength_in_cm * All.UnitVelocity_in_cm_per_s;
      break;

    case IO_SUBFIND_VMAX:
    case IO_SUBFIND_CMVEL:
    case IO_SUBFIND_VELDISP:
    case IO_SUBFIND_STELLARVELDISP:
    case IO_SUBFIND_STELLARVELDISPHALFPROJ:
      return All.UnitVelocity_in_cm_per_s;
      break;

    case IO_SUBFIND_PARTPOS:
      return All.UnitLength_in_cm;
      break;

    case IO_SUBFIND_PARTVEL:
      return All.UnitVelocity_in_cm_per_s;
      break;

    case IO_SUBFIND_PARTMASS:
      return All.UnitMass_in_g;
      break;

      //The FOF additional units - Alan
    case IO_M_MEAN:
    case IO_M_CRIT:
    case IO_M_TOPHAT:
      return All.UnitMass_in_g;
      break;

    case IO_R_MEAN:
    case IO_R_CRIT:
    case IO_R_TOPHAT:
      return All.UnitLength_in_cm;
      break;
#endif

    default:
      return 1;
      break;
    }
}

float get_aexp_scaling(enum iofields blocknr)
{
  switch (blocknr)
    {
    case IO_POS:
      return 1.;
      break;

    case IO_VEL:
      return 0.5;
      break;

    case IO_ACCEL:
      return 0;
      break;

    case IO_RHO:
    case IO_WRHO:
      return -3.;
      break;

    case IO_HSML:
      return 1.;
      break;

    case IO_AO:
      return 0.;
      break;

    case IO_MASS:
    case IO_BHMASS:
    case IO_BHMDOT:
    case IO_BG_INITIAL_MASS:
      return 0;
      break;

    case IO_U:
      return -3 * GAMMA_MINUS1;
      break;

    case IO_BH_ENERGY:
      return 0;
      break;

    case IO_BG_MAX_ENTROPY:
      return 0;
      break;

    case IO_BG_TEMP:
    case IO_BG_MAX_TEMPERATURE:
      return 0;
      break;

    case IO_POT:
    case IO_BG_METALLICITY_WEIGHTED_POTENTIAL:
      return -1;
      break;

    case IO_SFR:
      return 0;
      break;

    case IO_BG_STELLAR_AGE:
    case IO_BH_BIRTH_TIME:
    case IO_BG_TIME_MAX_ENTROPY:
    case IO_BG_TIME_MAX_TEMPERATURE:
      return 0;
      break;

#ifdef SUBFIND
      //The SubFind additional units - Alan
    case IO_SUBFIND_MASS:
    case IO_SUBFIND_MASSTYPE:
    case IO_SUBFIND_HALFMASS:
    case IO_SUBFIND_HALFMASSPROJ:
      return 0;
      break;

    case IO_SUBFIND_POS:
    case IO_SUBFIND_CMPOS:
      return 1;
      break;

    case IO_SUBFIND_PARTPOS:
      return 1;
      break;

    case IO_SUBFIND_PARTVEL:
      return -1;
      break;

    case IO_SUBFIND_SPIN:
    case IO_SUBFIND_PARTMASS:
      return 0;
      break;

    case IO_SUBFIND_VMAX:
    case IO_SUBFIND_VELDISP:
    case IO_SUBFIND_STELLARVELDISP:
    case IO_SUBFIND_STELLARVELDISPHALFPROJ:
    case IO_SUBFIND_CMVEL:
    case IO_SUBFIND_VMAXRAD:
      return 0;
      break;

      //The FOF additional units - Alan
    case IO_M_MEAN:
    case IO_M_CRIT:
    case IO_M_TOPHAT:
      return 0;
      break;

    case IO_R_MEAN:
    case IO_R_CRIT:
    case IO_R_TOPHAT:
      return 1;
      break;
#endif

    default:
      return 0;
      break;
    }
}

void get_var_description(enum iofields blocknr, char *buf)
{
  switch (blocknr)
    {
    case IO_POS:
      strcpy(buf, "Co-moving coordinates. Physical: r = a x = Coordinates h^-1 a U_L [cm]");
      break;

    case IO_VEL:
      strcpy(buf, "Co-moving velocities. Physical v_p = a dx/dt  = Velocities a^1/2 U_V [cm/s]");
      break;

    case IO_ACCEL:
      strcpy(buf,
	     "Co-moving accelerations. Physical a_p = a d(dx/dt)/dt + 2 H(a) a dx/dt = Accelerations h U_V^2 / U_L [cm/s^2]");
      break;

    case IO_RHO:
      strcpy(buf, "Co-moving mass densities. Physical rho = Densities h^2 a^-3 U_M U_L^-3 [g/cm^3]");
      break;

    case IO_WRHO:
      strcpy(buf, "Co-moving weighted mass densities. Physical rho = Densities h^2 a^-3 U_M U_L^-3 [g/cm^3]");
      break;

    case IO_HSML:
      strcpy(buf, "Co-moving smoothing length. Physical h = SmoothingLength h^-1 a U_L [cm]");
      break;

    case IO_MASS:
      strcpy(buf, "Particle mass. Physical m = Mass h^-1 U_M [g]");
      break;

    case IO_AO:
      strcpy(buf, "Number of timesteps since last output into an adaptive output file");
      break;

    case IO_U:
      strcpy(buf,
	     "Thermal energy per unit mass. Physical u = InternalEnergy a^(-3(gamma - 1)) U_V^2 [(cm/s)^2]");
      break;

    case IO_AGE:
      strcpy(buf, "Expansion factor a when star particle was born");
      break;

    case IO_BG_ONEOS:
      strcpy(buf,
	     "Flag for effective equation os state. 1 if currently on EoS, 0 if has never been on EoS, -ExpansionFactor if left the EoS at ExpansionFactor");
      break;

    case IO_BG_SNII_KINETIC_FEEDBACK:
      strcpy(buf,
	     "Flag for stellar winds. 0 if has never been a wind particle, -ExpansionFactor if last put in the wind at ExpansionFactor");
      break;

    case IO_BHMASS:
      strcpy(buf, "BH mass. Physical m = Mass h^-1 U_M [g]");
      break;

    case IO_BHMDOT:
      strcpy(buf, "BH accretion rate. Physical mdot = BH_Mdot h^-1 U_M /U_T [g/s]");
      break;

    case IO_BH_BIRTH_TIME:
      strcpy(buf, "BH scale factor of birth");
      break;

    case IO_BH_ENERGY:
      strcpy(buf, "BH energy reservoir.  Physical e = BH_Energy * U_E [cgs]");
      break;

    case IO_ID:
      strcpy(buf, "Unique particle identifier");
      break;

    case IO_POT:
      strcpy(buf, "Co-moving Potential. Physical phi = Potential a^-1 U_V [(cm/s)^2]");
      break;

    case IO_SFR:
      strcpy(buf, "Star formation rate. Physical sfr = StarformationRate SOLAR_MASS SEC_PER_YEAR^-1 [g/s]");
      break;

    case IO_BG_TEMP:
      strcpy(buf, "Temperature [K]");
      break;

    case IO_BG_STELLAR_AGE:
      strcpy(buf, "Stellar age [Gyr]");
      break;

    case IO_BG_MAX_ENTROPY:
      strcpy(buf,
	     "Maximum entropy. Physical entropy = MaximumEntropy h^(2 - 2 gamma) U_P U_rho^-gamma [g (cm/s)^2 (cm^-3)^-gamma]");
      break;

    case IO_BG_MAX_TEMPERATURE:
      strcpy(buf, "Maximum temperature [K]");
      break;

    case IO_BG_TIME_MAX_ENTROPY:
      strcpy(buf, "Expansion factor a when particle had highest entropy");
      break;

    case IO_BG_TIME_MAX_TEMPERATURE:
      strcpy(buf, "Expansion factor a when particle had highest temperature");
      break;

    case IO_BG_METALS:
      strcpy(buf, "Mass fractions of chemical elements");
      break;

    case IO_BG_METALS_SMOOTHED:
      strcpy(buf, "Smoothed mass fractions of chemical elements");
      break;

    case IO_BG_METALLICITY:
      strcpy(buf, "Mass fraction of elements heavier than Helium");
      break;

    case IO_BG_METALLICITY_SMOOTHED:
      strcpy(buf, "Smoothed mass fraction of elements heavier than Helium");
      break;

    case IO_BG_METALLICITY_WEIGHTED_REDSHIFT:
      strcpy(buf, "Metallicity weighted redshift");
      break;

    case IO_BG_METALLICITY_WEIGHTED_POTENTIAL:
      strcpy(buf,
	     "Metallicity weighted co-moving potential. Physical phi = MetallicityWeightedPotential a^-1 U_V [(cm/s)^2]");
      break;

    case IO_BG_INITIAL_MASS:
      strcpy(buf, "Star particle mass at formation time. Physical m = InitialMass h^-1 U_M [g]");
      break;

    case IO_BG_IRON_FROM_SNIA:
      strcpy(buf, "Mass fraction of Iron from SNIa");
      break;

    case IO_BG_IRON_FROM_SNIA_SMOOTHED:
      strcpy(buf, "Smoothed mass fraction of Iron from SNIa");
      break;

#ifdef SUBFIND
      //The SubFind additional units - Alan
    case IO_SUBFIND_MASS:
      strcpy(buf, "Subhalo mass. Physical m = Mass h^-1 U_M [g]");
      break;

    case IO_SUBFIND_CMVEL:
      strcpy(buf,
	     "Subhalo centre of mass velocity. Co-moving velocities. Physical v_p = a dx/dt  = Velocities a^1 U_V [cm/s]");
      break;

    case IO_SUBFIND_VELDISP:
      strcpy(buf,
	     "Subhalo 1D velocity dispersion. Co-moving velocities. Physical v_p = a dx/dt  = Velocities a^1 U_V [cm/s]");
      break;

    case IO_SUBFIND_VMAX:
      strcpy(buf, "Subhalo Maximum circular velocity. Physical v_p = Velocities U_V [cm/s]");
      break;

    case IO_SUBFIND_HALFMASS:
      strcpy(buf, "Subhalo Halfmass radius. Physical r = Length h^-1 U_L [cm]");
      break;

    case IO_SUBFIND_VMAXRAD:
      strcpy(buf,
	     "Subhalo radius at which Maximum circular velocity occurs. Physical r = Length h^-1 U_L [cm]");
      break;

    case IO_SUBFIND_SPIN:
      strcpy(buf,
	     "Subhalo specific angular momentum. Co-moving angular mtm. j = h^-1 a^2 U_L [cm] U_V [cm/s]");
      break;

    case IO_SUBFIND_POS:
      strcpy(buf,
	     "Subhalo Min Potential. Co-moving coordinates. Physical: r = a x = Coordinates h^-1 a U_L [cm]");
      break;

    case IO_SUBFIND_CMPOS:
      strcpy(buf,
	     "Subhalo Centre of Mass. Co-moving coordinates. Physical: r = a x = Coordinates h^-1 a U_L [cm]");
      break;

      //The FOF additional units - Alan
    case IO_M_MEAN:
      strcpy(buf, "Mass within overdensity of 200 times mean background. Physical m = Mass h^-1 U_M [g]");
      break;

    case IO_M_CRIT:
      strcpy(buf, "Mass within overdensity of 200 times critical background. Physical m = Mass h^-1 U_M [g]");
      break;

    case IO_M_TOPHAT:
      strcpy(buf, "Mass within overdensity of Delta_c times crit background. Physical m = Mass h^-1 U_M [g]");
      break;

    case IO_R_MEAN:
      strcpy(buf,
	     "Radius at which the overdensity is 200 times mean background. Co-moving coordinates. Physical: r = a x = Coordinates h^-1 a U_L [cm]");
      break;

    case IO_R_CRIT:
      strcpy(buf,
	     "Radius at which the overdensity is 200 times critical background. Co-moving coordinates. Physical: r = a x = Coordinates h^-1 a U_L [cm]");
      break;

    case IO_R_TOPHAT:
      strcpy(buf,
	     "Radius at which the overdensity is Delta_c times crit background. Co-moving coordinates. Physical: r = a x = Coordinates h^-1 a U_L [cm]");
      break;
#endif

    default:
      strcpy(buf, "");
      break;
    }
}

void get_group_name(enum iofields blocknr, char *buf)
{
  switch (blocknr)
    {
    case IO_BG_METALS:
      strcpy(buf, "ElementAbundance");
      break;

    case IO_BG_METALS_SMOOTHED:
      strcpy(buf, "SmoothedElementAbundance");
      break;

    default:
      strcpy(buf, "default");
      break;
    }
}

void get_dataset_name(enum iofields blocknr, char *buf)
{
  strcpy(buf, "default");

  switch (blocknr)
    {
    case IO_POS:
      strcpy(buf, "Coordinates");
      break;
    case IO_VEL:
      strcpy(buf, "Velocity");
      break;
    case IO_ID:
      strcpy(buf, "ParticleIDs");
      break;
    case IO_MASS:
      strcpy(buf, "Mass");
      break;
    case IO_AO:
      strcpy(buf, "StepsSinceLastOutput");
      break;
    case IO_U:
      strcpy(buf, "InternalEnergy");
      break;
    case IO_RHO:
      strcpy(buf, "Density");
      break;
    case IO_WRHO:
      strcpy(buf, "WeightedDensity");
      break;
    case IO_HSML:
      strcpy(buf, "SmoothingLength");
      break;
    case IO_SFR:
      strcpy(buf, "StarFormationRate");
      break;
    case IO_AGE:
      strcpy(buf, "StellarFormationTime");
      break;
    case IO_Z:
      strcpy(buf, "Metallicity");
      break;
    case IO_POT:
      strcpy(buf, "Potential");
      break;
    case IO_ACCEL:
      strcpy(buf, "Acceleration");
      break;
    case IO_DTENTR:
      strcpy(buf, "RateOfChangeOfEntropy");
      break;
    case IO_TSTP:
      strcpy(buf, "TimeStep");
      break;
    case IO_BFLD:
      strcpy(buf, "MagneticField");
      break;
    case IO_DBDT:
      strcpy(buf, "RateOfChangeOfMagneticField");
      break;
    case IO_DIVB:
      strcpy(buf, "DivergenceOfMagneticField");
      break;
    case IO_ABVC:
      strcpy(buf, "ArtificialViscosity");
      break;
    case IO_AMDC:
      strcpy(buf, "ArtificialMagneticDissipatio");
      break;
    case IO_PHI:
      strcpy(buf, "DivBcleaningFunctionPhi");
      break;
    case IO_COOLRATE:
      strcpy(buf, "CoolingRate");
      break;
    case IO_CONDRATE:
      strcpy(buf, "ConductionRate");
      break;
    case IO_BSMTH:
      strcpy(buf, "SmoothedMagneticField");
      break;
    case IO_DENN:
      strcpy(buf, "Denn");
      break;
    case IO_EGYPROM:
      strcpy(buf, "EnergyReservoirForFeeback");
      break;
    case IO_EGYCOLD:
      strcpy(buf, "EnergyReservoirForColdPhase");
      break;
    case IO_CR_C0:
      strcpy(buf, "CR_C0");
      break;
    case IO_CR_Q0:
      strcpy(buf, "CR_q0");
      break;
    case IO_CR_P0:
      strcpy(buf, "CR_P0");
      break;
    case IO_CR_E0:
      strcpy(buf, "CR_E0");
      break;
    case IO_CR_n0:
      strcpy(buf, "CR_n0");
      break;
    case IO_CR_ThermalizationTime:
      strcpy(buf, "CR_ThermalizationTime");
      break;
    case IO_CR_DissipationTime:
      strcpy(buf, "CR_DissipationTime");
      break;
    case IO_BHMASS:
      strcpy(buf, "BH_Mass");
      break;
    case IO_BH_ENERGY:
      strcpy(buf, "BH_Energy");
      break;
    case IO_BHMDOT:
      strcpy(buf, "BH_Mdot");
      break;
    case IO_BH_BIRTH_TIME:
      strcpy(buf, "BH_BirthTime");
      break;
    case IO_MACH:
      strcpy(buf, "MachNumber");
      break;
    case IO_DTENERGY:
      strcpy(buf, "DtEnergy");
      break;
    case IO_PRESHOCK_DENSITY:
      strcpy(buf, "Preshock_Density");
      break;
    case IO_PRESHOCK_ENERGY:
      strcpy(buf, "Preshock_Energy");
      break;
    case IO_PRESHOCK_XCR:
      strcpy(buf, "Preshock_XCR");
      break;
    case IO_DENSITY_JUMP:
      strcpy(buf, "DensityJump");
      break;
    case IO_ENERGY_JUMP:
      strcpy(buf, "EnergyJump");
      break;
    case IO_CRINJECT:
      strcpy(buf, "CR_DtE");
      break;
    case IO_Zs:
      strcpy(buf, "Zs");
      break;
    case IO_iMass:
      strcpy(buf, "SSPInitialMass");
      break;
    case IO_CLDX:
      strcpy(buf, "CloudFraction");
      break;
    case IO_BG_TEMP:
      strcpy(buf, "Temperature");
      break;
    case IO_BG_ONEOS:
      strcpy(buf, "OnEquationOfState");
      break;
    case IO_BG_SNII_KINETIC_FEEDBACK:
      strcpy(buf, "WindFlag");
      break;
      /* metals */
    case IO_BG_METALS:
      strcpy(buf, "ElementAbundance");
      break;
    case IO_BG_METALS_SMOOTHED:
      strcpy(buf, "SmoothedElementAbundance");
      break;
    case IO_BG_INITIAL_MASS:
      strcpy(buf, "InitialMass");
      break;
    case IO_BG_METALLICITY:
      strcpy(buf, "Metallicity");
      break;
    case IO_BG_METALLICITY_SMOOTHED:
      strcpy(buf, "SmoothedMetallicity");
      break;
    case IO_BG_IRON_FROM_SNIA:
      strcpy(buf, "IronFromSNIa");
      break;
    case IO_BG_IRON_FROM_SNIA_SMOOTHED:
      strcpy(buf, "SmoothedIronFromSNIa");
      break;
    case IO_BG_MAX_ENTROPY:
      strcpy(buf, "MaximumEntropy");
      break;
    case IO_BG_MAX_TEMPERATURE:
      strcpy(buf, "MaximumTemperature");
      break;
    case IO_BG_TIME_MAX_ENTROPY:
      strcpy(buf, "AExpMaximumEntropy");
      break;
    case IO_BG_TIME_MAX_TEMPERATURE:
      strcpy(buf, "AExpMaximumTemperature");
      break;
    case IO_BG_METALLICITY_WEIGHTED_REDSHIFT:
      strcpy(buf, "MetallicityWeightedRedshift");
      break;
    case IO_BG_METALLICITY_WEIGHTED_POTENTIAL:
      strcpy(buf, "MetallicityWeightedPotential");
      break;
    case IO_BG_STELLAR_AGE:
      strcpy(buf, "StellarAgeScaling");
      break;

#ifdef SUBFIND
      //The SubFind additional units - Alan
    case IO_SUBFIND_MASS:
      strcpy(buf, "Mass");
      break;

    case IO_SUBFIND_MASSTYPE:
      strcpy(buf, "MassType");
      break;

    case IO_SUBFIND_CMVEL:
      strcpy(buf, "CenterOfMassVelocity");
      break;

    case IO_SUBFIND_VELDISP:
      strcpy(buf, "SubVelDisp");
      break;

    case IO_SUBFIND_STELLARVELDISP:
      strcpy(buf, "SubStellarVelDisp");
      break;

    case IO_SUBFIND_STELLARVELDISPHALFPROJ:
      strcpy(buf, "SubStellarVelDispHalfProj");
      break;

    case IO_SUBFIND_VMAX:
      strcpy(buf, "SubVmax");
      break;

    case IO_SUBFIND_HALFMASS:
      strcpy(buf, "SubHalfMass");
      break;

    case IO_SUBFIND_HALFMASSPROJ:
      strcpy(buf, "SubHalfMassProj");
      break;

    case IO_SUBFIND_VMAXRAD:
      strcpy(buf, "SubVmaxRad");
      break;

    case IO_SUBFIND_SPIN:
      strcpy(buf, "SubSpin");
      break;

    case IO_SUBFIND_POS:
      strcpy(buf, "Position");
      break;

    case IO_SUBFIND_PARTPOS:
      strcpy(buf, "ParticlePosition");
      break;

    case IO_SUBFIND_PARTVEL:
      strcpy(buf, "ParticleVelocity");
      break;

    case IO_SUBFIND_PARTMASS:
      strcpy(buf, "ParticleMasses");
      break;

    case IO_SUBFIND_CMPOS:
      strcpy(buf, "CenterOfMass");
      break;

      //The FOF additional units - Alan
    case IO_M_MEAN:
      strcpy(buf, "Halo_M_Mean200");
      break;

    case IO_M_CRIT:
      strcpy(buf, "Halo_M_Crit200");
      break;

    case IO_M_TOPHAT:
      strcpy(buf, "Halo_M_TopHat200");
      break;

    case IO_R_MEAN:
      strcpy(buf, "Halo_R_Mean200");
      break;

    case IO_R_CRIT:
      strcpy(buf, "Halo_R_Crit200");
      break;

    case IO_R_TOPHAT:
      strcpy(buf, "Halo_R_TopHat200");
      break;
#endif


    case IO_BG_e:
      strcpy(buf, "ElectronNumberDensity");
      break;
    case IO_BG_HII:
      strcpy(buf, "MassFraction_H+");
      break;
    case IO_BG_Hminus:
      strcpy(buf, "MassFraction_H-");
      break;
    case IO_BG_H2I:
      strcpy(buf, "MassFraction_H2");
      break;
    case IO_BG_H2II:
      strcpy(buf, "MassFraction_H2+");
      break;
    case IO_BG_HeII:
      strcpy(buf, "MassFraction_He+");
      break;
    case IO_BG_HeIII:
      strcpy(buf, "MassFraction_He++");
      break;
    case IO_BG_DII:
      strcpy(buf, "MassFraction_D+");
      break;
    case IO_BG_HD:
      strcpy(buf, "MassFraction_HD");
      break;
/*     case IO_BG_DUST: */
/*       strcpy(buf, "Dusticity"); */
/*       break; */
/*     case IO_BG_DUST_SMOOTHED: */
/*       strcpy(buf, "SmoothedDusticity"); */
/*       break; */
/*     case IO_BG_DUST_SIZE: */
/*       strcpy(buf, "DustGrainSize"); */
/*       break; */
/*     case IO_LW_popII: */
/*       strcpy(buf, "LocalJ21FromPopII"); */
/*       break; */
/*     case IO_LW_popIII: */
/*       strcpy(buf, "LocalJ21FromPopIII"); */
/*       break; */

    case IO_LASTENTRY:
      endrun(218);
      break;
    }
}



/*! This function writes a snapshot file containing the data from processors
 *  'writeTask' to 'lastTask'. 'writeTask' is the one that actually writes.
 *  Each snapshot file contains a header first, then particle positions,
 *  velocities and ID's.  Then particle masses are written for those particle
 *  types with zero entry in MassTable.  After that, first the internal
 *  energies u, and then the density is written for the SPH particles.  If
 *  cooling is enabled, mean molecular weight and neutral hydrogen abundance
 *  are written for the gas particles. This is followed by the SPH smoothing
 *  length and further blocks of information, depending on included physics
 *  and compile-time flags.
 */
void write_file(char *fname, int writeTask, int lastTask)
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

#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0,
    hdf5_dataspace_memory, element_grp[6], selement_grp[6];
#ifdef BG_STELLAR_EVOLUTION
  int k;

  hid_t element_dset[BG_NELEMENTS], hdf5_element_dataspace_in_file[BG_NELEMENTS];

  hid_t selement_dset[BG_NELEMENTS], hdf5_selement_dataspace_in_file[BG_NELEMENTS];

  char mbuf[500];
#endif
  hid_t hdf5_paramgrp = 0, hdf5_unitsgrp = 0, hdf5_constgrp = 0, this_grp;

  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;

  herr_t hdf5_status;

  hsize_t dims[2], count[2], start[2];

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


  /* open file and write header */

  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == 3)
	{
#ifdef HAVE_HDF5
	  sprintf(buf, "%s.hdf5", fname);
	  hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	  hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);
	  hdf5_paramgrp = H5Gcreate(hdf5_file, "/Parameters", 0);
	  hdf5_unitsgrp = H5Gcreate(hdf5_file, "/Units", 0);
	  hdf5_constgrp = H5Gcreate(hdf5_file, "/Constants", 0);

	  for(type = 0; type < 6; type++)
	    {
	      element_grp[type] = -1;
	      selement_grp[type] = -1;
	      hdf5_grp[type] = -1;
	      if(header.npart[type] > 0)
		{
		  sprintf(buf, "/PartType%d", type);
		  hdf5_grp[type] = H5Gcreate(hdf5_file, buf, 0);
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
	  write_header_attributes_in_hdf5(hdf5_headergrp);
	  write_parameters_attributes_in_hdf5(hdf5_paramgrp);
	  write_units_attributes_in_hdf5(hdf5_unitsgrp);
	  write_constants_attributes_in_hdf5(hdf5_constgrp);

#endif
	}
      else
	{
	  if(!(fd = fopen(fname, "w")))
	    {
	      printf("can't open file `%s' for writing snapshot.\n", fname);
	      endrun(123);
	    }

	  if(All.SnapFormat == 2)
	    {
	      blksize = sizeof(int) + 4 * sizeof(char);
	      SKIP;
	      my_fwrite((void *) "HEAD", sizeof(char), 4, fd);
	      nextblock = sizeof(header) + 2 * sizeof(int);
	      my_fwrite(&nextblock, sizeof(int), 1, fd);
	      SKIP;
	    }

	  blksize = sizeof(header);
	  SKIP;
	  my_fwrite(&header, sizeof(header), 1, fd);
	  SKIP;
	}
    }

  ntask = lastTask - writeTask + 1;

  if(ThisTask == writeTask)
    printf("Starting variables output\n");

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

	  npart = get_particles_in_block(blocknr, &typelist[0]);

	  if(npart > 0)
	    {
	      if(ThisTask == writeTask)
		{

		  if(All.SnapFormat == 1 || All.SnapFormat == 2)
		    {
		      if(All.SnapFormat == 2)
			{
			  blksize = sizeof(int) + 4 * sizeof(char);
			  SKIP;
			  get_Tab_IO_Label(blocknr, label);
			  my_fwrite(label, sizeof(char), 4, fd);
			  nextblock = npart * bytes_per_blockelement + 2 * sizeof(int);
			  my_fwrite(&nextblock, sizeof(int), 1, fd);
			  SKIP;
			}

		      blksize = npart * bytes_per_blockelement;
		      SKIP;

		    }
		}

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
			    case 3:
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_INT);
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
				fill_write_buffer(blocknr, &offset, pc, type);

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

	      if(ThisTask == writeTask)
		{
		  if(All.SnapFormat == 1 || All.SnapFormat == 2)
		    SKIP;
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
      if(All.SnapFormat == 3)
	{
#ifdef HAVE_HDF5
	  for(type = 5; type >= 0; type--)
	    if(header.npart[type] > 0)
	      {
#ifdef BG_STELLAR_EVOLUTION
		if(element_grp[type] > 0)
		  H5Gclose(element_grp[type]);

		if(selement_grp[type] > 0)
		  H5Gclose(selement_grp[type]);
#endif
		H5Gclose(hdf5_grp[type]);
	      }
	  H5Gclose(hdf5_headergrp);
	  H5Gclose(hdf5_paramgrp);
	  H5Gclose(hdf5_unitsgrp);
	  H5Gclose(hdf5_constgrp);
	  H5Fclose(hdf5_file);
#endif
	}
      else
	fclose(fd);
    }
  if(ThisTask == writeTask)
    printf("Done variables output\n\n");
}


#ifdef HAVE_HDF5

void write_attributes_in_hdf5(enum iofields blocknr, hid_t dataset)
{
  int ndesc;

  float h_scaling, aexp_scaling;

  double units;

  char var_desc[NDESC];

  hid_t dataspace = 0, attribute = 0, type = 0;

  /* attribute for units */
  units = get_units(blocknr);
  dataspace = H5Screate(H5S_SCALAR);
  attribute = H5Acreate(dataset, "CGSConversionFactor", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_DOUBLE, &units);
  H5Aclose(attribute);
  H5Sclose(dataspace);

  /* attribute for h scaling */
  h_scaling = get_h_scaling(blocknr);
  dataspace = H5Screate(H5S_SCALAR);
  attribute = H5Acreate(dataset, "h-scale-exponent", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_FLOAT, &h_scaling);
  H5Aclose(attribute);
  H5Sclose(dataspace);

  /* attribute for a scaling */
  aexp_scaling = get_aexp_scaling(blocknr);
  dataspace = H5Screate(H5S_SCALAR);
  attribute = H5Acreate(dataset, "aexp-scale-exponent", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_FLOAT, &aexp_scaling);
  H5Aclose(attribute);
  H5Sclose(dataspace);

  /* attribute variable description */
  get_var_description(blocknr, var_desc);
  dataspace = H5Screate(H5S_SCALAR);
  ndesc = NDESC;
  type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, NDESC);
  attribute = H5Acreate(dataset, "VarDescription", type, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, type, var_desc);
  H5Aclose(attribute);
  H5Tclose(type);
  H5Sclose(dataspace);
}


void write_dummy_attributes_in_hdf5(char var_desc[NDESC], hid_t dataset)
{
  int ndesc;

  float fdummy = 0;

  double ddummy = 1;

  hid_t dataspace = 0, attribute = 0, type = 0;

  /* attribute for units */
  dataspace = H5Screate(H5S_SCALAR);
  attribute = H5Acreate(dataset, "CGSConversionFactor", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_DOUBLE, &ddummy);
  H5Aclose(attribute);
  H5Sclose(dataspace);

  /* attribute for h scaling */
  dataspace = H5Screate(H5S_SCALAR);
  attribute = H5Acreate(dataset, "h-scale-exponent", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_FLOAT, &fdummy);
  H5Aclose(attribute);
  H5Sclose(dataspace);

  /* attribute for a scaling */
  dataspace = H5Screate(H5S_SCALAR);
  attribute = H5Acreate(dataset, "aexp-scale-exponent", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_FLOAT, &fdummy);
  H5Aclose(attribute);
  H5Sclose(dataspace);

  /* attribute variable description */
  dataspace = H5Screate(H5S_SCALAR);
  type = H5Tcopy(H5T_C_S1);
  ndesc = NDESC;
  H5Tset_size(type, NDESC);
  attribute = H5Acreate(dataset, "VarDescription", type, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, type, var_desc);
  H5Aclose(attribute);
  H5Tclose(type);
  H5Sclose(dataspace);
}


void write_units_attributes_in_hdf5(hid_t handle)
{
  hid_t hdf5_dataspace, hdf5_attribute;

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "UnitLength_in_cm", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.UnitLength_in_cm);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "UnitMass_in_g", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.UnitMass_in_g);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(handle, "UnitVelocity_in_cm_per_s", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.UnitVelocity_in_cm_per_s);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "UnitDensity_in_cgs", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.UnitDensity_in_cgs);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "UnitEnergy_in_cgs", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.UnitEnergy_in_cgs);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "UnitPressure_in_cgs", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.UnitPressure_in_cgs);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "UnitTime_in_s", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.UnitTime_in_s);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
}


void write_parameters_attributes_in_hdf5(hid_t handle)
{
  hid_t hdf5_dataspace, hdf5_attribute, hdf5_subgrp;

#ifdef BG_STELLAR_EVOLUTION
  hid_t hdf5_datatype;

  hsize_t nelements = BG_NELEMENTS, string_length = STR_LENGTH, el_name_length = EL_NAME_LENGTH;

  int bg_nelements = BG_NELEMENTS;
#endif

  /* Create a group for writing wind paramaters */
#ifdef BG_SNII_KINETIC_FEEDBACK
  hdf5_subgrp = H5Gcreate(handle, "WindParameters", 0);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "SNII_WindOn", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &All.SNII_WindOn);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNII_WindIsotropicOn", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &All.SNII_WindIsotropicOn);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNII_WindSpeed_KMpS", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNII_WindSpeed_KMpS);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNII_WindMassLoading", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNII_WindMassLoading);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNII_WindDelay_YR", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNII_WindDelay_YR);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  H5Gclose(hdf5_subgrp);
#endif

#ifdef BG_STELLAR_EVOLUTION
  /* Create a group for writing elements paramaters */
  hdf5_subgrp = H5Gcreate(handle, "ChemicalElements", 0);

  /* Number of elements */
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "BG_NELEMENTS", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &bg_nelements);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  /* Element names */
  hdf5_datatype = H5Tcopy(H5T_C_S1);
  /*H5Tset_size(hdf5_datatype, H5T_VARIABLE); */
  H5Tset_size(hdf5_datatype, el_name_length);
  hdf5_dataspace = H5Screate_simple(1, &nelements, NULL);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "ElementNames", hdf5_datatype, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, hdf5_datatype, ElementNames);
  H5Aclose(hdf5_attribute);
  H5Tclose(hdf5_datatype);
  H5Sclose(hdf5_dataspace);

  /* initial abundances */
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "InitAbundance_Hydrogen", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.InitAbundance_Hydrogen);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "InitAbundance_Helium", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.InitAbundance_Helium);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "InitAbundance_Carbon", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.InitAbundance_Carbon);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "InitAbundance_Nitrogen", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.InitAbundance_Nitrogen);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "InitAbundance_Oxygen", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.InitAbundance_Oxygen);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "InitAbundance_Neon", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.InitAbundance_Neon);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "InitAbundance_Magnesium", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.InitAbundance_Magnesium);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "InitAbundance_Silicon", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.InitAbundance_Silicon);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "InitAbundance_Iron", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.InitAbundance_Iron);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "CalciumOverSilicon", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.CalciumOverSilicon);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SulphurOverSilicon", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SulphurOverSilicon);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  H5Gclose(hdf5_subgrp);

  /* Create a group for writing stellar evolution paramaters */
  hdf5_subgrp = H5Gcreate(handle, "StellarEvolutionParameters", 0);

  /* star formation parameters */
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SF_EOSGammaEffective", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SF_EOSGammaEffective);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SF_EOSEnergyAtThreshold_ERG", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SF_EOSEnergyAtThreshold_ERG);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SF_EOSMinPhysDens_HpCM3", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SF_EOSMinPhysDens_HpCM3);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SF_THRESH_MinPhysDens_HpCM3", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SF_THRESH_MinPhysDens_HpCM3);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SF_THRESH_MaxTemp_K", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SF_THRESH_MaxTemp_K);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SF_THRESH_MinOverDens", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SF_THRESH_MinOverDens);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SF_SchmidtLawCoeff_MSUNpYRpKPC2", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SF_SchmidtLawCoeff_MSUNpYRpKPC2);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SF_SchmidtLawExponent", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SF_SchmidtLawExponent);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SF_SchmidtLawCoeff_GpSpCM2", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SF_SchmidtLawCoeff_GpSpCM2);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  /* IMF parameters */
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(hdf5_datatype, string_length);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "IMF_Model", hdf5_datatype, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, hdf5_datatype, All.IMF_Model);
  H5Aclose(hdf5_attribute);
  H5Tclose(hdf5_datatype);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(hdf5_datatype, string_length);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "IMF_LifetimeModel", hdf5_datatype, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, hdf5_datatype, All.IMF_LifetimeModel);
  H5Aclose(hdf5_attribute);
  H5Tclose(hdf5_datatype);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "IMF_MinMass_MSUN", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.IMF_MinMass_MSUN);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "IMF_MaxMass_MSUN", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.IMF_MaxMass_MSUN);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  /* SNIa parameters */
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(hdf5_datatype, string_length);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "SNIa_Model", hdf5_datatype, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, hdf5_datatype, All.SNIa_Model);
  H5Aclose(hdf5_attribute);
  H5Tclose(hdf5_datatype);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNIa_Efficiency_fracwd", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNIa_Efficiency_fracwd);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "SNIa_Energy_ERG", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNIa_Energy_ERG);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "SNIa_MassTransferOn", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &All.SNIa_MassTransferOn);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNIa_EnergyTransferOn", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &All.SNIa_EnergyTransferOn);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  /* SNII parameters */
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "SNII_MassTransferOn", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &All.SNII_MassTransferOn);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNII_MinMass_MSUN", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNII_MinMass_MSUN);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNII_MaxMass_MSUN", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNII_MaxMass_MSUN);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNII_Factor_Hydrogen", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNII_Factor_Hydrogen);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNII_Factor_Helium", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNII_Factor_Helium);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNII_Factor_Carbon", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNII_Factor_Carbon);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNII_Factor_Nitrogen", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNII_Factor_Nitrogen);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNII_Factor_Oxygen", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNII_Factor_Oxygen);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "SNII_Factor_Neon", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNII_Factor_Neon);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNII_Factor_Magnesium", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNII_Factor_Magnesium);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "SNII_Factor_Silicon", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNII_Factor_Silicon);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "SNII_Factor_Iron", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNII_Factor_Iron);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  /* AGB parameters */
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "AGB_MassTransferOn", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &All.AGB_MassTransferOn);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  H5Gclose(hdf5_subgrp);
#endif

  /* Numerical parameters */
  hdf5_subgrp = H5Gcreate(handle, "NumericalParameters", 0);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "ComovingIntegrationOn", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &All.ComovingIntegrationOn);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "TypeOfTimestepCriterion", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &All.TypeOfTimestepCriterion);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "ErrTolIntAccuracy", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.ErrTolIntAccuracy);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "MaxSizeTimestep", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.MaxSizeTimestep);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "MinSizeTimestep", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.MinSizeTimestep);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "ErrTolTheta", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.ErrTolTheta);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "TypeOfOpeningCriterion", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &All.TypeOfOpeningCriterion);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "ErrTolForceAcc", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.ErrTolForceAcc);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "TreeDomainUpdateFrequency", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.TreeDomainUpdateFrequency);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "DesNumNgb", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &All.DesNumNgb);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "MaxNumNgbDeviation", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.MaxNumNgbDeviation);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "ArtBulkViscConst", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.ArtBulkViscConst);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
#ifdef BG_SFR
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "InitGasU_ERG", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.InitGasU_ERG);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "MinGasU_ERG", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.MinGasU_ERG);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
#endif
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "CourantFac", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.CourantFac);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "PartAllocFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.PartAllocFactor);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "TreeAllocFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.TreeAllocFactor);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_subgrp, "BufferSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BufferSize);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(hdf5_subgrp, "MinGasHsmlFractional", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.MinGasHsmlFractional);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  H5Gclose(hdf5_subgrp);
}


void write_constants_attributes_in_hdf5(hid_t handle)
{
  hid_t hdf5_dataspace, hdf5_attribute;

  double temp;

  temp = PI;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "PI", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = GAMMA;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "GAMMA", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = GRAVITY;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "GRAVITY", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = SOLAR_MASS;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "SOLAR_MASS", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = SOLAR_LUM;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "SOLAR_LUM", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = RAD_CONST;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "RAD_CONST", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = AVOGADRO;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "AVOGADRO", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = BOLTZMANN;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "BOLTZMANN", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = GAS_CONST;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "GAS_CONST", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = C;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "C", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = PLANCK;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "PLANCK", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = CM_PER_MPC;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "CM_PER_MPC", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = PROTONMASS;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "PROTONMASS", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = ELECTRONMASS;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "ELECTRONMASS", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = ELECTRONCHARGE;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "ELECTRONCHARGE", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = HUBBLE;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "HUBBLE", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = T_CMB0;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "T_CMB0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = SEC_PER_MEGAYEAR;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "SEC_PER_MEGAYEAR", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  temp = SEC_PER_YEAR;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "SEC_PER_YEAR", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
#ifdef BG_COOLING
  temp = STEFAN;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "STEFAN", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
#endif
  temp = THOMPSON;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "THOMPSON", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
#ifdef BG_COOLING
  temp = EV_TO_ERG;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "EV_TO_ERG", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
#endif
  temp = ZSOLAR;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Z_Solar", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &temp);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

}


void write_header_attributes_in_hdf5(hid_t handle)
{
  hsize_t adim[1] = { 6 };
  hid_t hdf5_dataspace, hdf5_attribute, hdf5_type;
 
  herr_t hdf5_status;

  /* RunLabel */
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_type = H5Tcopy(H5T_C_S1);
  hdf5_status = H5Tset_size(hdf5_type, 200);
  hdf5_attribute = H5Acreate(handle, "RunLabel", hdf5_type, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, hdf5_type, All.RunLabel);
  H5Aclose(hdf5_attribute);
  H5Tclose(hdf5_type);
  H5Sclose(hdf5_dataspace);
 
  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npart);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total_HighWord", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);


  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "ExpansionFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  double time_in_gyr;

  time_in_gyr = bg_get_elapsed_time(0, header.time, 1);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Time_GYR", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &time_in_gyr);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  /* this is not in the "standard" header!!!!!!!! */
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "OmegaBaryon", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.OmegaBaryon);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
}
#endif





/*! This catches I/O errors occuring for my_fwrite(). In this case we
 *  better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fwrite) on task=%d has occured: %s\n", ThisTask, strerror(errno));
      fflush(stdout);
      endrun(777);
    }
  return nwritten;
}


/*! This catches I/O errors occuring for fread(). In this case we
 *  better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      if(feof(stream))
	printf("I/O error (fread) on task=%d has occured: end of file\n", ThisTask);
      else
	printf("I/O error (fread) on task=%d has occured: %s\n", ThisTask, strerror(errno));
      printf("I/O error (fread) on task=%d has occured: %s\n", ThisTask, strerror(errno));
      fflush(stdout);
      endrun(778);
    }
  return nread;
}
