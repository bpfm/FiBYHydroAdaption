#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

#include "allvars.h"
#include "proto.h"
#include "bg_proto.h"

/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */

/*! This routine contains the main simulation loop that iterates over
 * single timesteps. The loop terminates when the cpu-time limit is
 * reached, when a `stop' file is found in the output directory, or
 * when the simulation ends because we arrived at TimeMax.
 */
void run(void)
{
  FILE *fd;

  int stopflag = 0;

  char buf[200], stopfname[200], contfname[200];


  sprintf(stopfname, "%sstop", All.OutputDir);
  sprintf(contfname, "%scont", All.OutputDir);
  unlink(contfname);


  CPU_Step[CPU_MISC] += measure_time();


  do				/* main loop */
    {

      find_next_sync_point_and_drift();	/* find next synchronization point and drift particles to this time.
					 * If needed, this function will also write an output file
					 * at the desired time.
					 */

      every_timestep_stuff();	/* write some info to log-files */

#ifdef COMPUTE_POTENTIAL_ENERGY
      if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics)
	All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TotNumPart * All.TreeDomainUpdateFrequency);
#endif

      domain_Decomposition();	/* do domain decomposition if needed */


      compute_accelerations(0);	/* compute accelerations for 
				 * the particles that are to be advanced  
				 */

      /* check whether we want a full energy statistics */
      if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics)
	{
#ifdef COMPUTE_POTENTIAL_ENERGY
	  compute_potential();
#endif
	  energy_statistics();	/* compute and output energy statistics */

	  All.TimeLastStatistics += All.TimeBetStatistics;
	}

#ifdef BG_EXTRA_ARRAYS
      bg_compute_extra_arrays();
#endif

      advance_and_find_timesteps();	/* 'kick' active particles in
					 * momentum space and compute new
					 * timesteps for them
					 */

      write_cpu_log();		/* produce some CPU usage info */

      All.NumCurrentTiStep++;

      /* Check whether we need to interrupt the run */
      if(ThisTask == 0)
	{
	  /* Is the stop-file present? If yes, interrupt the run. */
	  if((fd = fopen(stopfname, "r")))
	    {
	      fclose(fd);
	      stopflag = 1;
	      unlink(stopfname);
	    }

	  /* are we running out of CPU-time ? If yes, interrupt run. */
	  if(CPUThisRun > 0.85 * All.TimeLimitCPU)
	    {
	      printf("CPUThisRun = %g, TimeLimitCPU = %g\n", CPUThisRun, All.TimeLimitCPU);
	      printf("reaching time-limit. stopping.\n");
	      stopflag = 2;
	    }
	}

      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(stopflag)
	{
#ifdef BG_SNAPSHOT_AS_RESTART_FILE
	  savepositions(All.SnapshotFileCount - 1, stopflag);
#else
	  restart(0);		/* write restart file */
#endif
	  MPI_Barrier(MPI_COMM_WORLD);

	  if(stopflag == 2 && ThisTask == 0)
	    {
	      if((fd = fopen(contfname, "w")))
		fclose(fd);
	    }

	  if(stopflag == 2 && All.ResubmitOn && ThisTask == 0)
	    {
	      close_outputfiles();
	      sprintf(buf, "%s", All.ResubmitCommand);
#ifndef NOCALLSOFSYSTEM
	      int ret;

	      ret = system(buf);
#endif
	    }
	  return;
	}

      /* is it time to write a regular restart-file? (for security) */
      if(ThisTask == 0)
	{
	  if((CPUThisRun - All.TimeLastRestartFile) >= All.CpuTimeBetRestartFile)
	    {
	      All.TimeLastRestartFile = CPUThisRun;
	      stopflag = 3;
	    }
	  else
	    stopflag = 0;
	}

      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(stopflag == 3)
	{
#ifdef BG_SNAPSHOT_AS_RESTART_FILE
	  savepositions(All.SnapshotFileCount - 1, stopflag);
#else
	  restart(0);		/* write an occasional restart file */
#endif
	  stopflag = 0;
	}
    }
  while(All.Ti_Current < TIMEBASE && All.Time <= All.TimeMax);

#ifndef BG_SNAPSHOT_AS_RESTART_FILE
  restart(0);
#endif

  savepositions(All.SnapshotFileCount++, 0);	/* write a last snapshot
						 * file at final time (will
						 * be overwritten if
						 * All.TimeMax is increased
						 * and the run is continued)
						 */
}


/*! This function finds the next synchronization point of the system
 * (i.e. the earliest point of time any of the particles needs a force
 * computation), and drifts the system to this point of time.  If the
 * system dirfts over the desired time of a snapshot file, the
 * function will drift to this moment, generate an output, and then
 * resume the drift.
 */
void find_next_sync_point_and_drift(void)
{
  int n, i, prev, dt_bin, ti_next_for_bin, ti_next_kick, ti_next_kick_global;

  long long numforces2;

  double timeold;

  timeold = All.Time;

  /* find the next kick time */
  for(n = 0, ti_next_kick = TIMEBASE; n < TIMEBINS; n++)
    {
      if(TimeBinCount[n])
	{
	  if(n > 0)
	    {
	      dt_bin = (1 << n);
	      ti_next_for_bin = (All.Ti_Current / dt_bin) * dt_bin + dt_bin;	/* next kick time for this timebin */
	    }
	  else
	    {
	      dt_bin = 0;
	      ti_next_for_bin = All.Ti_Current;
	    }

	  if(ti_next_for_bin < ti_next_kick)
	    ti_next_kick = ti_next_for_bin;
	}
    }

  MPI_Allreduce(&ti_next_kick, &ti_next_kick_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  /* mark the bins that will be active */
  for(n = 1, TimeBinActive[0] = 1, NumForceUpdate = TimeBinCount[0]; n < TIMEBINS; n++)
    {
      dt_bin = (1 << n);

      if((ti_next_kick_global % dt_bin) == 0)
	{
	  TimeBinActive[n] = 1;
	  NumForceUpdate += TimeBinCount[n];
	}
      else
	TimeBinActive[n] = 0;
    }

  sumup_large_ints(1, &NumForceUpdate, &GlobNumForceUpdate);

  if(GlobNumForceUpdate == All.TotNumPart)
    Flag_FullStep = 1;
  else
    Flag_FullStep = 0;

  All.NumForcesSinceLastDomainDecomp += GlobNumForceUpdate;

#if defined(BG_THERMAL_FEEDBACK_MAKE_PARTICLES_ACTIVE) && defined(TIMESTEP_LIMITER)
  int lowest_occupied_bin;

  for(n = TIMEBINS - 1; n >= 0 ; n--)
    if(TimeBinCount[n])
      lowest_occupied_bin = n;

  MPI_Allreduce(&lowest_occupied_bin, &All.LowestOccupiedTimeBin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif

  while((ti_next_kick_global >= All.Ti_nextoutput && All.Ti_nextoutput >= 0) ||
#ifdef BG_OUTPUT_GRID
	(ti_next_kick_global >= All.Ti_nextgridoutput && All.Ti_nextgridoutput >= 0) ||
#endif
#ifdef FOF
	(ti_next_kick_global >= All.Ti_nextfof && All.Ti_nextfof >= 0) ||
#endif
	(ti_next_kick_global >= All.Ti_nextlineofsight && All.Ti_nextlineofsight > 0))
    {
      if(All.Ti_nextlineofsight > 0 &&
#ifdef BG_OUTPUT_GRID
	 All.Ti_nextlineofsight < All.Ti_nextgridoutput &&
#endif
#ifdef FOF
	 All.Ti_nextlineofsight < All.Ti_nextfof &&
#endif
	 All.Ti_nextlineofsight < All.Ti_nextoutput)
	{
	  move_particles(All.Ti_nextlineofsight);
	  All.Ti_Current = All.Ti_nextlineofsight;
	}
#ifdef BG_OUTPUT_GRID
      else if(All.Ti_nextgridoutput >= 0 &&
#ifdef FOF
	      All.Ti_nextgridoutput < All.Ti_nextfof &&
#endif
	      All.Ti_nextgridoutput < All.Ti_nextoutput)
	{
	  move_particles(All.Ti_nextgridoutput);
	  All.Ti_Current = All.Ti_nextgridoutput;
	}
#endif
#ifdef FOF
      else if(All.Ti_nextfof >= 0 && All.Ti_nextfof < All.Ti_nextoutput)
	{
	  move_particles(All.Ti_nextfof);
	  All.Ti_Current = All.Ti_nextfof;
	}
#endif
      else
	{
	  move_particles(All.Ti_nextoutput);
	  All.Ti_Current = All.Ti_nextoutput;
	}

      if(All.ComovingIntegrationOn)
	All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
      else
	All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;


      if(All.Ti_Current == All.Ti_nextoutput)
	{
	  CPU_Step[CPU_DRIFT] += measure_time();

#ifdef OUTPUTPOTENTIAL
#if !defined(EVALPOTENTIAL) || (defined(EVALPOTENTIAL) && defined(RECOMPUTE_POTENTIAL_ON_OUTPUT))
	  All.NumForcesSinceLastDomainDecomp =
	    (long long) (1 + All.TotNumPart * All.TreeDomainUpdateFrequency);
	  domain_Decomposition();
	  compute_potential();
#endif
#endif

	  savepositions(All.SnapshotFileCount++, 0);	/* write snapshot file */

	  All.Ti_nextoutput = find_next_outputtime(All.Ti_nextoutput + 1);
	}

#ifdef BG_OUTPUT_GRID
      if(All.Ti_Current == All.Ti_nextgridoutput)
	{
	  CPU_Step[CPU_DRIFT] += measure_time();

	  binning();		/* write grid file */
	  All.Ti_nextgridoutput = find_next_gridoutputtime(All.Ti_nextgridoutput + 1);

	  CPU_Step[CPU_BINNING] += measure_time();
	}
#endif

#ifdef OUTPUTLINEOFSIGHT
      if(All.Ti_Current == All.Ti_nextlineofsight)
	{
	  CPU_Step[CPU_DRIFT] += measure_time();

	  lineofsight_output();

	  All.Ti_nextlineofsight = find_next_losoutputtime(All.Ti_nextlineofsight + 1);

	  CPU_Step[CPU_LINEOFSIGHT] += measure_time();
	}
#endif

#ifdef FOF
      if(All.Ti_Current == All.Ti_nextfof)
	{
	  CPU_Step[CPU_DRIFT] += measure_time();

	  fof_fof(All.FOFFileCount++);

	  All.Ti_nextfof = find_next_fofoutputtime(All.Ti_nextfof + 1);

	  CPU_Step[CPU_FOF] += measure_time();
	}
#endif
    }

  All.Ti_Current = ti_next_kick_global;

  if(All.ComovingIntegrationOn)
    All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
  else
    All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

  All.TimeStep = All.Time - timeold;


  FirstActiveParticle = -1;

  for(n = 0, prev = -1; n < TIMEBINS; n++)
    {
      if(TimeBinActive[n])
	{
	  for(i = FirstInTimeBin[n]; i >= 0; i = NextInTimeBin[i])
	    {
	      if(prev == -1)
		FirstActiveParticle = i;

	      if(prev >= 0)
		NextActiveParticle[prev] = i;

	      prev = i;
	    }
	}
    }

  if(prev >= 0)
    NextActiveParticle[prev] = -1;


  /* drift the active particles, others will be drifted on the fly if needed */

  for(i = FirstActiveParticle, NumForceUpdate = 0; i >= 0; i = NextActiveParticle[i])
    {
      drift_particle(i, All.Ti_Current);

      NumForceUpdate++;
    }

  sumup_large_ints(1, &NumForceUpdate, &numforces2);
  if(GlobNumForceUpdate != numforces2)
    {
      printf("terrible\n");
      endrun(2);
    }


  CPU_Step[CPU_DRIFT] += measure_time();
}




int ShouldWeDoDynamicUpdate(void)
{
  int n, num, dt_bin, ti_next_for_bin, ti_next_kick, ti_next_kick_global;

  long long numforces;


  /* find the next kick time */
  for(n = 0, ti_next_kick = TIMEBASE; n < TIMEBINS; n++)
    {
      if(TimeBinCount[n])
	{
	  if(n > 0)
	    {
	      dt_bin = (1 << n);
	      ti_next_for_bin = (All.Ti_Current / dt_bin) * dt_bin + dt_bin;	/* next kick time for this timebin */
	    }
	  else
	    {
	      dt_bin = 0;
	      ti_next_for_bin = All.Ti_Current;
	    }

	  if(ti_next_for_bin < ti_next_kick)
	    ti_next_kick = ti_next_for_bin;
	}
    }

  MPI_Allreduce(&ti_next_kick, &ti_next_kick_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  /* count the particles that will be active */
  for(n = 1, num = TimeBinCount[0]; n < TIMEBINS; n++)
    {
      dt_bin = (1 << n);

      if((ti_next_kick_global % dt_bin) == 0)
	num += TimeBinCount[n];
    }

  sumup_large_ints(1, &num, &numforces);

  if(ThisTask == 0)
    printf("I'm guessing %d%09d particles to be active in the next step\n",
	   (int) (numforces / 1000000000), (int) (numforces % 1000000000));

  if((All.NumForcesSinceLastDomainDecomp + numforces) >= All.TreeDomainUpdateFrequency * All.TotNumPart)
    return 0;
  else
    return 1;
}



/*! this function returns the next output time that is equal or larger to
 *  ti_curr
 */
int find_next_outputtime(int ti_curr)
{
  int i, ti, ti_next, iter = 0;

  double next, time;

  ti_next = -1;

  if(All.OutputListOn)
    {
      for(i = 0; i < All.OutputListLength; i++)
	{
	  time = All.OutputListTimes[i];

	  if(time >= All.TimeBegin && time <= All.TimeMax)
	    {
	      if(All.ComovingIntegrationOn)
		ti = (int) (log(time / All.TimeBegin) / All.Timebase_interval);
	      else
		ti = (int) ((time - All.TimeBegin) / All.Timebase_interval);

	      if(ti >= ti_curr)
		{
		  if(ti_next == -1)
		    ti_next = ti;

		  if(ti_next > ti)
		    ti_next = ti;
		}
	    }
	}
    }
  else
    {
      if(All.ComovingIntegrationOn)
	{
	  if(All.TimeBetSnapshot <= 1.0)
	    {
	      printf("TimeBetSnapshot > 1.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}
      else
	{
	  if(All.TimeBetSnapshot <= 0.0)
	    {
	      printf("TimeBetSnapshot > 0.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}

      time = All.TimeOfFirstSnapshot;

      iter = 0;

      while(time < All.TimeBegin)
	{
	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(110);
	    }
	}

      while(time <= All.TimeMax)
	{
	  if(All.ComovingIntegrationOn)
	    ti = (int) (log(time / All.TimeBegin) / All.Timebase_interval);
	  else
	    ti = (int) ((time - All.TimeBegin) / All.Timebase_interval);

	  if(ti >= ti_curr)
	    {
	      ti_next = ti;
	      break;
	    }

	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(111);
	    }
	}
    }

  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE;	/* this will prevent any further output */

      if(ThisTask == 0)
	printf("\nThere is no valid time for a further snapshot file.\n");
    }
  else
    {
      if(All.ComovingIntegrationOn)
	next = All.TimeBegin * exp(ti_next * All.Timebase_interval);
      else
	next = All.TimeBegin + ti_next * All.Timebase_interval;

      if(ThisTask == 0)
	printf("\nSetting next time for snapshot file to Time_next= %g\n\n", next);
    }

  return ti_next;
}


#ifdef BG_OUTPUT_GRID
int find_next_gridoutputtime(int ti_curr)
{
  int ti, ti_next, iter = 0;

  double next, time;

  ti_next = -1;

  if(All.ComovingIntegrationOn)
    {
      if(All.TimeBetGridOutput <= 1.0)
	{
	  printf("TimeBetGridOutput > 1.0 required for your simulation.\n");
	  endrun(13123);
	}
    }
  else
    {
      if(All.TimeBetGridOutput <= 0.0)
	{
	  printf("TimeBetGridOutput > 0.0 required for your simulation.\n");
	  endrun(13123);
	}
    }

  time = All.TimeOfFirstGridOutput;

  iter = 0;

  while(time < All.TimeBegin)
    {
      if(All.ComovingIntegrationOn)
	time *= All.TimeBetGridOutput;
      else
	time += All.TimeBetGridOutput;

      iter++;

      if(iter > 1000000)
	{
	  printf("Can't determine next grid output time.\n");
	  endrun(110);
	}
    }

  while(time <= All.TimeMax)
    {
      if(All.ComovingIntegrationOn)
	ti = (int) (log(time / All.TimeBegin) / All.Timebase_interval);
      else
	ti = (int) ((time - All.TimeBegin) / All.Timebase_interval);

      if(ti >= ti_curr)
	{
	  ti_next = ti;
	  break;
	}

      if(All.ComovingIntegrationOn)
	time *= All.TimeBetGridOutput;
      else
	time += All.TimeBetGridOutput;

      iter++;

      if(iter > 1000000)
	{
	  printf("Can't determine next grid output time.\n");
	  endrun(111);
	}
    }

  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE;	/* this will prevent any further output */

      if(ThisTask == 0)
	printf("\nThere is no valid time for a further grid file.\n");
    }
  else
    {
      if(All.ComovingIntegrationOn)
	next = All.TimeBegin * exp(ti_next * All.Timebase_interval);
      else
	next = All.TimeBegin + ti_next * All.Timebase_interval;

      if(ThisTask == 0)
	printf("\nSetting next time for grid output file to Time_next= %g\n\n", next);
    }

  return ti_next;
}
#endif


#ifdef OUTPUTLINEOFSIGHT
int find_next_losoutputtime(int ti_curr)
{
  int ti, ti_next, iter = 0;

  double next, time;

  ti_next = -1;

  if(All.ComovingIntegrationOn)
    {
      if(All.TimeBetLineOfSight <= 1.0)
	{
	  printf("TimeBetLineOfSight > 1.0 required for your simulation.\n");
	  endrun(13123);
	}
    }
  else
    {
      if(All.TimeBetLineOfSight <= 0.0)
	{
	  printf("TimeBetLineOfSight > 0.0 required for your simulation.\n");
	  endrun(13123);
	}
    }

  time = All.TimeOfFirstLineOfSight;

  iter = 0;

  while(time < All.TimeBegin)
    {
      if(All.ComovingIntegrationOn)
	time *= All.TimeBetLineOfSight;
      else
	time += All.TimeBetLineOfSight;

      iter++;

      if(iter > 1000000)
	{
	  printf("Can't determine next line of sight output time.\n");
	  endrun(110);
	}
    }

  while(time <= All.TimeMax)
    {
      if(All.ComovingIntegrationOn)
	ti = (int) (log(time / All.TimeBegin) / All.Timebase_interval);
      else
	ti = (int) ((time - All.TimeBegin) / All.Timebase_interval);

      if(ti >= ti_curr)
	{
	  ti_next = ti;
	  break;
	}

      if(All.ComovingIntegrationOn)
	time *= All.TimeBetLineOfSight;
      else
	time += All.TimeBetLineOfSight;

      iter++;

      if(iter > 1000000)
	{
	  printf("Can't determine next line of sight output time.\n");
	  endrun(111);
	}
    }

  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE;	/* this will prevent any further output */

      if(ThisTask == 0)
	printf("\nThere is no valid time for a further line of sight file.\n");
    }
  else
    {
      if(All.ComovingIntegrationOn)
	next = All.TimeBegin * exp(ti_next * All.Timebase_interval);
      else
	next = All.TimeBegin + ti_next * All.Timebase_interval;

      if(ThisTask == 0)
	printf("\nSetting next time for line of sight file to Time_next= %g\n\n", next);
    }

  return ti_next;
}
#endif



#ifdef FOF
int find_next_fofoutputtime(int ti_curr)
{
  int ti, ti_next, iter = 0;

  double next, time;

  ti_next = -1;

  if(All.ComovingIntegrationOn)
    {
      if(All.TimeBetFOF <= 1.0)
	{
	  printf("TimeBetFOF > 1.0 required for your simulation.\n");
	  endrun(13123);
	}
    }
  else
    {
      if(All.TimeBetFOF <= 0.0)
	{
	  printf("TimeBetFOF > 0.0 required for your simulation.\n");
	  endrun(13123);
	}
    }

  time = All.TimeOfFirstFOF;

  iter = 0;

  while(time < All.TimeBegin)
    {
      if(All.ComovingIntegrationOn)
	time *= All.TimeBetFOF;
      else
	time += All.TimeBetFOF;

      iter++;

      if(iter > 1000000)
	{
	  printf("Can't determine next FOF output time.\n");
	  endrun(110);
	}
    }

  while(time <= All.TimeMax)
    {
      if(All.ComovingIntegrationOn)
	ti = (int) (log(time / All.TimeBegin) / All.Timebase_interval);
      else
	ti = (int) ((time - All.TimeBegin) / All.Timebase_interval);

      if(ti >= ti_curr)
	{
	  ti_next = ti;
	  break;
	}

      if(All.ComovingIntegrationOn)
	time *= All.TimeBetFOF;
      else
	time += All.TimeBetFOF;

      iter++;

      if(iter > 1000000)
	{
	  printf("Can't determine next FOF output time.\n");
	  endrun(111);
	}
    }

  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE;	/* this will prevent any further output */

      if(ThisTask == 0)
	printf("\nThere is no valid time for a further FOF file.\n");
    }
  else
    {
      if(All.ComovingIntegrationOn)
	next = All.TimeBegin * exp(ti_next * All.Timebase_interval);
      else
	next = All.TimeBegin + ti_next * All.Timebase_interval;

      if(ThisTask == 0)
	printf("\nSetting next time for FOF file to Time_next= %g\n\n", next);
    }

  return ti_next;
}
#endif




/*! This routine writes one line for every timestep to two log-files.
 * In FdInfo, we just list the timesteps that have been done, while in
 * FdCPU the cumulative cpu-time consumption in various parts of the
 * code is stored.
 *
 * Additionally, if adaptive timestepping has been requested then this
 * is carried out.
 */
void every_timestep_stuff(void)
{
  double z;

  int i, j;
  long long tot, tot_sph;
  long long tot_count[TIMEBINS], tot_count_sph[TIMEBINS];
  int *temp;

  /*
  MPI_Reduce(TimeBinCount, tot_count, TIMEBINS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(TimeBinCountSph, tot_count_sph, TIMEBINS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  */

  temp = (int *) mymalloc(NTask * TIMEBINS * sizeof(int));

  MPI_Allgather(TimeBinCount, TIMEBINS, MPI_INT, temp, TIMEBINS, MPI_INT, MPI_COMM_WORLD);
  for(i = 0; i < TIMEBINS; i++)
    {
      tot_count[i] = 0;
      for(j = 0; j < NTask; j++)
        tot_count[i] += temp[j * TIMEBINS + i];
    }

  MPI_Allgather(TimeBinCountSph, TIMEBINS, MPI_INT, temp, TIMEBINS, MPI_INT, MPI_COMM_WORLD);
  for(i = 0; i < TIMEBINS; i++)
    {
      tot_count_sph[i] = 0;
      for(j = 0; j < NTask; j++)
        tot_count_sph[i] += temp[j * TIMEBINS + i];
    }

  myfree(temp);



  if(ThisTask == 0)
    {
      if(All.ComovingIntegrationOn)
	{
	  z = 1.0 / (All.Time) - 1;
	  fprintf(FdInfo, "\nBegin Step %d, Time: %g, Redshift: %g, Nf = %d%09d, Systemstep: %g, Dloga: %g\n",
		  All.NumCurrentTiStep, All.Time, z,
		  (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000),
		  All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
	  printf("\nBegin Step %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", All.NumCurrentTiStep,
		 All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
	  fflush(FdInfo);
	}
      else
	{
	  fprintf(FdInfo, "\nBegin Step %d, Time: %g, Nf = %d%09d, Systemstep: %g\n", All.NumCurrentTiStep,
		  All.Time, (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000),
		  All.TimeStep);
	  printf("\nBegin Step %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time, All.TimeStep);
	  fflush(FdInfo);
	}

      printf("Occupied timebins: non-sph         sph       dt\n");
      for(i = TIMEBINS - 1, tot = tot_sph = 0; i >= 0; i--)
	if(tot_count_sph[i] > 0 || tot_count[i] > 0)
	  {
	    printf(" %c  bin=%2d      %d%09d  %d%09d   %6g\n",
		   TimeBinActive[i] ? 'X' : ' ',
		   i,
		   (int) ((tot_count[i] - tot_count_sph[i]) / 1000000000), (int) ((tot_count[i] - tot_count_sph[i]) % 1000000000),
		   (int) (tot_count_sph[i] / 1000000000), (int) (tot_count_sph[i] % 1000000000),
		   i > 0 ? (1 << i) * All.Timebase_interval : 0.0);
	    if(TimeBinActive[i])
	      {
		tot += tot_count[i];
		tot_sph += tot_count_sph[i];
	      }
	  }
      printf("               ------------------------\n");
#ifdef PMGRID
      if(All.PM_Ti_endstep == All.Ti_Current)
	printf("PM-Step. Total: %d%09d  %d%09d    Sum: %d%09d\n",
	       (int) ((tot - tot_sph) / 1000000000), (int) ((tot - tot_sph) % 1000000000),
               (int) (tot_sph / 1000000000), (int) (tot_sph % 1000000000),
               (int) (tot / 1000000000), (int) (tot % 1000000000));
      else
#endif
	printf("Total active:   %d%09d  %d%09d    Sum: %d%09d\n",
               (int) ((tot - tot_sph) / 1000000000), (int) ((tot - tot_sph) % 1000000000),
               (int) (tot_sph / 1000000000), (int) (tot_sph % 1000000000),
               (int) (tot / 1000000000), (int) (tot % 1000000000));
    }

  set_random_numbers();

#ifdef ADAPTIVE_OUTPUT
  if (All.Time > All.AO_TimeOfFirstOutput)
    write_adaptive_output();
#endif

}



void write_cpu_log(void)
{
  double max_CPU_Step[CPU_PARTS], avg_CPU_Step[CPU_PARTS], t0, t1, tsum;

  int i;

  CPU_Step[CPU_MISC] += measure_time();

  All.Cadj_Cpu += CPU_Step[CPU_TREEWALK1] + CPU_Step[CPU_TREEWALK2];

  for(i = 1, CPU_Step[0] = 0; i < CPU_PARTS; i++)
    CPU_Step[0] += CPU_Step[i];

  MPI_Reduce(CPU_Step, max_CPU_Step, CPU_PARTS, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(CPU_Step, avg_CPU_Step, CPU_PARTS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


  if(ThisTask == 0)
    {
      for(i = 0; i < CPU_PARTS; i++)
	avg_CPU_Step[i] /= NTask;

      put_symbol(0.0, 1.0, '#');

      for(i = 1, tsum = 0.0; i < CPU_PARTS; i++)
	{
	  if(max_CPU_Step[i] > 0)
	    {
	      t0 = tsum;
	      t1 = tsum + avg_CPU_Step[i] * (avg_CPU_Step[i] / max_CPU_Step[i]);
	      put_symbol(t0 / avg_CPU_Step[0], t1 / avg_CPU_Step[0], CPU_Symbol[i]);
	      tsum += t1 - t0;

	      t0 = tsum;
	      t1 = tsum + avg_CPU_Step[i] * ((max_CPU_Step[i] - avg_CPU_Step[i]) / max_CPU_Step[i]);
	      put_symbol(t0 / avg_CPU_Step[0], t1 / avg_CPU_Step[0], CPU_SymbolImbalance[i]);
	      tsum += t1 - t0;
	    }
	}

      put_symbol(tsum / max_CPU_Step[0], 1.0, '-');

      fprintf(FdBalance, "Step=%7d  sec=%10.3f  Nf=%10u  %s\n", All.NumCurrentTiStep, max_CPU_Step[0],
	      (unsigned int) GlobNumForceUpdate, CPU_String);
      fflush(FdBalance);
    }

  CPUThisRun += CPU_Step[0];

  for(i = 0; i < CPU_PARTS; i++)
    {
      All.CPU_Sum[i] += avg_CPU_Step[i];
      CPU_Step[i] = 0;
    }

  if(ThisTask == 0)
    {
      fprintf(FdCPU, "Step %d, Time: %g, CPUs: %d\n", All.NumCurrentTiStep, All.Time, NTask);
      fprintf(FdCPU,
	      "total         %10.2f  %5.1f%%\n"
	      "treegrav      %10.2f  %5.1f%%\n"
	      "   treebuild  %10.2f  %5.1f%%\n"
	      "   treeupdate %10.2f  %5.1f%%\n"
	      "   treewalk   %10.2f  %5.1f%%\n"
	      "   treecomm   %10.2f  %5.1f%%\n"
	      "   treeimbal  %10.2f  %5.1f%%\n"
	      "pmgrav        %10.2f  %5.1f%%\n"
	      "sph           %10.2f  %5.1f%%\n"
	      "   density    %10.2f  %5.1f%%\n"
	      "   denscomm   %10.2f  %5.1f%%\n"
	      "   densimbal  %10.2f  %5.1f%%\n"
	      "   hydrofrc   %10.2f  %5.1f%%\n"
	      "   hydcomm    %10.2f  %5.1f%%\n"
	      "   hydimbal   %10.2f  %5.1f%%\n"
	      "   hmaxupdate %10.2f  %5.1f%%\n"
	      "domain        %10.2f  %5.1f%%\n"
	      "potential     %10.2f  %5.1f%%\n"
	      "predict       %10.2f  %5.1f%%\n"
	      "kicks         %10.2f  %5.1f%%\n"
	      "i/o           %10.2f  %5.1f%%\n"
	      "peano         %10.2f  %5.1f%%\n"
	      "sfrcool       %10.2f  %5.1f%%\n"
	      "   cooling    %10.2f  %5.1f%%\n"
	      "   sfr        %10.2f  %5.1f%%\n"
	      "enrichment    %10.2f  %5.1f%%\n"
	      "   evolution  %10.2f  %5.1f%%\n"
	      "molecules     %10.2f  %5.1f%%\n"
	      "blackholes    %10.2f  %5.1f%%\n"
	      "fof           %10.2f  %5.1f%%\n"
	      "lineofsight   %10.2f  %5.1f%%\n"
	      "binning       %10.2f  %5.1f%%\n"
	      "misc          %10.2f  %5.1f%%\n",
	      All.CPU_Sum[CPU_ALL], 100.0,
	      All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2]
	      + All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]
	      + All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2]
	      + All.CPU_Sum[CPU_TREEBUILD] + All.CPU_Sum[CPU_TREEUPDATE]
	      + All.CPU_Sum[CPU_TREEMISC],
	      (All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2]
	       + All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]
	       + All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2]
	       + All.CPU_Sum[CPU_TREEBUILD] + All.CPU_Sum[CPU_TREEUPDATE]
	       + All.CPU_Sum[CPU_TREEMISC]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEBUILD],
	      (All.CPU_Sum[CPU_TREEBUILD]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEUPDATE],
	      (All.CPU_Sum[CPU_TREEUPDATE]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2],
	      (All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV],
	      (All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2],
	      (All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MESH],
	      (All.CPU_Sum[CPU_MESH]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_DENSCOMPUTE] + All.CPU_Sum[CPU_DENSWAIT]
	      + All.CPU_Sum[CPU_DENSCOMM] + All.CPU_Sum[CPU_DENSMISC]
	      + All.CPU_Sum[CPU_HYDCOMPUTE] + All.CPU_Sum[CPU_HYDWAIT] + All.CPU_Sum[CPU_TREEHMAXUPDATE]
	      + All.CPU_Sum[CPU_HYDCOMM] + All.CPU_Sum[CPU_HYDMISC],
	      (All.CPU_Sum[CPU_DENSCOMPUTE] + All.CPU_Sum[CPU_DENSWAIT]
	       + All.CPU_Sum[CPU_DENSCOMM] + All.CPU_Sum[CPU_DENSMISC]
	       + All.CPU_Sum[CPU_HYDCOMPUTE] + All.CPU_Sum[CPU_HYDWAIT] + All.CPU_Sum[CPU_TREEHMAXUPDATE]
	       + All.CPU_Sum[CPU_HYDCOMM] + All.CPU_Sum[CPU_HYDMISC]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_DENSCOMPUTE],
	      (All.CPU_Sum[CPU_DENSCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_DENSCOMM],
	      (All.CPU_Sum[CPU_DENSCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_DENSWAIT],
	      (All.CPU_Sum[CPU_DENSWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_HYDCOMPUTE],
	      (All.CPU_Sum[CPU_HYDCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_HYDCOMM],
	      (All.CPU_Sum[CPU_HYDCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_HYDWAIT],
	      (All.CPU_Sum[CPU_HYDWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEHMAXUPDATE],
	      (All.CPU_Sum[CPU_TREEHMAXUPDATE]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_DOMAIN],
	      (All.CPU_Sum[CPU_DOMAIN]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_POTENTIAL],
	      (All.CPU_Sum[CPU_POTENTIAL]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_DRIFT],
	      (All.CPU_Sum[CPU_DRIFT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TIMELINE],
	      (All.CPU_Sum[CPU_TIMELINE]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SNAPSHOT],
	      (All.CPU_Sum[CPU_SNAPSHOT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_PEANO],
	      (All.CPU_Sum[CPU_PEANO]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_COOLINGSFR], (All.CPU_Sum[CPU_COOLINGSFR]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_COOLING], (All.CPU_Sum[CPU_COOLING]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SFR], (All.CPU_Sum[CPU_SFR]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_ENRICH], (All.CPU_Sum[CPU_ENRICH]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_ENRICHSTEVOL], (All.CPU_Sum[CPU_ENRICHSTEVOL]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MOLECULES], (All.CPU_Sum[CPU_MOLECULES]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_BLACKHOLES], (All.CPU_Sum[CPU_BLACKHOLES]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_FOF], (All.CPU_Sum[CPU_FOF]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_LINEOFSIGHT], (All.CPU_Sum[CPU_LINEOFSIGHT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_BINNING], (All.CPU_Sum[CPU_BINNING]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MISC], (All.CPU_Sum[CPU_MISC]) / All.CPU_Sum[CPU_ALL] * 100);
      fprintf(FdCPU, "\n");
      fflush(FdCPU);
    }
}


void put_symbol(double t0, double t1, char c)
{
  int i, j;

  i = (int) (t0 * CPU_STRING_LEN + 0.5);
  j = (int) (t1 * CPU_STRING_LEN);

  if(i < 0)
    i = 0;
  if(j >= CPU_STRING_LEN)
    j = CPU_STRING_LEN;

  while(i <= j)
    CPU_String[i++] = c;

  CPU_String[CPU_STRING_LEN] = 0;
}


/*! This routine first calls a computation of various global
 * quantities of the particle distribution, and then writes some
 * statistics about the energies in the various particle components to
 * the file FdEnergy.
 */
void energy_statistics(void)
{
  compute_global_quantities_of_system();

  if(ThisTask == 0)
    {
      fprintf(FdEnergy,
	      "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      All.Time, SysState.EnergyInt, SysState.EnergyPot, SysState.EnergyKin, SysState.EnergyIntComp[0],
	      SysState.EnergyPotComp[0], SysState.EnergyKinComp[0], SysState.EnergyIntComp[1],
	      SysState.EnergyPotComp[1], SysState.EnergyKinComp[1], SysState.EnergyIntComp[2],
	      SysState.EnergyPotComp[2], SysState.EnergyKinComp[2], SysState.EnergyIntComp[3],
	      SysState.EnergyPotComp[3], SysState.EnergyKinComp[3], SysState.EnergyIntComp[4],
	      SysState.EnergyPotComp[4], SysState.EnergyKinComp[4], SysState.EnergyIntComp[5],
	      SysState.EnergyPotComp[5], SysState.EnergyKinComp[5], SysState.MassComp[0],
	      SysState.MassComp[1], SysState.MassComp[2], SysState.MassComp[3], SysState.MassComp[4],
	      SysState.MassComp[5]);

      fflush(FdEnergy);
    }
}
