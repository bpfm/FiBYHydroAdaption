#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "kernel.h"

#include "bg_proto.h"

/*! \file timestep.c 
 *  \brief routines for 'kicking' particles in
 *  momentum space and assigning new timesteps
 */

static double fac1, fac2, fac3, hubble_a, atime, a3inv;

static double dt_displacement = 0;

#ifdef PMGRID
static double dt_gravkickA, dt_gravkickB;
#endif


static int HighestOccupiedTimeBin, LowestOccupiedTimeBin;

int timestep_limiter_ngb_treefind(MyDouble searchcenter[3], MyFloat hsml, int target,
				  int *startnode, int mode, int *nexport, int *nsend_local);

void timestep_limiter(void);


#define EXTRASPACEFAC 2

#ifdef TIMESTEP_LIMITER
static struct tstep_limiter_in
{
  MyDouble Pos[3];
  MyFloat Hsml;

  short int TimeBin;

  unsigned int ID;
  int NodeList[NODELISTLENGTH];
} *DtLimiterDataIn, *DtLimiterDataGet;
#endif

/*! This function advances the system in momentum space, i.e. it does apply the 'kick' operation after the
 *  forces have been computed. Additionally, it assigns new timesteps to particles. At start-up, a
 *  half-timestep is carried out, as well as at the end of the simulation. In between, the half-step kick that
 *  ends the previous timestep and the half-step kick for the new timestep are combined into one operation.
 */
void advance_and_find_timesteps(void)
{
  int i, ti_step, ti_step_old, ti_min, tend, tstart, bin, binold, prev, next;

  double aphys;

#ifdef PMGRID
  int j, dt_step;

  double dt_gravkick, dt_hydrokick;
#endif
#ifdef MAKEGLASS
  double disp, dispmax, globmax, dmean, fac, disp2sum, globdisp2sum;
#endif
#if defined(TIME_DEP_ART_VISC) || defined(TIME_DEP_MAGN_DISP)
  double dmin1, dmin2;
#endif


  CPU_Step[CPU_MISC] += measure_time();

  if(All.ComovingIntegrationOn)
    {
      fac1 = 1 / (All.Time * All.Time);
      fac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
      fac3 = pow(All.Time, 3 * (1 - GAMMA) / 2.0);
      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time)
#ifdef DARKENERGY
	+ DarkEnergy_a(All.Time);
#else
	+ All.OmegaLambda;
#endif
      hubble_a = All.Hubble * sqrt(hubble_a);
      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;
    }
  else
    fac1 = fac2 = fac3 = hubble_a = a3inv = atime = 1;

  if(Flag_FullStep || dt_displacement == 0)
    find_dt_displacement_constraint(hubble_a * atime * atime);

#ifdef PMGRID
  if(All.ComovingIntegrationOn)
    dt_gravkickB = get_gravkick_factor(All.PM_Ti_begstep, All.Ti_Current) -
      get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
  else
    dt_gravkickB = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;
#endif


#ifdef MAKEGLASS
  for(i = 0, dispmax = 0, disp2sum = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].g.GravAccel[j] *= -1;
#ifdef PMGRID
	  P[i].GravPM[j] *= -1;
	  P[i].g.GravAccel[j] += P[i].GravPM[j];
	  P[i].GravPM[j] = 0;
#endif
	}

      disp = sqrt(P[i].g.GravAccel[0] * P[i].g.GravAccel[0] +
		  P[i].g.GravAccel[1] * P[i].g.GravAccel[1] + P[i].g.GravAccel[2] * P[i].g.GravAccel[2]);

      disp *= 2.0 / (3 * All.Hubble * All.Hubble);

      disp2sum += disp * disp;

      if(disp > dispmax)
	dispmax = disp;
    }

  MPI_Allreduce(&dispmax, &globmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&disp2sum, &globdisp2sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  dmean = pow(P[0].Mass / (All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);

  if(globmax > dmean)
    fac = dmean / globmax;
  else
    fac = 1.0;

  if(ThisTask == 0)
    {
      printf("\nglass-making:  dmean= %g  global disp-maximum= %g  rms= %g\n\n",
	     dmean, globmax, sqrt(globdisp2sum / All.TotNumPart));
      fflush(stdout);
    }

  for(i = 0, dispmax = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].Vel[j] = 0;
	  P[i].Pos[j] += fac * P[i].g.GravAccel[j] * 2.0 / (3 * All.Hubble * All.Hubble);
	  P[i].g.GravAccel[j] = 0;
	}
    }
#endif




  All.DoDynamicUpdate = ShouldWeDoDynamicUpdate();

  /* Now assign new timesteps and kick */

  if(All.DoDynamicUpdate)
    {
      GlobFlag++;
      DomainNumChanged = 0;
      DomainList = (int *) mymalloc(NTopleaves * sizeof(int));
      if(ThisTask == 0)
	printf("kicks will prepare for dynamic update of tree\n");
    }

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      ti_step = get_timestep(i, &aphys, 0);

      /* make it a power 2 subdivision */
      ti_min = TIMEBASE;
      while(ti_min > ti_step)
	ti_min >>= 1;
      ti_step = ti_min;

      bin = get_timestep_bin(ti_step);
      binold = P[i].TimeBin;

      /*
      if(bin > binold)
	{
	  if(TimeBinActive[bin] == 0)
	    {
	      bin = binold;
	      ti_step = bin ? (1 << bin) : 0;
	    }
	}
      */

      if(bin > binold)   /* timestep wants to increase */
        {
          while(TimeBinActive[bin] == 0 && bin > binold)    /* make sure the new step is synchronized */
            bin--;

          ti_step = bin ? (1 << bin) : 0;
        }

      if(All.Ti_Current >= TIMEBASE)	/* we here finish the last timestep. */
	{
	  ti_step = 0;
	  bin = 0;
	}

      if((TIMEBASE - All.Ti_Current) < ti_step)	/* check that we don't run beyond the end */
	{
	  endrun(888);		/* should not happen */
	  ti_step = TIMEBASE - All.Ti_Current;
	  ti_min = TIMEBASE;
	  while(ti_min > ti_step)
	    ti_min >>= 1;
	  ti_step = ti_min;
	}

      if(bin != binold)
	{
	  TimeBinCount[binold]--;
	  if(P[i].Type == 0)
	    {
	      TimeBinCountSph[binold]--;
#ifdef BG_SFR
	      TimeBinSfr[binold] -= SphP[i].Sfr;
	      TimeBinSfr[bin] += SphP[i].Sfr;
#ifdef BG_POPIII
#ifdef BG_METALSMOOTHING
	      if(SphP[i].MetallicitySmoothed < All.POPIII_MetallicityThreshold)
		{
		  TimeBinPOPIIISfr[binold] -= SphP[i].Sfr;
		  TimeBinPOPIIISfr[bin] += SphP[i].Sfr;
		}
#else
	      if(SphP[i].Metallicity < All.POPIII_MetallicityThreshold)
		{
		  TimeBinPOPIIISfr[binold] -= SphP[i].Sfr;
		  TimeBinPOPIIISfr[bin] += SphP[i].Sfr;
		}
#endif
#endif
#endif
	    }

#ifdef BLACK_HOLES
	  if(P[i].Type == 5)
	    {
	      TimeBin_BH_mass[binold] -= P[i].BH_Mass;
	      TimeBin_BH_dynamicalmass[binold] -= P[i].Mass;
	      TimeBin_BH_Mdot[binold] -= P[i].BH_Mdot;
	      TimeBin_BH_Medd[binold] -= P[i].BH_Mdot / P[i].BH_Mass;
	      TimeBin_BH_mass[bin] += P[i].BH_Mass;
	      TimeBin_BH_dynamicalmass[bin] += P[i].Mass;
	      TimeBin_BH_Mdot[bin] += P[i].BH_Mdot;
	      TimeBin_BH_Medd[bin] += P[i].BH_Mdot / P[i].BH_Mass;
	    }
#endif

	  prev = PrevInTimeBin[i];
	  next = NextInTimeBin[i];

	  if(FirstInTimeBin[binold] == i)
	    FirstInTimeBin[binold] = next;
	  if(LastInTimeBin[binold] == i)
	    LastInTimeBin[binold] = prev;
	  if(prev >= 0)
	    NextInTimeBin[prev] = next;
	  if(next >= 0)
	    PrevInTimeBin[next] = prev;

	  if(TimeBinCount[bin] > 0)
	    {
	      PrevInTimeBin[i] = LastInTimeBin[bin];
	      NextInTimeBin[LastInTimeBin[bin]] = i;
	      NextInTimeBin[i] = -1;
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

	  P[i].TimeBin = bin;
	}

      ti_step_old = binold ? (1 << binold) : 0;

      tstart = P[i].Ti_begstep + ti_step_old / 2;	/* midpoint of old step */
      tend = P[i].Ti_begstep + ti_step_old + ti_step / 2;	/* midpoint of new step */

      P[i].Ti_begstep += ti_step_old;

      do_the_kick(i, tstart, tend, P[i].Ti_begstep);
    }

  if(All.DoDynamicUpdate)
    {
      force_finish_kick_nodes();
      myfree(DomainList);
    }



#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
      All.DoDynamicUpdate = 0;

      ti_step = TIMEBASE;
      while(ti_step > (dt_displacement / All.Timebase_interval))
	ti_step >>= 1;

      if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  /* we only increase if an integer number of steps will bring us to the end */
	  if(((TIMEBASE - All.PM_Ti_endstep) % ti_step) > 0)
	    ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;	/* leave at old step */
	}

      if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	ti_step = 0;

      tstart = (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2;
      tend = All.PM_Ti_endstep + ti_step / 2;

      if(All.ComovingIntegrationOn)
	dt_gravkick = get_gravkick_factor(tstart, tend);
      else
	dt_gravkick = (tend - tstart) * All.Timebase_interval;

      All.PM_Ti_begstep = All.PM_Ti_endstep;
      All.PM_Ti_endstep = All.PM_Ti_begstep + ti_step;

      if(All.ComovingIntegrationOn)
	dt_gravkickB = -get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
      else
	dt_gravkickB =
	  -((All.PM_Ti_begstep + All.PM_Ti_endstep) / 2 - All.PM_Ti_begstep) * All.Timebase_interval;

      for(i = 0; i < NumPart; i++)
	{
	  for(j = 0; j < 3; j++)	/* do the kick */
	    P[i].Vel[j] += P[i].GravPM[j] * dt_gravkick;

	  if(P[i].Type == 0)
	    {
	      dt_step = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0);

	      if(All.ComovingIntegrationOn)
		{
		  dt_gravkickA = get_gravkick_factor(P[i].Ti_begstep, All.Ti_Current) -
		    get_gravkick_factor(P[i].Ti_begstep, P[i].Ti_begstep + dt_step / 2);
		  dt_hydrokick = get_hydrokick_factor(P[i].Ti_begstep, All.Ti_Current) -
		    get_hydrokick_factor(P[i].Ti_begstep, P[i].Ti_begstep + dt_step / 2);
		}
	      else
		dt_gravkickA = dt_hydrokick =
		  (All.Ti_Current - (P[i].Ti_begstep + dt_step / 2)) * All.Timebase_interval;

	      for(j = 0; j < 3; j++)
		SphP[i].VelPred[j] = P[i].Vel[j]
		  + P[i].g.GravAccel[j] * dt_gravkickA
		  + SphP[i].a.HydroAccel[j] * dt_hydrokick + P[i].GravPM[j] * dt_gravkickB;
	    }
	}
    }
#endif

  int lowest_occupied_bin, highest_occupied_bin, n;

  for(n = TIMEBINS - 1; n >=0; n--)
    if(TimeBinCount[n])
      lowest_occupied_bin = n;

  MPI_Allreduce(&lowest_occupied_bin, &LowestOccupiedTimeBin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  for(n = 0; n < TIMEBINS - 1; n++)
    if(TimeBinCount[n])
      highest_occupied_bin = n;

  MPI_Allreduce(&highest_occupied_bin, &HighestOccupiedTimeBin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

#ifdef TIMESTEP_LIMITER
  if(HighestOccupiedTimeBin - LowestOccupiedTimeBin > DT_LIMITER_FACTOR_LOG2)
    {
      if(ThisTask == 0)
        printf("Start time step limiter...\n");

      timestep_limiter();   /* apply modified Saitoh & Makino's time step limiter scheme */

      if(ThisTask == 0)
        printf("Time step limiter done.\n");
    }
  else
    {
      if(ThisTask == 0)
        printf("Skip time-step limiter: HighestOccupiedTimeBin=%d, LowestOccupiedTimeBin=%d\n",
               HighestOccupiedTimeBin, LowestOccupiedTimeBin);
    }
#endif /* End of TIMESTEP_LIMITER conditional */

  CPU_Step[CPU_TIMELINE] += measure_time();
}



void do_the_kick(int i, int tstart, int tend, int tcurrent)
{
  int j;

  MyFloat dv[3];

  double minentropy;

  double dt_entr, dt_gravkick, dt_hydrokick, dt_gravkick2, dt_hydrokick2, dt_entr2;

  if(All.ComovingIntegrationOn)
    {
      dt_entr = (tend - tstart) * All.Timebase_interval;
      dt_entr2 = (tend - tcurrent) * All.Timebase_interval;
      dt_gravkick = get_gravkick_factor(tstart, tend);
      dt_hydrokick = get_hydrokick_factor(tstart, tend);
      dt_gravkick2 = get_gravkick_factor(tcurrent, tend);
      dt_hydrokick2 = get_hydrokick_factor(tcurrent, tend);
    }
  else
    {
      dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;
      dt_gravkick2 = dt_hydrokick2 = dt_entr2 = (tend - tcurrent) * All.Timebase_interval;
    }


  /* do the kick */

  for(j = 0; j < 3; j++)
    {
      dv[j] = P[i].g.GravAccel[j] * dt_gravkick;
      P[i].Vel[j] += dv[j];
    }

  if(P[i].Type == 0)		/* SPH stuff */
    {
      for(j = 0; j < 3; j++)
	{
	  dv[j] += SphP[i].a.HydroAccel[j] * dt_hydrokick;
	  P[i].Vel[j] += SphP[i].a.HydroAccel[j] * dt_hydrokick;

	  SphP[i].VelPred[j] =
	    P[i].Vel[j] - dt_gravkick2 * P[i].g.GravAccel[j] - dt_hydrokick2 * SphP[i].a.HydroAccel[j];
#ifdef PMGRID
	  SphP[i].VelPred[j] += P[i].GravPM[j] * dt_gravkickB;
#endif
	}

      /* In case of cooling, we prevent that the entropy (and
         hence temperature decreases by more than a factor 0.5 */

      if(SphP[i].e.DtEntropy * dt_entr > -0.5 * SphP[i].Entropy)
	SphP[i].Entropy += SphP[i].e.DtEntropy * dt_entr;
      else
	SphP[i].Entropy *= 0.5;

      if(All.MinEgySpec)
	{
	  minentropy = All.MinEgySpec * GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);
	  if(SphP[i].Entropy < minentropy)
	    {
	      SphP[i].Entropy = minentropy;
	      SphP[i].e.DtEntropy = 0;
	    }
	}

      /* In case the timestep increases in the new step, we
         make sure that we do not 'overcool'. */

      dt_entr = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) / 2 * All.Timebase_interval;
      if(SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr < 0.5 * SphP[i].Entropy)
	SphP[i].e.DtEntropy = -0.5 * SphP[i].Entropy / dt_entr;
    }

  if(All.DoDynamicUpdate)
    force_kick_node(i, dv);
}



/*! This function normally (for flag==0) returns the maximum allowed timestep of a particle, expressed in
 *  terms of the integer mapping that is used to represent the total simulated timespan. The physical
 *  acceleration is returned in aphys. The latter is used in conjunction with the PSEUDOSYMMETRIC integration
 *  option, which also makes of the second function of get_timestep. When it is called with a finite timestep
 *  for flag, it returns the physical acceleration that would lead to this timestep, assuming timestep
 *  criterion 0.
 */
int get_timestep(int p,		/*!< particle index */
		 double *aphys,	/*!< acceleration (physical units) */
		 int flag	/*!< either 0 for normal operation, or finite timestep to get corresponding
				   aphys */ )
{
  double ax, ay, az, ac;

  double csnd = 0, dt = 0, dt_courant = 0;

  int ti_step;

#ifdef BLACK_HOLES
  double dt_accr, u_to_temp_fac;

  /* Conversion between internal energy per unit mass and temperature for a fully ionized plasma
     of primordial composition */
  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 *
    All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#endif

#ifdef BG_MOL_NETWORK
  double dt_chem;
#endif

  double dmin1, dmin2;

  if(flag == 0)
    {
      ax = fac1 * P[p].g.GravAccel[0];
      ay = fac1 * P[p].g.GravAccel[1];
      az = fac1 * P[p].g.GravAccel[2];

#ifdef PMGRID
      ax += fac1 * P[p].GravPM[0];
      ay += fac1 * P[p].GravPM[1];
      az += fac1 * P[p].GravPM[2];
#endif

      if(P[p].Type == 0)
	{
	  ax += fac2 * SphP[p].a.HydroAccel[0];
	  ay += fac2 * SphP[p].a.HydroAccel[1];
	  az += fac2 * SphP[p].a.HydroAccel[2];
	}

      ac = sqrt(ax * ax + ay * ay + az * az);	/* this is now the physical acceleration */
      *aphys = ac;
    }
  else
    ac = *aphys;

  if(ac == 0)
    ac = 1.0e-30;

  switch (All.TypeOfTimestepCriterion)
    {
    case 0:
      if(flag > 0)
	{
	  dt = flag * All.Timebase_interval;
	  dt /= hubble_a;	/* convert dloga to physical timestep  */
	  ac = 2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / (dt * dt);
	  *aphys = ac;
	  return flag;
	}
      if(P[p].Type == 0)
	{
#ifdef NOGRAVITY
	  dt = All.CourantFac * sqrt(2 * atime * PPP[p].Hsml * HSML_FACTOR_INV / ac);
#else
	  dt = All.CourantFac * sqrt(2 * atime *
		    DMIN(PPP[p].Hsml / 2.8 * HSML_FACTOR_INV, All.SofteningTable[P[p].Type]) / ac);
#endif
	}
      else
	dt = sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac);

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
      if(P[p].Type == 0)
	dt = sqrt(2 * All.ErrTolIntAccuracy * atime * PPP[p].Hsml / 2.8 / ac);
#else
      if(P[p].Type == 0)
	dt =
	  sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] *
	       pow(P[p].Mass / All.ReferenceGasMass, 1.0 / 3) / ac);
#endif
#endif
      break;

    default:
      endrun(888);
      break;
    }


  if(P[p].Type == 0)
    {
#ifndef PRESSURE_ENTROPY_SPH
      csnd = sqrt(GAMMA * SphP[p].Pressure / SphP[p].d.Density);
#else
      csnd = sqrt(GAMMA * SphP[p].Pressure / SphP[p].cky.WeightedDensity);
#endif

/* #ifdef ALTERNATIVE_VISCOUS_TIMESTEP */
/*       double dmax1, dmax2; */

/*       if(All.ComovingIntegrationOn) */
/* 	dt_courant = All.CourantFac * All.Time * DMAX(PPP[p].Hsml, All.SofteningTable[0]) / (fac3 * csnd); */
/*       else */
/* 	dt_courant = All.CourantFac * DMAX(PPP[p].Hsml, All.SofteningTable[0]) / csnd; */

/*       if(dt_courant > 2 * All.CourantFac * SphP[p].MinViscousDt) */
/* 	dt_courant = 2 * All.CourantFac * SphP[p].MinViscousDt; */
/* #else */
      if(All.ComovingIntegrationOn)
	dt_courant = 2 * All.CourantFac * All.Time * PPP[p].Hsml * HSML_FACTOR_INV / (fac3 * SphP[p].MaxSignalVel);
      else
	dt_courant = 2 * All.CourantFac * PPP[p].Hsml * HSML_FACTOR_INV / SphP[p].MaxSignalVel;
/* #endif */

      if(dt_courant < dt)
	dt = dt_courant;

#ifdef MYFALSE
      dt_viscous = All.CourantFac * SphP[p].MaxViscStep / hubble_a;	/* to convert dloga to physical dt */

      if(dt_viscous < dt)
	dt = dt_viscous;
#endif

    }

#ifdef BLACK_HOLES
  if(P[p].Type == 5)
    {
      if(P[p].BH_Mdot > 0 && P[p].BH_Mass > 0)
	{
	  dt_accr = 0.1 * (All.BlackHoleMinHeatTemp * All.BlackHoleNumberOfNeighboursToHeat / u_to_temp_fac) /
	    (All.BlackHoleFeedbackFactor * All.BlackHoleRadiativeEfficiency * P[p].BH_Mdot * pow(C / All.UnitVelocity_in_cm_per_s, 2));

	  /* dt_accr = 0.25 * P[p].BH_Mass / P[p].BH_Mdot; */
	  if(dt_accr < dt)
	    dt = dt_accr;
	}
    }
#endif

#ifdef BG_MOL_NETWORK
#ifdef BG_MOL_NETWORK_TIMESTEP
  if(P[p].Type == 0)
    {
      dt_chem = SphP[p].Chem_dt / All.UnitTime_in_s;
      dt_chem *= All.HubbleParam;

      if(dt_chem < dt)
	dt = dt_chem;
    }
#endif
#endif


  /* convert the physical timestep to dloga if needed. Note: If comoving integration has not been selected,
     hubble_a=1.
   */
  dt *= hubble_a;

#ifdef ONLY_PM
  dt = All.MaxSizeTimestep;
#endif



  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;


  if(dt >= dt_displacement)
    dt = dt_displacement;

  if(dt < All.MinSizeTimestep)
    {
#ifndef NOSTOP_WHEN_BELOW_MINTIMESTEP
      printf("warning: Timestep wants to be below the limit `MinSizeTimestep'\n");

      if(P[p].Type == 0)
	{
	  printf
	    ("Part-ID=%d  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxcsnd=%g dt0=%g eps=%g\n",
	     (int) P[p].ID, dt, dt_courant * hubble_a, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], PPP[p].Hsml,
	     csnd,
	     sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac) * hubble_a,
	     All.SofteningTable[P[p].Type]);

#ifdef NS_TIMESTEP
	  printf
	    ("Part-ID=%d  dt_NS=%g  A=%g  rho=%g  dotAvisc=%g  dtold=%g, meanpath=%g \n",
	     (int) P[p].ID, dt_NS * hubble_a, SphP[p].Entropy, SphP[p].d.Density,
	     SphP[p].ViscEntropyChange, (P[p].TimeBin ? (1 << P[p].TimeBin) : 0) * All.Timebase_interval,
	     All.IonMeanFreePath *
	     pow((SphP[p].Entropy * pow(SphP[p].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1),
		 2.0) / SphP[p].d.Density);

	  printf("Stressd=(%g|%g|%g) \n", SphP[p].u.s.StressDiag[0], SphP[p].u.s.StressDiag[1],
		 SphP[p].u.s.StressDiag[2]);
	  printf("Stressoffd=(%g|%g|%g) \n", SphP[p].u.s.StressOffDiag[0], SphP[p].u.s.StressOffDiag[1],
		 SphP[p].u.s.StressOffDiag[2]);
#endif


	}
      else
	{
	  printf("Part-ID=%d  dt=%g ac=%g xyz=(%g|%g|%g)\n", (int) P[p].ID, dt, ac, P[p].Pos[0], P[p].Pos[1],
		 P[p].Pos[2]);
	}
      fflush(stdout);
      endrun(888);
#endif
      dt = All.MinSizeTimestep;
    }

  ti_step = (int) (dt / All.Timebase_interval);



  if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
      printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
	     "We better stop.\n"
	     "Task=%d Part-ID=%d dt=%g dtc=%g dtv=%g dtdis=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g) tree=(%g|%g%g)\n\n",
	     ThisTask, (int) P[p].ID, dt, dt_courant, dt, dt_displacement,
	     All.Timebase_interval, ti_step, ac,
	     P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].g.GravAccel[0], P[p].g.GravAccel[1],
	     P[p].g.GravAccel[2]);
#ifdef PMGRID
      printf("pm_force=(%g|%g|%g)\n", P[p].GravPM[0], P[p].GravPM[1], P[p].GravPM[2]);
#endif
      if(P[p].Type == 0)
	printf("hydro-frc=(%g|%g|%g) dens=%g hsml=%g\n", SphP[p].a.HydroAccel[0], SphP[p].a.HydroAccel[1],
	       SphP[p].a.HydroAccel[2], SphP[p].d.Density, PPP[p].Hsml);

      fflush(stdout);
      endrun(818);
    }

  return ti_step;
}


/*! This function computes an upper limit ('dt_displacement') to the global timestep of the system based on
 *  the rms velocities of particles. For cosmological simulations, the criterion used is that the rms
 *  displacement should be at most a fraction MaxRMSDisplacementFac of the mean particle separation. Note that
 *  the latter is estimated using the assigned particle masses, separately for each particle type. If comoving
 *  integration is not used, the function imposes no constraint on the timestep.
 */
void find_dt_displacement_constraint(double hfac /*!<  should be  a^2*H(a)  */ )
{
  int i, j, type, *temp;

  int count[6];

  long long count_sum[6];

  double v[6], v_sum[6], mim[6], min_mass[6];

  double dt, dmean, asmth = 0;

  dt_displacement = All.MaxSizeTimestep;

  if(All.ComovingIntegrationOn)
    {
      for(type = 0; type < 6; type++)
	{
	  count[type] = 0;
	  v[type] = 0;
	  mim[type] = 1.0e30;
	}

      for(i = 0; i < NumPart; i++)
	{
	  v[P[i].Type] += P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2];
	  if(mim[P[i].Type] > P[i].Mass)
	    mim[P[i].Type] = P[i].Mass;
	  count[P[i].Type]++;
	}

      MPI_Allreduce(v, v_sum, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(mim, min_mass, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      temp = (int *) mymalloc(NTask * 6 * sizeof(int));
      MPI_Allgather(count, 6, MPI_INT, temp, 6, MPI_INT, MPI_COMM_WORLD);
      for(i = 0; i < 6; i++)
	{
	  count_sum[i] = 0;
	  for(j = 0; j < NTask; j++)
	    count_sum[i] += temp[j * 6 + i];
	}
      myfree(temp);

#ifdef BG_SFR
      /* add star and gas particles together to treat them on equal footing, using the original gas particle
         spacing. */
      v_sum[0] += v_sum[4];
      count_sum[0] += count_sum[4];
      v_sum[4] = v_sum[0];
      count_sum[4] = count_sum[0];
#ifdef BLACK_HOLES
      v_sum[0] += v_sum[5];
      count_sum[0] += count_sum[5];
      v_sum[5] = v_sum[0];
      count_sum[5] = count_sum[0];
      min_mass[5] = min_mass[0];
#endif
#endif

      for(type = 0; type < 6; type++)
	{
	  if(count_sum[type] > 0)
	    {
	      if(type == 0 || (type == 4 && All.StarformationOn))
		dmean =
		  pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);
	      else
		dmean =
		  pow(min_mass[type] /
		      ((All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);

#ifdef BLACK_HOLES
	      if(type == 5)
		dmean =
		  pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);
#endif

	      dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum[type] / count_sum[type]);

#ifdef PMGRID
	      asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
	      if(((1 << type) & (PLACEHIGHRESREGION)))
		asmth = All.Asmth[1];
#endif
	      if(asmth < dmean)
		dt = All.MaxRMSDisplacementFac * hfac * asmth / sqrt(v_sum[type] / count_sum[type]);
#endif

	      if(ThisTask == 0)
		printf("type=%d  dmean=%g asmth=%g minmass=%g a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",
		       type, dmean, asmth, min_mass[type], All.Time, sqrt(v_sum[type] / count_sum[type]), dt);

	      if(dt < dt_displacement)
		dt_displacement = dt;
	    }
	}

      if(ThisTask == 0)
	printf("displacement time constraint: %g  (%g)\n", dt_displacement, All.MaxSizeTimestep);
    }
}

int get_timestep_bin(int ti_step)
{
  int bin = -1;

  if(ti_step == 0)
    return 0;

  if(ti_step == 1)
    {
      printf("time-step of integer size 1 not allowed\n");
      endrun(112313);
    }

  while(ti_step)
    {
      bin++;
      ti_step >>= 1;
    }

  return bin;
}


#ifdef BG_THERMAL_FEEDBACK_MAKE_PARTICLES_ACTIVE
void synchronize_particle(int i)
{
  int j, no, ti_step_old, ti_step_new, ti_prev_for_bin, ti_next_for_bin;
  int bin = 0, binold, dt_bin;

  double dt_entr, dt_gravkick, dt_hydrokick;
  double time0, time1;

  MyFloat dp[3];
#ifdef SCALAR_FIELD
  MyFloat dp_dm[3];
#endif
  MyDouble dv[3];

  /* the particle is on an active time bin, the new time step and proper kick
     will be computed in timestep.c */
  if(TimeBinActive[P[i].TimeBin])
    return;

  if(!(SphP[i].FlagFeedbackDtLimiter & (1 << BITFLAG_DT_LIMITER)))
    {
      /* store last active time */
      SphP[i].Ti_begstep = P[i].Ti_begstep;

      /* flag particle for SFR calculation */
      SphP[i].FlagFeedbackDtLimiter |= (1 << BITFLAG_DT_LIMITER);

      /*
      printf("[LIMITER, TASK %d]     flagged particle, ID=%d, flag=%d\n",
             ThisTask, P[i].ID, SphP[i].FlagFeedbackDtLimiter);
      */
    }

  /* find the largest active bin with particles in it */
/*   for(j = 0; j < TIMEBINS; j++) */
/*     if(TimeBinActive[j]) */
/*       if(j > bin) */
/* 	bin = j; */

  bin = All.LowestOccupiedTimeBin;

  /* actual particle time bin */
  binold = P[i].TimeBin;

  if(binold > 0)
    {
      dt_bin = (1 << binold);
      ti_prev_for_bin = (All.Ti_Current / dt_bin) * dt_bin;
      ti_next_for_bin = (All.Ti_Current / dt_bin) * dt_bin + dt_bin;
    }
  else
    {
      printf("[SYNC %2d] bin=%d (!!!!!)\n", ThisTask, bin);
      endrun(5172);
    }

  ti_step_old = binold ? (1 << binold) : 0;
  ti_step_new = bin ? (1 << bin) : 0;

  time0 = P[i].Ti_begstep + ti_step_old / 2;
  time1 = All.Ti_Current - ti_step_new / 2;

  if(All.ComovingIntegrationOn)
    {
      dt_entr = (time1 - time0) * All.Timebase_interval;
      dt_gravkick = get_gravkick_factor(time0, time1);
      dt_hydrokick = get_hydrokick_factor(time0, time1);
    }
  else
    dt_entr = dt_gravkick = dt_hydrokick = (time1 - time0) * All.Timebase_interval;

  /* update the entropy at half step */
  dt_entr = time1 - time0;
  dt_entr *= All.Timebase_interval;

  if(SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr < 0)
    {
      printf("[SYNC %2d] PID=%lld, entropy=%g, dtentropy=%g, dt_entr=%g\n",
	     ThisTask, P[i].ID, SphP[i].Entropy, SphP[i].e.DtEntropy, dt_entr);
      endrun(5173);
    }

  SphP[i].Entropy += SphP[i].e.DtEntropy * dt_entr;

  /* update velocity at half step */
  /*
#ifdef PMGRID
  for(j = 0; j < 3; j++)
    dv[j] =
      (P[i].g.GravAccel[j] + P[i].GravPM[j]) * dt_gravkick + SphP[i].a.HydroAccel[j] * dt_hydrokick;
#else
  for(j = 0; j < 3; j++)
    dv[j] = P[i].g.GravAccel[j] * dt_gravkick + SphP[i].a.HydroAccel[j] * dt_hydrokick;
#endif
  */

  /* update node momentum. if not done, node would miss some momentum because
     of the artificial drifting of the velocity of the particle in it */
  /*
  if(All.DoDynamicUpdate)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].Vel[j] += dv[j];

	  dp[j] = P[i].Mass * dv[j];
#ifdef SCALARFIELD
	  if(P[i].Type != 0)
	    dp_dm[j] = P[i].Mass * dv[j];
	  else
	    dp_dm[j] = 0;
#endif 
	}

      no = Father[i];

      while(no >= 0)
	{
	  force_drift_node(no, All.Ti_Current);

	  for(j = 0; j < 3; j++)
	    {
	      Extnodes[no].dp[j] += dp[j];
#ifdef SCALARFIELD
	      Extnodes[no].dp_dm[j] += dp_dm[j];
#endif
	    }

	  no = Nodes[no].u.d.father;
	}
    }
  */

  /* new bin for the particle */
  P[i].TimeBin = bin;

  /* new ti_current and ti_begstep for the particle */
  P[i].Ti_current = All.Ti_Current;
  P[i].Ti_begstep = All.Ti_Current - ti_step_new;

  /* the particle is made active... update the following */
  NumForceUpdate++;
}
#endif


#ifdef TIMESTEP_LIMITER
int timestep_limiter_ngb_treefind(MyDouble searchcenter[3], MyFloat hsml, int target,
				  int *startnode, int mode, int *nexport, int *nsend_local)
{
  int numngb, no, p, task, nexport_save;

  struct NODE *current;

  MyDouble dx, dy, dz, dist;

#ifdef PERIODIC
  MyDouble xtmp;
#endif


  nexport_save = *nexport;

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(P[p].Type > 0)
	    continue;

	  dist = hsml;
	  dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  Ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(123127);

	      if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
		{
		  Exportflag[task] = target;
		  Exportnodecount[task] = NODELISTLENGTH;
		}

	      if(Exportnodecount[task] == NODELISTLENGTH)
		{
		  if(*nexport >= All.BunchSize)
		    {
		      /* out if buffer space. Need to discard work for this particle and interrupt */
		      *nexport = nexport_save;
		      if(nexport_save == 0)
			endrun(14003);	/* in this case, the buffer is too small to process even a single particle */
		      for(task = 0; task < NTask; task++)
			nsend_local[task] = 0;
		      for(no = 0; no < nexport_save; no++)
			nsend_local[DataIndexTable[no].Task]++;
		      return -1;
		    }
		  Exportnodecount[task] = 0;
		  Exportindex[task] = *nexport;
		  DataIndexTable[*nexport].Task = task;
		  DataIndexTable[*nexport].Index = target;
		  DataIndexTable[*nexport].IndexGet = *nexport;

		  *nexport = *nexport + 1;
		  nsend_local[task]++;
		}

	      DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
		DomainNodeIndex[no - (All.MaxPart + MaxNodes)];

	      if(Exportnodecount[task] < NODELISTLENGTH)
		DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;

	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
		  return numngb;
		}
	    }

	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;;
	  dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return numngb;
}


/* Main driver routine that deals with neighbour finding and parallelization 
   (export to other CPUs if needed) */
void timestep_limiter(void)
{
  int i, j, ndone, ndone_flag, dummy, iter;
  int ngrp, sendTask, recvTask, nexport, nimport;
  int place;
  int imax1, imax2;
  MPI_Status status;

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     (1 + EXTRASPACEFAC) * sizeof(struct tstep_limiter_in)));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));

  iter = 0;
  i = FirstActiveParticle;	/* begin with this index */

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */

      DtLimiterDataGet =
	(struct tstep_limiter_in *) mymalloc(EXTRASPACEFAC * All.BunchSize * sizeof(struct tstep_limiter_in));

      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	if(P[i].Type == 0)
	  if(P[i].TimeBin < HighestOccupiedTimeBin - DT_LIMITER_FACTOR_LOG2)
	    if(timestep_limiter_evaluate(i, 0, &nexport, Send_count) < 0)
	      break;

#ifdef MYSORT
      mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#else
      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#endif


      MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);


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

      if(nimport > EXTRASPACEFAC * All.BunchSize)
	{
	  printf("On task=%d, we need to do a realloc unfortunately... nimport=%d  fac*Bunchsize=%d\n",
		 ThisTask, nimport, EXTRASPACEFAC * All.BunchSize);
	  fflush(stdout);
	  DtLimiterDataGet =
	    (struct tstep_limiter_in *) myrealloc(DtLimiterDataGet,
						  IMAX(nimport, EXTRASPACEFAC * All.BunchSize) *
						  sizeof(struct tstep_limiter_in));
	}

      DtLimiterDataIn = (struct tstep_limiter_in *) mymalloc(nexport * sizeof(struct tstep_limiter_in));


      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  DtLimiterDataIn[j].Pos[0] = P[place].Pos[0];
	  DtLimiterDataIn[j].Pos[1] = P[place].Pos[1];
	  DtLimiterDataIn[j].Pos[2] = P[place].Pos[2];

	  DtLimiterDataIn[j].Hsml = PPP[place].Hsml;
	  DtLimiterDataIn[j].TimeBin = P[place].TimeBin;

	  DtLimiterDataIn[j].ID = P[place].ID;

	  memcpy(DtLimiterDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	}


      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&DtLimiterDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct tstep_limiter_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &DtLimiterDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct tstep_limiter_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
		}
	    }
	}


      myfree(DtLimiterDataIn);

      for(j = 0; j < nimport; j++)
	timestep_limiter_evaluate(j, 1, &dummy, &dummy);

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      myfree(DtLimiterDataGet);

      iter++;
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
}


int timestep_limiter_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int startnode;
  unsigned int id;

  int numngb_inbox, listindex = 0;

  MyDouble *pos;

  int j, k, n, bin;
  int bin_new, bin_size_new;
  int bin_old, bin_size_old;
  int ti_begstep_new, ti_begstep_old, ti_endstep_old;

  double dx, dy, dz, h, h2, r2;

  int tstart, tend;

  double dt_entr;
  double dt_gravkick;
  double dt_hydrokick;

  MyDouble dv[3];


  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
      bin = P[target].TimeBin;
      id = P[target].ID;
    }
  else
    {
      pos = DtLimiterDataGet[target].Pos;
      h = DtLimiterDataGet[target].Hsml;
      bin = DtLimiterDataGet[target].TimeBin;
      id = DtLimiterDataGet[target].ID;
    }


  h2 = h * h;

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = DtLimiterDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  /* limited time bin and its size */
  bin_new = bin + DT_LIMITER_FACTOR_LOG2;  /* new limited time bin */
  bin_size_new = (1 << bin_new);        /* new limited time bin size */

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb_inbox = timestep_limiter_ngb_treefind(&pos[0], h, target, &startnode, mode, nexport, nsend_local);

	  if(numngb_inbox < 0)
	    return -1;

	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = Ngblist[n];

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
		  if(!TimeBinActive[P[j].TimeBin])
		    {
		      /* current end time for particle j */
		      ti_endstep_old = P[j].Ti_begstep + (1 << P[j].TimeBin);

		      if(ti_endstep_old > All.Ti_Current + bin_size_new)
			{
			  /*
			  printf("[TASK %2d] INACTIVE: bin_old=%d, bin_new=%d, t_end=%d, t_cur=%d\n",
				 ThisTask, P[j].TimeBin, bin_new, ti_endstep_old, All.Ti_Current);
			  */

                          if(!(SphP[j].FlagFeedbackDtLimiter & (1 << BITFLAG_DT_LIMITER)))
                            {
			      /* store last active time */
                              SphP[j].Ti_begstep = P[j].Ti_begstep;

			      /* flag the particle for SFR calculation */
                              SphP[j].FlagFeedbackDtLimiter |= (1 << BITFLAG_DT_LIMITER);

			      /*
                              printf("[LIMITER, TASK %d]     flagged particle, ID=%d, flag=%d\n",
                                     ThisTask, P[j].ID, SphP[j].FlagFeedbackDtLimiter);
			      */
                            }

			  bin_old = P[j].TimeBin;
			  bin_size_old = (1 << bin_old);

			  ti_begstep_old = P[j].Ti_begstep;

			  /* here look for the next synchronization time *after* the current time going in steps
			     of bin_size_new. it could be that the current time is at synchronization
			     but that would require making the particle active, something that does not seem
			     necessary in Saitoh's scheme */
			  for(k = 0; ti_begstep_old + k * bin_size_new <= All.Ti_Current; k++);

			  ti_begstep_new = ti_begstep_old + (k - 1) * bin_size_new;

			  /* it's needed to synchronize mid-step quantities to the new fictious time-step.
			     the shrinked time-step is bin_size_new, while the old time-step bin_size_old */
			  tstart = ti_begstep_old + bin_size_old / 2;
			  tend = ti_begstep_new + bin_size_new / 2;

			  if(All.ComovingIntegrationOn)
			    {
			      dt_entr = (tend - tstart) * All.Timebase_interval;
			      dt_gravkick = get_gravkick_factor(tstart, tend);
			      dt_hydrokick = get_hydrokick_factor(tstart, tend);
			    }
			  else
			    {
			      dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;
			    }

			  for(k = 0; k < 3; k++)
			    {
			      dv[k] = P[j].g.GravAccel[k] * dt_gravkick;
			      P[j].Vel[k] += dv[k];
			    }

			  if(P[j].Type == 0)        /* SPH stuff */
			    {
			      /* note that VelPred is already at the right time */
			      for(k = 0; k < 3; k++)
				P[j].Vel[k] += SphP[j].a.HydroAccel[k] * dt_hydrokick;

			      SphP[j].Entropy += SphP[j].e.DtEntropy * dt_entr;
			    }

			  P[j].TimeBin = bin_new;     /* new time bin */
			  P[j].Ti_begstep = ti_begstep_new;    /* new beginning of step */

			  //NTouched++;

			  TimeBinCount[bin_old]--;
			  if(P[j].Type == 0)
			    {
			      TimeBinCountSph[bin_old]--;
#ifdef BG_SFR
			      TimeBinSfr[bin_old] -= SphP[j].Sfr;
			      TimeBinSfr[bin_new] += SphP[j].Sfr;
#endif
			    }

			  int prev, next;

			  prev = PrevInTimeBin[j];
			  next = NextInTimeBin[j];

			  if(FirstInTimeBin[bin_old] == j)
			    FirstInTimeBin[bin_old] = next;
			  if(LastInTimeBin[bin_old] == j)
			    LastInTimeBin[bin_old] = prev;
			  if(prev >= 0)
			    NextInTimeBin[prev] = next;
			  if(next >= 0)
			    PrevInTimeBin[next] = prev;

			  if(TimeBinCount[bin_new] > 0)
			    {
			      PrevInTimeBin[j] = LastInTimeBin[bin_new];
			      NextInTimeBin[LastInTimeBin[bin_new]] = j;
			      NextInTimeBin[j] = -1;
			      LastInTimeBin[bin_new] = j;
			    }
			  else
			    {
			      FirstInTimeBin[bin_new] = LastInTimeBin[bin_new] = j;
			      PrevInTimeBin[j] = NextInTimeBin[j] = -1;
			    }

			  TimeBinCount[bin_new]++;
			  if(P[j].Type == 0)
			    TimeBinCountSph[bin_new]++;
			}
		    }
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = DtLimiterDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  return 0;
}
#endif /* TIMESTEP_LIMITER */
