#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>

#include "allvars.h"
#include "proto.h"

#if defined(BG_COOLING) || defined(BG_STELLAR_EVOLUTION)
#include "bg_proto.h"
#include "bg_yields.h"
#endif

/*! \file init.c
 *  \brief code for initialisation of a simulation from initial conditions
 */


/*! This function reads the initial conditions, and allocates storage for the
 *  tree(s). Various variables of the particle data are initialised and An
 *  intial domain decomposition is performed. If SPH particles are present,
 *  the inial SPH smoothing lengths are determined.
 */
void init(void)
{
  int i, j;

  double a3, atime;

#ifdef FOF
  double temp_time;
#endif

#ifdef BLACK_HOLES
  int count_holes = 0;
#endif
#if defined(BG_COOLING) || defined(BG_STELLAR_EVOLUTION)
  int index;

  double mass_fractions[BG_NELEMENTS];
#endif

  All.Time = All.TimeBegin;

  if(RestartFlag == 3 && RestartSnapNum < 0)
    {
      if(ThisTask == 0)
	printf("Need to give the snapshot number if FOF/SUBFIND is selected for output\n");
      endrun(0);
    }

#ifdef ADAPTIVE_OUTPUT
  for(i = 0; i < NumPart; i++)
    {
      P[i].StepsSinceLastOutput=0;
    }
  All.NeedFileRefresh = 1;
#endif

  switch (All.ICFormat)
    {
    case 1:
    case 2:
    case 3:
      read_ic(All.InitCondFile);
      break;
    default:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }

  All.Time = All.TimeBegin;

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      All.Ti_Current = 0;
      a3 = All.Time * All.Time * All.Time;
      atime = All.Time;
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      All.Ti_Current = 0;
      a3 = 1;
      atime = 1;
    }

  set_softenings();

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;
#ifdef FOF
  All.FOFFileCount = 0;
#endif

  if(RestartFlag == 2)
    All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;

#ifdef FOF
    /* Determine FOF group number */
    temp_time = All.TimeOfFirstFOF;
    while(temp_time < All.Time) {
      All.FOFFileCount += 1;
      if (All.ComovingIntegrationOn) {
	temp_time *= All.TimeBetFOF;
      } else {
	temp_time += All.TimeBetFOF;
      }
    }
#endif

  All.TotNumOfForces = 0;
  All.NumForcesSinceLastDomainDecomp = 0;

  All.TopNodeAllocFactor = 0.008;
  All.TreeAllocFactor = 0.7;

  All.Cadj_Cost = 1.0e-30;
  All.Cadj_Cpu = 1.0e-3;

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      check_omega();

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;

  if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
    }

#ifdef BG_STELLAR_EVOLUTION
  if(strcmp(All.SNIa_Model, "GR83") == 0)
    All.SNIa_Mode = 0;
  else if(strcmp(All.SNIa_Model, "Gaussian") == 0)
    All.SNIa_Mode = 1;
  else if(strcmp(All.SNIa_Model, "Efolding") == 0)
    All.SNIa_Mode = 2;
  else if(strcmp(All.SNIa_Model, "GaussianEfolding") == 0)
    All.SNIa_Mode = 3;
  else
    {
      printf("[init]: wrong SNIa model... check your parameter file!\n");
      printf("        valid models are: GR83, Gaussian, Efolding, GaussianEfolding\n");
      endrun(-1);
    }
#endif

#if defined(BG_COOLING) || defined(BG_SFR)
  index = element_index("Hydrogen");
  if(index > -1)
    mass_fractions[index] = All.InitAbundance_Hydrogen;

  index = element_index("Helium");
  if(index > -1)
    mass_fractions[index] = All.InitAbundance_Helium;

  index = element_index("Helium");
  if(index > -1)
    mass_fractions[index] = All.InitAbundance_Helium;

  index = element_index("Carbon");
  if(index > -1)
    mass_fractions[index] = All.InitAbundance_Carbon;

  index = element_index("Nitrogen");
  if(index > -1)
    mass_fractions[index] = All.InitAbundance_Nitrogen;

  index = element_index("Oxygen");
  if(index > -1)
    mass_fractions[index] = All.InitAbundance_Oxygen;

  index = element_index("Neon");
  if(index > -1)
    mass_fractions[index] = All.InitAbundance_Neon;

  index = element_index("Magnesium");
  if(index > -1)
    mass_fractions[index] = All.InitAbundance_Magnesium;

  index = element_index("Silicon");
  if(index > -1)
    mass_fractions[index] = All.InitAbundance_Silicon;

  index = element_index("Iron");
  if(index > -1)
    mass_fractions[index] = All.InitAbundance_Iron;
#endif

  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
	P[i].g.GravAccel[j] = 0;
#ifdef PMGRID
      for(j = 0; j < 3; j++)
	P[i].GravPM[j] = 0;
#endif
      P[i].Ti_begstep = 0;
      P[i].Ti_current = 0;
      P[i].TimeBin = 0;
      P[i].OldAcc = 0;
      P[i].GravCost = 1;
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY)
      if(RestartFlag != 3)
	P[i].p.Potential = 0;
#endif

      if(P[i].Type == 0 && RestartFlag == 0)
	{
#if (defined(BG_SFR) && defined(BG_STELLAR_EVOLUTION)) || defined(BG_COOLING)
	  /* metals */
	  for(j = 0; j < BG_NELEMENTS; j++)
	    SphP[i].Metals[j] = mass_fractions[j] * P[i].Mass;
#endif
#ifdef BG_SFR
#ifdef BG_STELLAR_EVOLUTION
	  //	  SphP[i].Metallicity = 0;
	  SphP[i].Metallicity = 1 - mass_fractions[element_index("Hydrogen")] -
	    mass_fractions[element_index("Helium")];
#ifdef BG_METALSMOOTHING
	  SphP[i].MetallicitySmoothed = 0;
/* #ifdef BG_DUST */
/*           SphP[i].DusticitySmoothed = 0; */
/* #endif */
#endif
#ifdef BG_Z_WEIGHTED_REDSHIFT
	  SphP[i].MetallicityWeightedRedshift = 0;
#endif
#ifdef BG_SNIA_IRON
	  SphP[i].IronFromSNIa = 0;
#endif
/* #ifdef BG_DUST */
/*           SphP[i].Dusticity = 0; */
/* #ifdef BG_DUST_DESTRUCTION_SUBLIMATION */
/* 	  SphP[i].DustGrainSize = 0; */
/* #endif */
/* #endif */
#endif /* BG_STELLAR_EVOLUTION */

	  SphP[i].OnEOS = 0;
#endif /* BG_SFR */

#ifdef BG_EXTRA_ARRAYS
	  SphP[i].MaximumEntropy = 0;
	  SphP[i].MaximumTemperature = 0;
	  SphP[i].TimeMaximumEntropy = 0;
	  SphP[i].TimeMaximumTemperature = 0;
#endif

#if defined(BG_MOL_NETWORK) || defined(BG_MOL_COOLING)
	  if(P[i].Type == 0 && RestartFlag == 0)
	    {
	      SphP[i].x_Hp = All.InitAbundance_Hp;
	      SphP[i].x_Hm = All.InitAbundance_Hm;
	      SphP[i].x_H2 = All.InitAbundance_H2;
	      SphP[i].x_H2p = All.InitAbundance_H2p;
	      SphP[i].x_Hep = All.InitAbundance_Hep;
	      SphP[i].x_Hepp = All.InitAbundance_Hepp;
	      SphP[i].x_Dp = All.InitAbundance_Dp;
	      SphP[i].x_HD = All.InitAbundance_HD;
	      SphP[i].n_e = All.InitNumberDensity_e;
	    }
#endif
	}

#ifdef BG_SFR
#ifdef BG_STELLAR_EVOLUTION
      if(P[i].Type == 4 && RestartFlag == 0)
	StarP[P[i].StarID].SNIaRate = 0;

      if(P[i].Type == 4 && RestartFlag == 2)	/* restart from snapshot */
	{
	  double hubble_a, time_hubble_a, dt, dtime, dtime_in_Gyr;

	  double age_of_star_in_Gyr;

	  MyFloat MetalsReleased[BG_NELEMENTS];

	  MyFloat DustMassReleased, MetalMassReleased, EnergyReleased;

#ifdef BG_SNIA_IRON
	  MyFloat IronFromSNIa;
#endif
	  /* double NumberOfSNIa; */

	  if(All.ComovingIntegrationOn)
	    {
	      hubble_a = All.Hubble * sqrt(All.Omega0 /
					   (All.Time * All.Time * All.Time)
					   + (1 - All.Omega0 - All.OmegaLambda) /
					   (All.Time * All.Time) + All.OmegaLambda);
	      time_hubble_a = All.Time * hubble_a;
	    }
	  else
	    time_hubble_a = 1;

	  /* use a dummy time-step, no time-step defined yet */
	  dt = 1e-5;

	  if(All.ComovingIntegrationOn)
	    dtime = All.Time * dt / time_hubble_a;
	  else
	    dtime = dt;

	  dtime *= All.UnitTime_in_s / All.HubbleParam;	/* convert to  seconds */
	  dtime_in_Gyr = dtime / SEC_PER_MEGAYEAR / 1000;

	  age_of_star_in_Gyr = bg_get_elapsed_time(StarP[P[i].StarID].StarBirthTime, All.Time, 1);	/* note: this is valid for a flat universe! */

	  /* this function encapsulates the stellar evolution. It returns how much
	     mass in each element is returned during stellar evolution of the stellar
	     population between [age_of_star_in_Myr, age_of_star_in_Myr + dtime].
	     It also returns the amount of energy released by SNIa */
#ifdef BG_SNIA_IRON
	  bg_stellar_evolution(age_of_star_in_Gyr - dtime_in_Gyr, dtime_in_Gyr, i, MetalsReleased, &MetalMassReleased,
			       &DustMassReleased, &IronFromSNIa, &EnergyReleased);
	  /* &IronFromSNIa, &EnergyReleased, &NumberOfSNIa); */
#else
	  bg_stellar_evolution(age_of_star_in_Gyr - dtime_in_Gyr, dtime_in_Gyr, i, MetalsReleased, &MetalMassReleased,
			       &DustMassReleased, &EnergyReleased);
	  /* &EnergyReleased, &NumberOfSNIa); */
#endif
	  StarP[P[i].StarID].SNIaRate = NumberOfSNIa / All.HubbleParam / (dtime_in_Gyr * 1.e9);
	}
#endif /* BG_STELLAR_EVOLUTION */

#ifdef BG_SNII_KINETIC_FEEDBACK
      if(P[i].Type == 0 && RestartFlag == 0)
	SphP[i].WindFlag = 0;
#endif
#endif /* BG_SFR */

#ifdef BLACK_HOLES
      if(P[i].Type == 5)
	{
	  count_holes++;

	  if(RestartFlag == 0)
	    P[i].BH_Mass = All.OrigGasMass * All.SeedBHMassOverGasMass;
	}
#endif
    }

#ifdef BG_STELLAR_EVOLUTION
  output_snia_data();
#endif

#ifdef BLACK_HOLES
  MPI_Allreduce(&count_holes, &All.TotBHs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  for(i = 0; i < TIMEBINS; i++)
    TimeBinActive[i] = 1;

  reconstruct_timebins();

#ifdef PMGRID
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      for(j = 0; j < 3; j++)
	{
	  SphP[i].VelPred[j] = P[i].Vel[j];
	  SphP[i].a.HydroAccel[j] = 0;
	}

      SphP[i].e.DtEntropy = 0;

      if(RestartFlag == 0)
	{
#ifndef READ_HSML
	  PPP[i].Hsml = 0;
#endif
	  SphP[i].d.Density = -1;
	  SphP[i].v.DivVel = 0;
	}
#ifdef BG_SFR
      SphP[i].Sfr = 0;
#endif

#ifdef TIMESTEP_LIMITER
      SphP[i].FlagFeedbackDtLimiter = 0;
#endif
/* #if defined(BH_THERMALFEEDBACK) || defined(BG_SNII_THERMAL_FEEDBACK) */
/*       SphP[i].i.Injected_Energy = 0; */
/* #endif */
    }

  All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TotNumPart * All.TreeDomainUpdateFrequency);

  Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */

  TreeReconstructFlag = 1;

  domain_Decomposition();	/* do initial domain decomposition (gives equal numbers of particles) */

  set_softenings();

  if(RestartFlag == 3)
    {
#ifdef FOF
      fof_fof(RestartSnapNum);
#endif
      endrun(0);
    }

  /* will build tree */
  ngb_treebuild();

  All.Ti_Current = 0;

#ifdef PRESSURE_ENTROPY_SPH
  for(i = 0; i < N_gas; i++)
    SphP[i].EntropyVarPred = pow(SphP[i].Entropy, GAMMA_INV);
#endif

  setup_smoothinglengths();

  /* in case of pressure-entropy SPH, we need another iteration here
   * we got an estimate of the density, but the pressure is still wrong
   * if the ICs contain internal energies */

#ifdef PRESSURE_ENTROPY_SPH
  if(header.flag_entropy_instead_u == 0)
    {
      for(i = 0; i < N_gas; i++)
	{
	  SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].cky.WeightedDensity / a3, GAMMA_MINUS1);
	  SphP[i].EntropyVarPred = pow(SphP[i].Entropy, GAMMA_INV);
	}

      density(); /* this is to recompute WeighteDensity with the correct EntropyVarPred */
    }
#endif

  /* at this point, the entropy variable actually contains the 
   * internal energy, read in from the initial conditions file. 
   * Once the density has been computed, we can convert to entropy.
   */

#ifndef PRESSURE_ENTROPY_SPH
  for(i = 0; i < N_gas; i++)
    {
      if(header.flag_entropy_instead_u == 0)
	{
	  if(ThisTask == 0 && i == 0)
	    printf("Converting u -> entropy !\n");
	  SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].d.Density / a3, GAMMA_MINUS1);
	}
    }
#endif

  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      SphP[i].EntropyPred = SphP[i].Entropy;
      SphP[i].e.DtEntropy = 0;
      SphP[i].v.DivVel = 0;
    }
}


/*! This routine computes the mass content of the box and compares it to the
 * specified value of Omega-matter.  If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;

  int i;

  for(i = 0; i < NumPart; i++)
    mass += P[i].Mass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega =
    masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

  if(fabs(omega - All.Omega0) > 1.0e-3)
    {
      if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, All.Omega0);
	  printf("\nI better stop.\n");

	  fflush(stdout);
	}
      endrun(1);
    }
}



/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void setup_smoothinglengths(void)
{
  int i, no, p;

  if(RestartFlag == 0)
    {
      for(i = 0; i < N_gas; i++)
	{
	  no = Father[i];

	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }

#ifndef READ_HSML
#ifndef TWODIMS
	  PPP[i].Hsml =
	    pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else
	  PPP[i].Hsml =
	    pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif
#endif
	  if(PPP[i].Hsml > 10.0 * All.SofteningTable[0])
	    PPP[i].Hsml = All.SofteningTable[0];

	    PPP[i].Hsml = 56.4 * All.SofteningTable[0];
	}
    }

#ifdef BLACK_HOLES
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
	if(P[i].Type == 5)
	  PPP[i].Hsml = All.SofteningTable[5];
    }
#endif

#ifdef BG_SFR
  if(RestartFlag == 0)
    {
      for(i = 0; i < NumPart; i++)
	if(P[i].Type == 4)
	  {
	    no = Father[i];

	    while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	      {
		p = Nodes[no].u.d.father;

		if(p < 0)
		  break;

		no = p;
	      }
	    PPP[i].Hsml =
	      pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
	  }
    }
#endif

  density();

}
