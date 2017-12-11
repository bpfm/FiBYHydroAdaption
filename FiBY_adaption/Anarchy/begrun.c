#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"

#ifdef BG_COOLING
#include "bg_cooling.h"
#endif

#if defined(BG_COOLING) || defined (BG_STELLAR_EVOLUTION)
#include "bg_proto.h"
#include "bg_vars.h"
#endif

#if defined(BG_MOL_NETWORK) && defined(BG_MOL_NETWORK_TABULATED_RATES)
#include "bg_molecules.h"
#endif

/*! \file begrun.c
 *  \brief initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameterfile is read in and parsed, the initial
 *  conditions or restart files are read, and global variables are initialized
 *  to their proper values.
 */



/*! This function performs the initial set-up of the simulation. First, the
 *  parameterfile is set, then routines for setting units, reading
 *  ICs/restart-files are called, auxialiary memory is allocated, etc.
 */
void begrun(void)
{
  struct global_data_all_processes all;

#ifdef BG_SNII_THERMAL_FEEDBACK
  double snii_energy;
/* #if defined(BG_POPIII) && defined(BG_POPIII_THERMAL_FEEDBACK) */
/*   double pisn_energy; */
/* #endif */
#endif

  if(ThisTask == 0)
    {
      /*    printf("\nThis is P-Gadget, version `%s', svn-revision `%s'.\n", GADGETVERSION, svn_version()); */
      printf("\nThis is P-Gadget, version `%s'.\n", GADGETVERSION);
      printf("\nRunning on %d processors.\n", NTask);
    }

#if defined(X86FIX) && defined(SOFTDOUBLEDOUBLE)
  x86_fix();			/* disable 80bit treatment of internal FPU registers in favour of proper IEEE 64bit double precision arithmetic */
#endif

#ifdef PEDANTIC_MEMORY_HANDLER
  mymalloc_init();
#ifdef USE_HDF5_FIX
  hdf5_mymalloc_init();
#endif
#endif

  read_parameter_file(ParameterFile);	/* ... read in parameters for this run */


/* #ifdef FORCE_DOMAIN_DECOMPOSITION */
/*   PartsPerSecond = 1E20; */
/* #endif */


#ifdef DEBUG
  write_pid_file();
  enable_core_dumps_and_fpu_exceptions();
#endif

  set_units();

#if defined(BG_COOLING) || defined(BG_STELLAR_EVOLUTION)
  InitChemistry();

  if(ThisTask == 0)
    {
      printf("Chemistry initialized.\n");
      fflush(stdout);
    }
#endif

#if defined(BG_MOL_NETWORK) && defined(BG_MOL_NETWORK_TABULATED_RATES)
  init_rate_coefficient_tables();

  if(ThisTask == 0)
    {
      printf("Molecular network rate coefficient tables initialized.\n");
      fflush(stdout);
    }
#endif

#if defined(REPOSITION_ON_POTMIN) && !defined(COMPUTE_POTENTIAL_ENERGY)
  printf("Specified REPOSITION_ON_POTMIN without COMPUTE_POTENTIAL_ENERGY.  Not possible!\n");
  endrun(1);
#endif

#ifdef BG_STELLAR_EVOLUTION
  init_imf();

/* #ifdef BG_POPIII */
/*   init_popiii_imf(); */
/* #endif */

  if(ThisTask == 0)
    {
      printf("IMF initialized.\n");
      fflush(stdout);
    }

  init_yields();

  if(ThisTask == 0)
    {
      printf("Yields initialized.\n");
      fflush(stdout);
    }
#endif

#ifdef BG_SNII_THERMAL_FEEDBACK
  /* SNII energy */
  snii_energy = integrate_imf(log10(All.SNII_MinMass_MSUN), log10(All.SNII_MaxMass_MSUN), 0.0, 0) *
    All.SNII_Energy_ERG * All.UnitMass_in_g / SOLAR_MASS / All.UnitEnergy_in_cgs;

  All.SNII_AvailableEnergyPerUnitMass = snii_energy;
  All.SNII_EnergyPerUnitMass = All.SNII_EnergyPerUnitMass_ERGperG * All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  if(ThisTask == 0)
    {
      printf("SNII energy initialized.\n");
      printf("Number of SNII = %f /solar_mass.\n",
	     integrate_imf(log10(All.SNII_MinMass_MSUN), log10(All.SNII_MaxMass_MSUN), 0.0, 0));
      printf("Energy from SNII = %e erg/solar_mass, [%e code units].\n",
	     All.SNII_AvailableEnergyPerUnitMass * All.UnitEnergy_in_cgs * SOLAR_MASS / All.UnitMass_in_g,
	     All.SNII_AvailableEnergyPerUnitMass);
      printf("Energy given = %e erg/solar_mass, [%e code units].\n",
	     All.SNII_EnergyPerUnitMass_ERGperG * SOLAR_MASS, All.SNII_EnergyPerUnitMass);
      fflush(stdout);
    }

/* #if defined(BG_POPIII) && defined(BG_POPIII_THERMAL_FEEDBACK) */
/*   /\* SNII/PISN energy *\/ */
/*   snii_energy = integrate_popiii_imf(log10(All.POPIII_SNII_MinMass_MSUN), log10(All.POPIII_SNII_MaxMass_MSUN), 0.0, 0) * */
/*     All.POPIII_SNII_Energy_ERG * All.UnitMass_in_g / SOLAR_MASS / All.UnitEnergy_in_cgs; */

/*   pisn_energy = integrate_popiii_imf(log10(All.POPIII_PISN_MinMass_MSUN), log10(All.POPIII_PISN_MaxMass_MSUN), 0.0, 0) * */
/*     All.POPIII_PISN_Energy_ERG * All.UnitMass_in_g / SOLAR_MASS / All.UnitEnergy_in_cgs; */

/*   All.POPIII_AvailableEnergyPerUnitMass = snii_energy + pisn_energy; */
/*   All.POPIII_EnergyPerUnitMass = All.POPIII_EnergyPerUnitMass_ERGperG * All.UnitMass_in_g / All.UnitEnergy_in_cgs; */

/*   if(ThisTask == 0) */
/*     { */
/*       printf("\nPOPIII energy initialized.\n"); */
/*       printf("Number of SNII = %f /solar_mass.\n", */
/*              integrate_popiii_imf(log10(All.POPIII_SNII_MinMass_MSUN), log10(All.POPIII_SNII_MaxMass_MSUN), 0.0, 0)); */
/*       printf("Number of PISN = %f /solar_mass.\n", */
/*              integrate_popiii_imf(log10(All.POPIII_PISN_MinMass_MSUN), log10(All.POPIII_PISN_MaxMass_MSUN), 0.0, 0)); */
/*       printf("Energy from POPIII = %e erg/solar_mass, [%e code units], SNII/PISN = %e.\n", */
/*              All.POPIII_AvailableEnergyPerUnitMass * All.UnitEnergy_in_cgs * SOLAR_MASS / All.UnitMass_in_g, */
/*              All.POPIII_AvailableEnergyPerUnitMass, snii_energy / pisn_energy); */
/*       printf("Energy given = %e erg/solar_mass, [%e code units].\n", */
/* 	     All.POPIII_EnergyPerUnitMass_ERGperG * SOLAR_MASS, All.POPIII_EnergyPerUnitMass); */
/*       fflush(stdout); */
/*     } */
/* #endif */
#endif


#ifdef BG_COOLING
  All.Time = All.TimeBegin;
  InitCool();

  int z_index;

  float d_z;

  if(ThisTask == 0)
    printf("\nLoading cooling tables...");

  get_redshift_index(1 / All.Time - 1, &z_index, &d_z);
  LoadCoolingTables(z_index);

  if(ThisTask == 0)
    printf(" done\n");

#ifdef BG_COOLING_SHIELDING
  All.UVBG_PhysDensThresh = All.UVBG_PhysDensThresh_HpCM3 *
    PROTONMASS / HYDROGEN_MASSFRAC / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);

  LoadCollisionalCoolingTables();
#endif
#endif


#ifdef USE_HDF5_FIX
  H5close();
  hdf5_memory_cleanup();
#endif

#ifdef PERIODIC
  ewald_init();
#endif

#ifdef PERIODIC
  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
  fwa_init();
#endif
#endif

#ifdef TIME_DEP_ART_VISC
  All.ViscSource = All.ViscSource0 / log((GAMMA + 1) / (GAMMA - 1));
  All.DecayTime = 1 / All.DecayLength * sqrt((GAMMA - 1) / 2 * GAMMA);
#endif

  open_log_files();

#ifdef BG_SFR
  write_log_files_header();
#endif

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42);	/* start-up seed */

#ifdef PMGRID
  if(RestartFlag != 3)
    long_range_init();
#endif

  All.TimeLastRestartFile = CPUThisRun;

  if(RestartFlag == 0 || RestartFlag >= 2)
    {
      set_random_numbers();

      init();			/* ... read in initial model */
    }
  else
    {
      all = All;		/* save global variables. (will be read from restart file) */

      restart(RestartFlag);	/* ... read restart file. Note: This also resets 
				   all variables in the struct `All'. 
				   However, during the run, some variables in the parameter
				   file are allowed to be changed, if desired. These need to 
				   copied in the way below.
				   Note:  All.PartAllocFactor is treated in restart() separately.  
				 */

      All.MinSizeTimestep = all.MinSizeTimestep;
      All.MaxSizeTimestep = all.MaxSizeTimestep;
      All.BufferSize = all.BufferSize;
      All.TimeLimitCPU = all.TimeLimitCPU;
      All.ResubmitOn = all.ResubmitOn;
      All.TimeBetSnapshot = all.TimeBetSnapshot;
      All.TimeBetStatistics = all.TimeBetStatistics;
      All.CpuTimeBetRestartFile = all.CpuTimeBetRestartFile;
      All.ErrTolIntAccuracy = all.ErrTolIntAccuracy;
      All.MaxRMSDisplacementFac = all.MaxRMSDisplacementFac;

      All.ErrTolForceAcc = all.ErrTolForceAcc;
      All.TypeOfTimestepCriterion = all.TypeOfTimestepCriterion;
      All.TypeOfOpeningCriterion = all.TypeOfOpeningCriterion;
      All.NumFilesWrittenInParallel = all.NumFilesWrittenInParallel;
      All.TreeDomainUpdateFrequency = all.TreeDomainUpdateFrequency;

      All.OutputListOn = all.OutputListOn;
      All.CourantFac = all.CourantFac;

      All.OutputListLength = all.OutputListLength;
      memcpy(All.OutputListTimes, all.OutputListTimes, sizeof(double) * All.OutputListLength);

#ifdef BG_COOLING
      All.MetDepCoolingOn = all.MetDepCoolingOn;
      All.REION_H_ZCenter = all.REION_H_ZCenter;
#endif

#ifdef BG_SFR
      if(ThisTask == 0)
	printf("EoS Gamma Effective: was %e, now is %e\n",
	       All.SF_EOSGammaEffective, all.SF_EOSGammaEffective);

      All.SF_EOSGammaEffective = all.SF_EOSGammaEffective;

      if(ThisTask == 0)
	printf("EoS Threshold Energy: was %e, now is %e\n",
	       All.SF_EOSEnergyAtThreshold_ERG, all.SF_EOSEnergyAtThreshold_ERG);

      All.SF_EOSEnergyAtThreshold_ERG = all.SF_EOSEnergyAtThreshold_ERG;

      if(ThisTask == 0)
	printf("Schmidt Law Normalization: was %e, now is %e\n",
	       All.SF_SchmidtLawCoeff_MSUNpYRpKPC2, all.SF_SchmidtLawCoeff_MSUNpYRpKPC2);

      All.SF_SchmidtLawCoeff_MSUNpYRpKPC2 = all.SF_SchmidtLawCoeff_MSUNpYRpKPC2;

      set_units_sfr();
#endif

#ifdef BG_STELLAR_EVOLUTION
      strcpy(All.IMF_Model, all.IMF_Model);
      strcpy(All.SNIa_Model, all.SNIa_Model);
      All.SNIa_MassTransferOn = all.SNIa_MassTransferOn;
      All.SNIa_EnergyTransferOn = all.SNIa_EnergyTransferOn;
      All.AGB_MassTransferOn = all.AGB_MassTransferOn;
#endif

#ifdef BG_SNII_THERMAL_FEEDBACK
      All.SNII_EnergyTransferOn = all.SNII_EnergyTransferOn;
      All.SNII_Energy_ERG = all.SNII_Energy_ERG;
      All.SNII_EnergyFraction = all.SNII_EnergyFraction;
      All.SNII_EnergyPerUnitMass_ERGperG = all.SNII_EnergyPerUnitMass_ERGperG;
#endif

/* #if defined(BG_POPIII) && defined(BG_POPIII_THERMAL_FEEDBACK) */
/*       All.POPIII_EnergyTransferOn = all.POPIII_EnergyTransferOn; */
/*       All.POPIII_EnergyFraction = all.POPIII_EnergyFraction; */
/*       All.POPIII_EnergyPerUnitMass_ERGperG = all.POPIII_EnergyPerUnitMass_ERGperG; */

/*       All.POPIII_SNII_Energy_ERG = all.POPIII_SNII_Energy_ERG; */
/*       All.POPIII_PISN_Energy_ERG = all.POPIII_PISN_Energy_ERG; */
/* #endif */

#ifdef BG_SNII_KINETIC_FEEDBACK
      All.SNII_WindOn = all.SNII_WindOn;
      All.SNII_WindIsotropicOn = all.SNII_WindIsotropicOn;
      All.SNII_WindMassLoading = all.SNII_WindMassLoading;
      All.SNII_WindSpeed_KMpS = all.SNII_WindSpeed_KMpS;
#ifdef BG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING
      All.SNII_WindDecouplingTime_YR = all.SNII_WindDecouplingTime_YR;
#endif
#endif

#if defined(BG_SNII_KINETIC_FEEDBACK) || defined(BG_SNII_THERMAL_FEEDBACK)
      All.SNII_WindDelay_YR = all.SNII_WindDelay_YR;
#endif

#ifdef TIME_DEP_ART_VISC
      All.ViscSource = all.ViscSource;
      All.ViscSource0 = all.ViscSource0;
      All.DecayTime = all.DecayTime;
      All.DecayLength = all.DecayLength;
      All.AlphaMin = all.AlphaMin;
#endif


#ifdef DARKENERGY
      All.DarkEnergyParam = all.DarkEnergyParam;
#endif

      strcpy(All.ResubmitCommand, all.ResubmitCommand);
      strcpy(All.OutputListFilename, all.OutputListFilename);
      strcpy(All.OutputDir, all.OutputDir);
      strcpy(All.RestartFile, all.RestartFile);
      strcpy(All.EnergyFile, all.EnergyFile);
      strcpy(All.InfoFile, all.InfoFile);
      strcpy(All.CpuFile, all.CpuFile);
      strcpy(All.TimingsFile, all.TimingsFile);
      strcpy(All.SnapshotFileBase, all.SnapshotFileBase);

      if(All.TimeMax != all.TimeMax)
	readjust_timebase(All.TimeMax, all.TimeMax);

#ifdef NO_TREEDATA_IN_RESTART
      /* if this is not activated, the tree was stored in the restart-files,
         which also allocated the storage for it already */

      /* ensures that domain reconstruction will be done and new tree will be constructed */
      All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TotNumPart * All.TreeDomainUpdateFrequency);
#endif
    }

#ifdef PMGRID
  long_range_init_regionsize();
#endif

  reconstruct_timebins();


  if(All.ComovingIntegrationOn)
    init_drift_table();

  if(RestartFlag == 2)
    {
      All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 100);
#ifdef BG_OUTPUT_GRID
      All.Ti_nextgridoutput = find_next_gridoutputtime(All.Ti_Current + 100);
#endif
#ifdef OUTPUTLINEOFSIGHT
      All.Ti_nextlineofsight = find_next_losoutputtime(All.Ti_Current + 100);
#endif
#ifdef FOF
      All.Ti_nextfof = find_next_fofoutputtime(All.Ti_Current + 100);
#endif
    }
  else
    {
      All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);
#ifdef BG_OUTPUT_GRID
      All.Ti_nextgridoutput = find_next_gridoutputtime(All.Ti_Current);
#endif
#ifdef OUTPUTLINEOFSIGHT
      All.Ti_nextlineofsight = find_next_losoutputtime(All.Ti_Current);
#endif
#ifdef FOF
      All.Ti_nextfof = find_next_fofoutputtime(All.Ti_Current);
#endif
    }

  All.TimeLastRestartFile = CPUThisRun;
}




/*! Computes conversion factors between internal code units and the
 *  cgs-system.
 */
void set_units(void)
{
#ifdef STATICNFW
  double Mtot;
#endif

  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;

  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  if(ThisTask == 0)
    {
      printf("\nHubble (internal units) = %g\n", All.Hubble);
      printf("G (internal units) = %g\n", All.G);
      printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
      printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
      printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
      printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
      printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
      printf("\n");
    }

#if defined(BG_SFR) || defined(BG_COOLING)
  All.MinEgySpec = All.MinGasU_ERG * All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#endif

#ifdef BG_SFR
  set_units_sfr();
#endif

#ifdef STATICNFW
  R200 = pow(NFW_M200 * All.G / (100 * All.Hubble * All.Hubble), 1.0 / 3);
  Rs = R200 / NFW_C;
  Dc = 200.0 / 3 * NFW_C * NFW_C * NFW_C / (log(1 + NFW_C) - NFW_C / (1 + NFW_C));
  RhoCrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  V200 = 10 * All.Hubble * R200;
  if(ThisTask == 0)
    printf("V200= %g\n", V200);

  fac = 1.0;
  Mtot = enclosed_mass(R200);
  if(ThisTask == 0)
    printf("M200= %g\n", Mtot);
  fac = V200 * V200 * V200 / (10 * All.G * All.Hubble) / Mtot;
  Mtot = enclosed_mass(R200);
  if(ThisTask == 0)
    printf("M200= %g\n", Mtot);
#endif
}

#ifdef STATICNFW
/*! auxiliary function for static NFW potential
 */
double enclosed_mass(double R)
{
  /* Eps is in units of Rs !!!! */

  if(R > Rs * NFW_C)
    R = Rs * NFW_C;

  return fac * 4 * M_PI * RhoCrit * Dc *
    (-(Rs * Rs * Rs * (1 - NFW_Eps + log(Rs) - 2 * NFW_Eps * log(Rs) + NFW_Eps * NFW_Eps * log(NFW_Eps * Rs)))
     / ((NFW_Eps - 1) * (NFW_Eps - 1)) +
     (Rs * Rs * Rs *
      (Rs - NFW_Eps * Rs - (2 * NFW_Eps - 1) * (R + Rs) * log(R + Rs) +
       NFW_Eps * NFW_Eps * (R + Rs) * log(R + NFW_Eps * Rs))) / ((NFW_Eps - 1) * (NFW_Eps - 1) * (R + Rs)));
}
#endif



/*!  This function opens various log-files that report on the status and
 *   performance of the simulstion. On restart from restart-files
 *   (start-option 1), the code will append to these files.
 */
void open_log_files(void)
{
  char mode[2], buf[200];

  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");

#if defined(BLACK_HOLES) && defined(VERBOSE_LOGFILES)
  /* Note: This is done by everyone */
  sprintf(buf, "%sblackhole_details_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBlackHolesDetails = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%sblackhole_heating_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBlackHolesFeedback = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }


#endif

#ifdef BG_SNII_THERMAL_FEEDBACK
  /* Note: This is done by everyone */
#ifdef BG_VERBOSE_LOGFILES
  sprintf(buf, "%sthermal_feedback_%d.txt", All.OutputDir, ThisTask);
  if(!(FdThermalFeedback = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif
#endif

  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

#ifdef BLACK_HOLES
#ifdef FOF
  sprintf(buf, "%s/seeds.txt", All.OutputDir);
  if(!(FdBlackHolesSeeds = fopen(buf, mode)))
   {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
   }
#endif
#ifdef BH_DEBUG
  sprintf(buf, "%s/blackhole_energy.txt", All.OutputDir);
  if(!(FdBlackHolesEnergy = fopen(buf, mode)))
   {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
   }
#endif

#endif


  sprintf(buf, "%s%s", All.OutputDir, All.CpuFile);
  if(!(FdCPU = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.InfoFile);
  if(!(FdInfo = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.EnergyFile);
  if(!(FdEnergy = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.TimingsFile);
  if(!(FdTimings = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "balance.txt");
  if(!(FdBalance = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  fprintf(FdBalance, "\n");
  fprintf(FdBalance, "Treewalk1      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWALK1],
	  CPU_SymbolImbalance[CPU_TREEWALK1]);
  fprintf(FdBalance, "Treewalk2      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWALK2],
	  CPU_SymbolImbalance[CPU_TREEWALK2]);
  fprintf(FdBalance, "Treewait1      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWAIT1],
	  CPU_SymbolImbalance[CPU_TREEWAIT1]);
  fprintf(FdBalance, "Treewait2      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWAIT2],
	  CPU_SymbolImbalance[CPU_TREEWAIT2]);
  fprintf(FdBalance, "Treesend       = '%c' / '%c'\n", CPU_Symbol[CPU_TREESEND],
	  CPU_SymbolImbalance[CPU_TREESEND]);
  fprintf(FdBalance, "Treerecv       = '%c' / '%c'\n", CPU_Symbol[CPU_TREERECV],
	  CPU_SymbolImbalance[CPU_TREERECV]);
  fprintf(FdBalance, "Treebuild      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEBUILD],
	  CPU_SymbolImbalance[CPU_TREEBUILD]);
  fprintf(FdBalance, "Treeupdate     = '%c' / '%c'\n", CPU_Symbol[CPU_TREEUPDATE],
	  CPU_SymbolImbalance[CPU_TREEUPDATE]);
  fprintf(FdBalance, "Treehmaxupdate = '%c' / '%c'\n", CPU_Symbol[CPU_TREEHMAXUPDATE],
	  CPU_SymbolImbalance[CPU_TREEHMAXUPDATE]);
  fprintf(FdBalance, "Treemisc =       '%c' / '%c'\n", CPU_Symbol[CPU_TREEMISC],
	  CPU_SymbolImbalance[CPU_TREEMISC]);
  fprintf(FdBalance, "Domain decomp  = '%c' / '%c'\n", CPU_Symbol[CPU_DOMAIN],
	  CPU_SymbolImbalance[CPU_DOMAIN]);
  fprintf(FdBalance, "Density compute= '%c' / '%c'\n", CPU_Symbol[CPU_DENSCOMPUTE],
	  CPU_SymbolImbalance[CPU_DENSCOMPUTE]);
  fprintf(FdBalance, "Density imbal  = '%c' / '%c'\n", CPU_Symbol[CPU_DENSWAIT],
	  CPU_SymbolImbalance[CPU_DENSWAIT]);
  fprintf(FdBalance, "Density commu  = '%c' / '%c'\n", CPU_Symbol[CPU_DENSCOMM],
	  CPU_SymbolImbalance[CPU_DENSCOMM]);
  fprintf(FdBalance, "Density misc   = '%c' / '%c'\n", CPU_Symbol[CPU_DENSMISC],
	  CPU_SymbolImbalance[CPU_DENSMISC]);
  fprintf(FdBalance, "Hydro compute  = '%c' / '%c'\n", CPU_Symbol[CPU_HYDCOMPUTE],
	  CPU_SymbolImbalance[CPU_HYDCOMPUTE]);
  fprintf(FdBalance, "Hydro imbalance= '%c' / '%c'\n", CPU_Symbol[CPU_HYDWAIT],
	  CPU_SymbolImbalance[CPU_HYDWAIT]);
  fprintf(FdBalance, "Hydro comm     = '%c' / '%c'\n", CPU_Symbol[CPU_HYDCOMM],
	  CPU_SymbolImbalance[CPU_HYDCOMM]);
  fprintf(FdBalance, "Hydro misc     = '%c' / '%c'\n", CPU_Symbol[CPU_HYDMISC],
	  CPU_SymbolImbalance[CPU_HYDMISC]);
  fprintf(FdBalance, "Drifts         = '%c' / '%c'\n", CPU_Symbol[CPU_DRIFT], CPU_SymbolImbalance[CPU_DRIFT]);
  fprintf(FdBalance, "Blackhole      = '%c' / '%c'\n", CPU_Symbol[CPU_BLACKHOLES],
	  CPU_SymbolImbalance[CPU_BLACKHOLES]);
  fprintf(FdBalance, "Kicks          = '%c' / '%c'\n", CPU_Symbol[CPU_TIMELINE],
	  CPU_SymbolImbalance[CPU_TIMELINE]);
  fprintf(FdBalance, "Potential      = '%c' / '%c'\n", CPU_Symbol[CPU_POTENTIAL],
	  CPU_SymbolImbalance[CPU_POTENTIAL]);
  fprintf(FdBalance, "PM             = '%c' / '%c'\n", CPU_Symbol[CPU_MESH], CPU_SymbolImbalance[CPU_MESH]);
  fprintf(FdBalance, "Peano-Hilbert  = '%c' / '%c'\n", CPU_Symbol[CPU_PEANO], CPU_SymbolImbalance[CPU_PEANO]);
  fprintf(FdBalance, "Cooling & SFR  = '%c' / '%c'\n", CPU_Symbol[CPU_COOLINGSFR],
	  CPU_SymbolImbalance[CPU_COOLINGSFR]);
  fprintf(FdBalance, "Cooling        = '%c' / '%c'\n", CPU_Symbol[CPU_COOLING],
	  CPU_SymbolImbalance[CPU_COOLING]);
  fprintf(FdBalance, "Star formation = '%c' / '%c'\n", CPU_Symbol[CPU_SFR],
	  CPU_SymbolImbalance[CPU_SFR]);
  fprintf(FdBalance, "Enrichment     = '%c' / '%c'\n", CPU_Symbol[CPU_ENRICH],
	  CPU_SymbolImbalance[CPU_ENRICH]);
  fprintf(FdBalance, "Stellar evol.  = '%c' / '%c'\n", CPU_Symbol[CPU_ENRICHSTEVOL],
	  CPU_SymbolImbalance[CPU_ENRICHSTEVOL]);
  fprintf(FdBalance, "Snapshot dump  = '%c' / '%c'\n", CPU_Symbol[CPU_SNAPSHOT],
	  CPU_SymbolImbalance[CPU_SNAPSHOT]);
  fprintf(FdBalance, "FoF            = '%c' / '%c'\n", CPU_Symbol[CPU_FOF], CPU_SymbolImbalance[CPU_FOF]);
  fprintf(FdBalance, "Lineofsight    = '%c' / '%c'\n", CPU_Symbol[CPU_LINEOFSIGHT],
	  CPU_SymbolImbalance[CPU_LINEOFSIGHT]);
  fprintf(FdBalance, "Binning        = '%c' / '%c'\n", CPU_Symbol[CPU_BINNING],
	  CPU_SymbolImbalance[CPU_BINNING]);
  fprintf(FdBalance, "Miscellaneous  = '%c' / '%c'\n", CPU_Symbol[CPU_MISC], CPU_SymbolImbalance[CPU_MISC]);
  fprintf(FdBalance, "\n");

#ifdef BG_SFR
  sprintf(buf, "%s%s", All.OutputDir, "sfr.txt");
  if(!(FdSfr = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
/* #ifdef BG_POPIII */
/*   sprintf(buf, "%s%s", All.OutputDir, "popiii_sfr.txt"); */
/*   if(!(FdPOPIIISfr = fopen(buf, mode))) */
/*     { */
/*       printf("error in opening file '%s'\n", buf); */
/*       endrun(1); */
/*     } */
/* #endif */
#endif

#ifdef BG_STELLAR_EVOLUTION
  sprintf(buf, "%s%s", All.OutputDir, "metals_tot.txt");
  if(!(FdMetTot = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  sprintf(buf, "%s%s", All.OutputDir, "metals_gas.txt");
  if(!(FdMetGas = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  sprintf(buf, "%s%s", All.OutputDir, "metals_stars.txt");
  if(!(FdMetStars = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  sprintf(buf, "%s%s", All.OutputDir, "metals_sf.txt");
  if(!(FdMetSF = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  sprintf(buf, "%s%s", All.OutputDir, "SNIa.txt");
  if(!(FdSNIa = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK
  sprintf(buf, "%s%s", All.OutputDir, "kinetic_feedback.txt");
  if(!(FdKineticFeedback = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef BLACK_HOLES
  sprintf(buf, "%s%s", All.OutputDir, "blackholes.txt");
  if(!(FdBlackHoles = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef FORCETEST
  if(RestartFlag == 0)
    {
      sprintf(buf, "%s%s", All.OutputDir, "forcetest.txt");
      if(!(FdForceTest = fopen(buf, "w")))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
      fclose(FdForceTest);
    }
#endif
}



void write_log_files_header(void)
{
#ifdef BG_STELLAR_EVOLUTION
  int i;
#endif

  if(RestartFlag == 0)
    {
#ifdef BG_SNII_THERMAL_FEEDBACK
#ifdef BG_VERBOSE_LOGFILES
      fprintf(FdThermalFeedback, "%12s%14s%14s%13s%13s%13s%13s%13s%13s\n",
	      "#       Time", "StarID", "ParticleID", "E_old", "E_new", "T_old", "T_new", "m_star", "m_gas");
#endif
#endif

      if(ThisTask == 0)
	{
#ifdef BG_STELLAR_EVOLUTION
	  fprintf(FdMetGas, "%12s", "#       Time");
	  fprintf(FdMetStars, "%12s", "#       Time");
	  fprintf(FdMetSF, "%12s", "#       Time");
	  fprintf(FdMetTot, "%12s", "#       Time");

	  for(i = 0; i < BG_NELEMENTS; i++)
	    {
	      fprintf(FdMetGas, "%13s", ElementNames[i]);
	      fprintf(FdMetStars, "%13s", ElementNames[i]);
	      fprintf(FdMetSF, "%13s", ElementNames[i]);
	      fprintf(FdMetTot, "%13s", ElementNames[i]);
	    }

	  fprintf(FdMetGas, "%13s\n", "Total Mass");
	  fprintf(FdMetStars, "%13s\n", "Total Mass");
	  fprintf(FdMetSF, "%13s\n", "Total Mass");
	  fprintf(FdMetTot, "%13s\n", "Total Mass");

	  fprintf(FdSNIa, "%12s%13s\n", "#       Time", "num. SNIa/yr");
#endif

#ifdef BG_SFR
	  fprintf(FdSfr, "%12s%14s%14s%14s%14s\n", "#       Time", "m_* [code]", "sfr [Msun/yr]",
		  "SFR [Msun/yr]", "M_* [code]");
/* #ifdef BG_POPIII */
/* 	  fprintf(FdPOPIIISfr, "%12s%14s%14s%14s%14s\n", "#       Time", "m_* [code]", "sfr [Msun/yr]", */
/* 		  "SFR [Msun/yr]", "M_* [code]"); */
/* #endif */
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK
	  fprintf(FdKineticFeedback, "%12s%13s%13s\n", "#       Time", "ExpMassLoad", "ActMassLoad");
#endif
	}
    }
}




/*!  This function closes the global log-files.
 */
void close_outputfiles(void)
{
#if defined(BLACK_HOLES) && defined(VERBOSE_LOGFILES)
  fclose(FdBlackHolesDetails);	/* needs to be done by everyone */
#endif

#ifdef BG_SNII_THERMAL_FEEDBACK
#ifdef BG_VERBOSE_LOGFILES
  fclose(FdThermalFeedback);	/* needs to be done by everyone */
#endif
#endif

  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);
  fclose(FdBalance);

#ifdef BG_SFR
  fclose(FdSfr);
/* #ifdef BG_POPIII */
/*   fclose(FdPOPIIISfr); */
/* #endif */
#endif

#ifdef BG_STELLAR_EVOLUTION
  fclose(FdMetGas);
  fclose(FdMetStars);
  fclose(FdMetSF);
  fclose(FdMetTot);
  fclose(FdSNIa);
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK
  fclose(FdKineticFeedback);
#endif

#ifdef BLACK_HOLES
  fclose(FdBlackHoles);
#ifdef FOF
  fclose(FdBlackHolesSeeds);
#endif
#endif
}





/*! This function parses the parameterfile in a simple way.  Each paramater is
 *  defined by a keyword (`tag'), and can be either of type douple, int, or
 *  character string.  The routine makes sure that each parameter appears
 *  exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
void read_parameter_file(char *fname)
{
#define REAL 1
#define STRING 2
#define INT 3
#define LABEL 4
#define MAXTAGS 300

  FILE *fd, *fdout;

  char buf[200], buf1[200], buf2[200], buf3[400];

  int i, j, nt, itag;

  int id[MAXTAGS];

  void *addr[MAXTAGS];

  char tag[MAXTAGS][50];

  int pnum, errorFlag = 0;

  char *label, mytag[300];

  All.StarformationOn = 0;	/* defaults */


  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(int) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(float) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(double) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }


  if(ThisTask == 0)		/* read parameter file on process 0 */
    {
      nt = 0;

      strcpy(tag[nt], "InitCondFile");
      addr[nt] = All.InitCondFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "RunLabel");
      addr[nt] = All.RunLabel;
      id[nt++] = LABEL;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnapshotFileBase");
      addr[nt] = All.SnapshotFileBase;
      id[nt++] = STRING;

      strcpy(tag[nt], "EnergyFile");
      addr[nt] = All.EnergyFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "CpuFile");
      addr[nt] = All.CpuFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "InfoFile");
      addr[nt] = All.InfoFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "TimingsFile");
      addr[nt] = All.TimingsFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "RestartFile");
      addr[nt] = All.RestartFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "ResubmitCommand");
      addr[nt] = All.ResubmitCommand;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListFilename");
      addr[nt] = All.OutputListFilename;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListOn");
      addr[nt] = &All.OutputListOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Omega0");
      addr[nt] = &All.Omega0;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaBaryon");
      addr[nt] = &All.OmegaBaryon;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaLambda");
      addr[nt] = &All.OmegaLambda;
      id[nt++] = REAL;

      strcpy(tag[nt], "HubbleParam");
      addr[nt] = &All.HubbleParam;
      id[nt++] = REAL;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &All.BoxSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "PeriodicBoundariesOn");
      addr[nt] = &All.PeriodicBoundariesOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeOfFirstSnapshot");
      addr[nt] = &All.TimeOfFirstSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "CpuTimeBetRestartFile");
      addr[nt] = &All.CpuTimeBetRestartFile;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetStatistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBegin");
      addr[nt] = &All.TimeBegin;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeMax");
      addr[nt] = &All.TimeMax;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetSnapshot");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitLength_in_cm");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitMass_in_g");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = REAL;

      strcpy(tag[nt], "TreeDomainUpdateFrequency");
      addr[nt] = &All.TreeDomainUpdateFrequency;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolIntAccuracy");
      addr[nt] = &All.ErrTolIntAccuracy;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolTheta");
      addr[nt] = &All.ErrTolTheta;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolForceAcc");
      addr[nt] = &All.ErrTolForceAcc;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinGasHsmlFractional");
      addr[nt] = &All.MinGasHsmlFractional;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSizeTimestep");
      addr[nt] = &All.MaxSizeTimestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinSizeTimestep");
      addr[nt] = &All.MinSizeTimestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxRMSDisplacementFac");
      addr[nt] = &All.MaxRMSDisplacementFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "ArtBulkViscConst");
      addr[nt] = &All.ArtBulkViscConst;
      id[nt++] = REAL;

      strcpy(tag[nt], "CourantFac");
      addr[nt] = &All.CourantFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "DesNumNgb");
      addr[nt] = &All.DesNumNgb;
      id[nt++] = INT;

#ifdef SUBFIND
      strcpy(tag[nt], "ErrTolThetaSubfind");
      addr[nt] = &All.ErrTolThetaSubfind;
      id[nt++] = REAL;

      strcpy(tag[nt], "DesLinkNgb");
      addr[nt] = &All.DesLinkNgb;
      id[nt++] = INT;

      if(All.DesLinkNgb > GROUP_MIN_LEN)
	{
	  printf
	    ("Due to the infinite loop problem in subfind_process_group_serial when All.DesLinkNgb > GROUP_MIN_LEN it has proved necessary to set DesLinkNgb equal to GROUP_MIN_LEN - Alan \n");
	  All.DesLinkNgb = GROUP_MIN_LEN;
	  printf("DesLinkNgb = %d\n", All.DesLinkNgb);
	  endrun(10);
	}
#endif

      strcpy(tag[nt], "MaxNumNgbDeviation");
      addr[nt] = &All.MaxNumNgbDeviation;
      id[nt++] = REAL;

      strcpy(tag[nt], "ComovingIntegrationOn");
      addr[nt] = &All.ComovingIntegrationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "ICFormat");
      addr[nt] = &All.ICFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "SnapFormat");
      addr[nt] = &All.SnapFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesPerSnapshot");
      addr[nt] = &All.NumFilesPerSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = INT;

      strcpy(tag[nt], "ResubmitOn");
      addr[nt] = &All.ResubmitOn;
      id[nt++] = INT;

      strcpy(tag[nt], "CoolingOn");
      addr[nt] = &All.CoolingOn;
      id[nt++] = INT;

      strcpy(tag[nt], "StarformationOn");
      addr[nt] = &All.StarformationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "StellarEvolutionOn");
      addr[nt] = &All.StellarEvolutionOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfTimestepCriterion");
      addr[nt] = &All.TypeOfTimestepCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfOpeningCriterion");
      addr[nt] = &All.TypeOfOpeningCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeLimitCPU");
      addr[nt] = &All.TimeLimitCPU;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningHalo");
      addr[nt] = &All.SofteningHalo;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningDisk");
      addr[nt] = &All.SofteningDisk;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBulge");
      addr[nt] = &All.SofteningBulge;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningGas");
      addr[nt] = &All.SofteningGas;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningStars");
      addr[nt] = &All.SofteningStars;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBndry");
      addr[nt] = &All.SofteningBndry;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningHaloMaxPhys");
      addr[nt] = &All.SofteningHaloMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningDiskMaxPhys");
      addr[nt] = &All.SofteningDiskMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBulgeMaxPhys");
      addr[nt] = &All.SofteningBulgeMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningGasMaxPhys");
      addr[nt] = &All.SofteningGasMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningStarsMaxPhys");
      addr[nt] = &All.SofteningStarsMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBndryMaxPhys");
      addr[nt] = &All.SofteningBndryMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "BufferSize");
      addr[nt] = &All.BufferSize;
      id[nt++] = INT;

      strcpy(tag[nt], "PartAllocFactor");
      addr[nt] = &All.PartAllocFactor;
      id[nt++] = REAL;

#ifdef BG_SFR
      strcpy(tag[nt], "StarAllocFactor");
      addr[nt] = &All.StarAllocFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "GasFraction");
      addr[nt] = &All.GasFraction;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "GravityConstantInternal");
      addr[nt] = &All.GravityConstantInternal;
      id[nt++] = REAL;

#if defined(BG_SFR) || defined(BG_COOLING)
      strcpy(tag[nt], "InitGasU_ERG");
      addr[nt] = &All.InitGasU_ERG;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinGasU_ERG");
      addr[nt] = &All.MinGasU_ERG;
      id[nt++] = REAL;
#endif

#ifdef OUTPUTLINEOFSIGHT
      strcpy(tag[nt], "TimeOfFirstLineOfSight");
      addr[nt] = &All.TimeOfFirstLineOfSight;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetLineOfSight");
      addr[nt] = &All.TimeBetLineOfSight;
      id[nt++] = REAL;
#endif

#ifdef ADAPTIVE_OUTPUT
      strcpy(tag[nt], "AO_SubSamplingFactor");
      addr[nt] = &All.AO_SubSamplingFactor;
      id[nt++] = INT;

      strcpy(tag[nt], "AO_TimeOfFirstOutput");
      addr[nt] = &All.AO_TimeOfFirstOutput;
      id[nt++] = REAL;
#endif

#ifdef FOF
      strcpy(tag[nt], "TimeOfFirstFOF");
      addr[nt] = &All.TimeOfFirstFOF;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetFOF");
      addr[nt] = &All.TimeBetFOF;
      id[nt++] = REAL;
#endif

#ifdef BG_OUTPUT_GRID
      strcpy(tag[nt], "TimeOfFirstGridOutput");
      addr[nt] = &All.TimeOfFirstGridOutput;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetGridOutput");
      addr[nt] = &All.TimeBetGridOutput;
      id[nt++] = REAL;
#endif

#ifdef BG_SFR
      strcpy(tag[nt], "Generations");
      addr[nt] = &All.Generations;
      id[nt++] = INT;

      strcpy(tag[nt], "SF_EOSGammaEffective");
      addr[nt] = &All.SF_EOSGammaEffective;
      id[nt++] = REAL;

      strcpy(tag[nt], "SF_EOSMinPhysDens_HpCM3");
      addr[nt] = &All.SF_EOSMinPhysDens_HpCM3;
      id[nt++] = REAL;

      strcpy(tag[nt], "SF_EOSMinOverDens");
      addr[nt] = &All.SF_EOSMinOverDens;
      id[nt++] = REAL;

      strcpy(tag[nt], "SF_EOSEnergyAtThreshold_ERG");
      addr[nt] = &All.SF_EOSEnergyAtThreshold_ERG;
      id[nt++] = REAL;

      strcpy(tag[nt], "SF_SchmidtLawCoeff_MSUNpYRpKPC2");
      addr[nt] = &All.SF_SchmidtLawCoeff_MSUNpYRpKPC2;
      id[nt++] = REAL;
#ifdef BG_DOUBLE_IMF
      strcpy(tag[nt], "SF_SchmidtLawCoeff1_MSUNpYRpKPC2");
      addr[nt] = &All.SF_SchmidtLawCoeff1_MSUNpYRpKPC2;
      id[nt++] = REAL;
#endif
      strcpy(tag[nt], "SF_SchmidtLawExponent");
      addr[nt] = &All.SF_SchmidtLawExponent;
      id[nt++] = REAL;

      strcpy(tag[nt], "SF_EOSTempThreshMargin_DEX");
      addr[nt] = &All.SF_EOSTempThreshMargin_DEX;
      id[nt++] = REAL;
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK
      strcpy(tag[nt], "SNII_WindOn");
      addr[nt] = &All.SNII_WindOn;
      id[nt++] = INT;

      strcpy(tag[nt], "SNII_WindIsotropicOn");
      addr[nt] = &All.SNII_WindIsotropicOn;
      id[nt++] = INT;

      strcpy(tag[nt], "SNII_WindSpeed_KMpS");
      addr[nt] = &All.SNII_WindSpeed_KMpS;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_WindMassLoading");
      addr[nt] = &All.SNII_WindMassLoading;
      id[nt++] = REAL;
#ifdef BG_DOUBLE_IMF
      strcpy(tag[nt], "SNII_WindMassLoading1");
      addr[nt] = &All.SNII_WindMassLoading1;
      id[nt++] = REAL;
#endif
#ifdef BG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING
      strcpy(tag[nt], "SNII_WindDecouplingTime_YR");
      addr[nt] = &All.SNII_WindDecouplingTime_YR;
      id[nt++] = REAL;
#endif
#endif

#if defined(BG_SNII_KINETIC_FEEDBACK) || defined(BG_SNII_THERMAL_FEEDBACK)
      strcpy(tag[nt], "SNII_WindDelay_YR");
      addr[nt] = &All.SNII_WindDelay_YR;
      id[nt++] = REAL;
#endif

#if defined(BG_STELLAR_EVOLUTION) || defined(BG_COOLING)
      strcpy(tag[nt], "InitAbundance_Hydrogen");
      addr[nt] = &All.InitAbundance_Hydrogen;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_Helium");
      addr[nt] = &All.InitAbundance_Helium;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_Carbon");
      addr[nt] = &All.InitAbundance_Carbon;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_Nitrogen");
      addr[nt] = &All.InitAbundance_Nitrogen;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_Oxygen");
      addr[nt] = &All.InitAbundance_Oxygen;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_Neon");
      addr[nt] = &All.InitAbundance_Neon;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_Magnesium");
      addr[nt] = &All.InitAbundance_Magnesium;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_Silicon");
      addr[nt] = &All.InitAbundance_Silicon;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_Iron");
      addr[nt] = &All.InitAbundance_Iron;
      id[nt++] = REAL;

      strcpy(tag[nt], "CalciumOverSilicon");
      addr[nt] = &All.CalciumOverSilicon;
      id[nt++] = REAL;

      strcpy(tag[nt], "SulphurOverSilicon");
      addr[nt] = &All.SulphurOverSilicon;
      id[nt++] = REAL;
#endif

#ifdef BG_COOLING
      strcpy(tag[nt], "REION_H_ZCenter");
      addr[nt] = &All.REION_H_ZCenter;
      id[nt++] = REAL;

      strcpy(tag[nt], "REION_H_ZSigma");
      addr[nt] = &All.REION_H_ZSigma;
      id[nt++] = REAL;

      strcpy(tag[nt], "REION_H_Heating_EVpH");
      addr[nt] = &All.REION_H_Heating_EVpH;
      id[nt++] = REAL;

      strcpy(tag[nt], "REION_He_ZCenter");
      addr[nt] = &All.REION_He_ZCenter;
      id[nt++] = REAL;

      strcpy(tag[nt], "REION_He_ZSigma");
      addr[nt] = &All.REION_He_ZSigma;
      id[nt++] = REAL;

      strcpy(tag[nt], "REION_He_Heating_EVpH");
      addr[nt] = &All.REION_He_Heating_EVpH;
      id[nt++] = REAL;

      strcpy(tag[nt], "CoolTablePath");
      addr[nt] = &All.CoolTablePath;
      id[nt++] = STRING;

      strcpy(tag[nt], "MetDepCoolingOn");
      addr[nt] = &All.MetDepCoolingOn;
      id[nt++] = INT;

#ifdef BG_COOLING_SHIELDING
      strcpy(tag[nt], "UVBG_PhysDensThresh_HpCM3");
      addr[nt] = &All.UVBG_PhysDensThresh_HpCM3;
      id[nt++] = REAL;
#endif
#endif

#ifdef BG_STELLAR_EVOLUTION
      strcpy(tag[nt], "IMF_Model");
      addr[nt] = &All.IMF_Model;
      id[nt++] = STRING;

      strcpy(tag[nt], "IMF_LifetimeModel");
      addr[nt] = &All.IMF_LifetimeModel;
      id[nt++] = STRING;

      strcpy(tag[nt], "IMF_Exponent");
      addr[nt] = &All.IMF_Exponent;
      id[nt++] = REAL;
#ifdef BG_DOUBLE_IMF
      strcpy(tag[nt], "IMF_Exponent1");
      addr[nt] = &All.IMF_Exponent1;
      id[nt++] = REAL;

      strcpy(tag[nt], "IMF_PhysDensThresh_HpCM3");
      addr[nt] = &All.IMF_PhysDensThresh_HpCM3;
      id[nt++] = REAL;
#endif
      strcpy(tag[nt], "IMF_MinMass_MSUN");
      addr[nt] = &All.IMF_MinMass_MSUN;
      id[nt++] = REAL;

      strcpy(tag[nt], "IMF_MaxMass_MSUN");
      addr[nt] = &All.IMF_MaxMass_MSUN;
      id[nt++] = REAL;

      strcpy(tag[nt], "AGB_MassTransferOn");
      addr[nt] = &All.AGB_MassTransferOn;
      id[nt++] = INT;

      strcpy(tag[nt], "SNIa_Model");
      addr[nt] = &All.SNIa_Model;
      id[nt++] = STRING;

      strcpy(tag[nt], "SNIa_Efficiency_fracwd");
      addr[nt] = &All.SNIa_Efficiency_fracwd;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNIa_MassTransferOn");
      addr[nt] = &All.SNIa_MassTransferOn;
      id[nt++] = INT;

      strcpy(tag[nt], "SNIa_EnergyTransferOn");
      addr[nt] = &All.SNIa_EnergyTransferOn;
      id[nt++] = INT;

      strcpy(tag[nt], "SNIa_Energy_ERG");
      addr[nt] = &All.SNIa_Energy_ERG;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_MassTransferOn");
      addr[nt] = &All.SNII_MassTransferOn;
      id[nt++] = INT;

      strcpy(tag[nt], "SNII_Factor_Hydrogen");
      addr[nt] = &All.SNII_Factor_Hydrogen;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_Factor_Helium");
      addr[nt] = &All.SNII_Factor_Helium;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_Factor_Carbon");
      addr[nt] = &All.SNII_Factor_Carbon;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_Factor_Nitrogen");
      addr[nt] = &All.SNII_Factor_Nitrogen;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_Factor_Oxygen");
      addr[nt] = &All.SNII_Factor_Oxygen;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_Factor_Neon");
      addr[nt] = &All.SNII_Factor_Neon;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_Factor_Magnesium");
      addr[nt] = &All.SNII_Factor_Magnesium;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_Factor_Silicon");
      addr[nt] = &All.SNII_Factor_Silicon;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_Factor_Iron");
      addr[nt] = &All.SNII_Factor_Iron;
      id[nt++] = REAL;

      strcpy(tag[nt], "YieldTablePath");
      addr[nt] = &All.YieldTablePath;
      id[nt++] = STRING;
#endif

#ifdef BG_SFR
#ifdef BG_SNII_THERMAL_FEEDBACK
      strcpy(tag[nt], "SNII_EnergyTransferOn");
      addr[nt] = &All.SNII_EnergyTransferOn;
      id[nt++] = INT;

      strcpy(tag[nt], "SNII_EnergyPerUnitMass_ERGperG");
      addr[nt] = &All.SNII_EnergyPerUnitMass_ERGperG;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_EnergyFraction");
      addr[nt] = &All.SNII_EnergyFraction;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_Energy_ERG");
      addr[nt] = &All.SNII_Energy_ERG;
      id[nt++] = REAL;
#endif

#if defined(BG_STELLAR_EVOLUTION) || defined(BG_SNII_THERMAL_FEEDBACK)
      strcpy(tag[nt], "SNII_MinMass_MSUN");
      addr[nt] = &All.SNII_MinMass_MSUN;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_MaxMass_MSUN");
      addr[nt] = &All.SNII_MaxMass_MSUN;
      id[nt++] = REAL;
#endif

/* #if defined(BG_POPIII) && defined(BG_POPIII_THERMAL_FEEDBACK) */
/*       strcpy(tag[nt], "POPIII_EnergyTransferOn"); */
/*       addr[nt] = &All.POPIII_EnergyTransferOn; */
/*       id[nt++] = INT; */

/*       strcpy(tag[nt], "POPIII_EnergyPerUnitMass_ERGperG"); */
/*       addr[nt] = &All.POPIII_EnergyPerUnitMass_ERGperG; */
/*       id[nt++] = REAL; */

/*       strcpy(tag[nt], "POPIII_EnergyFraction"); */
/*       addr[nt] = &All.POPIII_EnergyFraction; */
/*       id[nt++] = REAL; */

/*       strcpy(tag[nt], "POPIII_SNII_Energy_ERG"); */
/*       addr[nt] = &All.POPIII_SNII_Energy_ERG; */
/*       id[nt++] = REAL; */

/*       strcpy(tag[nt], "POPIII_PISN_Energy_ERG"); */
/*       addr[nt] = &All.POPIII_PISN_Energy_ERG; */
/*       id[nt++] = REAL; */

/*       strcpy(tag[nt], "POPIII_SNII_MinMass_MSUN"); */
/*       addr[nt] = &All.POPIII_SNII_MinMass_MSUN; */
/*       id[nt++] = REAL; */

/*       strcpy(tag[nt], "POPIII_SNII_MaxMass_MSUN"); */
/*       addr[nt] = &All.POPIII_SNII_MaxMass_MSUN; */
/*       id[nt++] = REAL; */

/*       strcpy(tag[nt], "POPIII_PISN_MinMass_MSUN"); */
/*       addr[nt] = &All.POPIII_PISN_MinMass_MSUN; */
/*       id[nt++] = REAL; */

/*       strcpy(tag[nt], "POPIII_PISN_MaxMass_MSUN"); */
/*       addr[nt] = &All.POPIII_PISN_MaxMass_MSUN; */
/*       id[nt++] = REAL; */
/* #endif */
#endif /* BG_SFR */


#if defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS_HSML)
      strcpy(tag[nt], "ReferenceGasMass");
      addr[nt] = &All.ReferenceGasMass;
      id[nt++] = REAL;
#endif


#ifdef BLACK_HOLES
      strcpy(tag[nt], "SeedBHMassOverGasMass");
      addr[nt] = &All.SeedBHMassOverGasMass;
      id[nt++] = REAL;

#ifdef FOF
      strcpy(tag[nt], "MinFoFSizeForNewSeed");
      addr[nt] = &All.MinFoFSizeForNewSeed;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "BlackHoleAccretionFactor");
      addr[nt] = &All.BlackHoleAccretionFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleAccretionSlope");
      addr[nt] = &All.BlackHoleAccretionSlope;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleFeedbackFactor");
      addr[nt] = &All.BlackHoleFeedbackFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleRadiativeEfficiency");
      addr[nt] = &All.BlackHoleRadiativeEfficiency;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleEddingtonFactor");
      addr[nt] = &All.BlackHoleEddingtonFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleMinHeatTemp");
      addr[nt] = &All.BlackHoleMinHeatTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleNumberOfNeighboursToHeat");
      addr[nt] = &All.BlackHoleNumberOfNeighboursToHeat;
      id[nt++] = REAL;

#endif


#ifdef BG_SFR
      strcpy(tag[nt], "SF_THRESH_MaxPhysDensOn");
      addr[nt] = &All.SF_THRESH_MaxPhysDensOn;
      id[nt++] = INT;

      strcpy(tag[nt], "SF_THRESH_MinOverDens");
      addr[nt] = &All.SF_THRESH_MinOverDens;
      id[nt++] = REAL;

      strcpy(tag[nt], "SF_THRESH_MinPhysDens_HpCM3");
      addr[nt] = &All.SF_THRESH_MinPhysDens_HpCM3;
      id[nt++] = REAL;

      strcpy(tag[nt], "SF_THRESH_MaxPhysDens_HpCM3");
      addr[nt] = &All.SF_THRESH_MaxPhysDens_HpCM3;
      id[nt++] = REAL;

      strcpy(tag[nt], "SF_THRESH_MaxTemp_K");
      addr[nt] = &All.SF_THRESH_MaxTemp_K;
      id[nt++] = REAL;

/* #ifdef BG_POPIII */
/*       strcpy(tag[nt], "POPIII_MassTransferOn"); */
/*       addr[nt] = &All.POPIII_MassTransferOn; */
/*       id[nt++] = INT; */

/*       strcpy(tag[nt], "POPIII_MetallicityThreshold_SOLAR"); */
/*       addr[nt] = &All.POPIII_MetallicityThreshold_SOLAR; */
/*       id[nt++] = REAL; */

/*       strcpy(tag[nt], "POPIII_IMF_Exponent"); */
/*       addr[nt] = &All.POPIII_IMF_Exponent; */
/*       id[nt++] = REAL; */

/*       strcpy(tag[nt], "POPIII_IMF_MinMass_MSUN"); */
/*       addr[nt] = &All.POPIII_IMF_MinMass_MSUN; */
/*       id[nt++] = REAL; */

/*       strcpy(tag[nt], "POPIII_IMF_MaxMass_MSUN"); */
/*       addr[nt] = &All.POPIII_IMF_MaxMass_MSUN; */
/*       id[nt++] = REAL; */

/*       strcpy(tag[nt], "POPIII_MinMass_MSUN"); */
/*       addr[nt] = &All.POPIII_MinMass_MSUN; */
/*       id[nt++] = REAL; */

/*       strcpy(tag[nt], "POPIII_MaxMass_MSUN"); */
/*       addr[nt] = &All.POPIII_MaxMass_MSUN; */
/*       id[nt++] = REAL; */
/* #endif */
#endif


#if defined(BG_MOL_COOLING) || defined(BG_MOL_NETWORK)
      strcpy(tag[nt], "InitAbundance_Hp");
      addr[nt] = &All.InitAbundance_Hp;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_Hm");
      addr[nt] = &All.InitAbundance_Hm;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_H2");
      addr[nt] = &All.InitAbundance_H2;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_H2p");
      addr[nt] = &All.InitAbundance_H2p;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_Hep");
      addr[nt] = &All.InitAbundance_Hep;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_Hepp");
      addr[nt] = &All.InitAbundance_Hepp;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_D");
      addr[nt] = &All.InitAbundance_D;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_Dp");
      addr[nt] = &All.InitAbundance_Dp;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitAbundance_HD");
      addr[nt] = &All.InitAbundance_HD;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitNumberDensity_e");
      addr[nt] = &All.InitNumberDensity_e;
      id[nt++] = REAL;
#endif

/* #ifdef BG_DUST */
/* #if defined(BG_DUST_DESTRUCTION_SUBLIMATION) || defined(BG_DUST_DESTRUCTION_SPUTTERING) */
/*       strcpy(tag[nt], "InitDustGrainSize_nM"); */
/*       addr[nt] = &All.InitDustGrainSize_nM; */
/*       id[nt++] = REAL; */
/* #endif */
/* #endif */

#ifdef DARKENERGY
      strcpy(tag[nt], "DarkEnergyParam");
      addr[nt] = &All.DarkEnergyParam;
      id[nt++] = REAL;
#endif

#ifdef RESCALEVINI
      strcpy(tag[nt], "VelIniScale");
      addr[nt] = &All.VelIniScale;
      id[nt++] = REAL;
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyFile");
      addr[nt] = All.DarkEnergyFile;
      id[nt++] = STRING;
#endif
#endif

#ifdef TIME_DEP_ART_VISC
      strcpy(tag[nt], "ViscositySourceScaling");
      addr[nt] = &All.ViscSource0;
      id[nt++] = REAL;

      strcpy(tag[nt], "ViscosityDecayLength");
      addr[nt] = &All.DecayLength;
      id[nt++] = REAL;

      strcpy(tag[nt], "ViscosityAlphaMin");
      addr[nt] = &All.AlphaMin;
      id[nt++] = REAL;
#endif


      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	      while(!feof(fd))
		{
		  char *ret;

		  *buf = 0;
		  ret = fgets(buf, 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;

		  if(buf1[0] == '%')
		    continue;

		  for(i = 0, j = -1; i < nt; i++)
		    if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			strcpy(mytag, tag[j]);
			tag[i][0] = 0;
			break;
		      }

		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case REAL:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy((char *) addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
			case LABEL:
			  label = buf;
			  /* remove name of tag */
			  for(itag = 0; itag < strlen(mytag); itag++)
			    label++;
			  /* remove leading blanks */
			  while(*label == ' ')
			    label++;
			  strcpy((char *) addr[j], label);
			  fprintf(fdout, "%s\n", buf);
			  break;
			}
		    }
		  else
		    {
		      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			      fname, buf1);
		      errorFlag = 1;
		    }
		}
	      fclose(fd);
	      fclose(fdout);

	      i = strlen(All.OutputDir);
	      if(i > 0)
		if(All.OutputDir[i - 1] != '/')
		  strcat(All.OutputDir, "/");

	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
	      sprintf(buf3, "cp %s %s", buf1, buf2);
#ifndef NOCALLSOFSYSTEM
	      int ret;

	      ret = system(buf3);
#endif
	    }
	}
      else
	{
	  printf("Parameter file %s not found.\n", fname);
	  errorFlag = 1;
	}


      for(i = 0; i < nt; i++)
	{
	  if(*tag[i])
	    {
	      printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	      errorFlag = 1;
	    }
	}

      if(All.OutputListOn && errorFlag == 0)
	errorFlag += read_outputlist(All.OutputListFilename);
      else
	All.OutputListLength = 0;
    }

  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(errorFlag)
    {
      MPI_Finalize();
      exit(0);
    }

  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);



  for(pnum = 0; All.NumFilesWrittenInParallel > (1 << pnum); pnum++);

  if(All.NumFilesWrittenInParallel != (1 << pnum))
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be a power of 2\n");
      endrun(0);
    }

  if(All.NumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
      endrun(0);
    }

  if(NTask < All.NumFilesPerSnapshot)
    {
      if(ThisTask == 0)
	printf("Fatal error.\nNumber of processors must be larger or equal than All.NumFilesPerSnapshot.\n");
      endrun(0);
    }
  if(All.SnapFormat < 1 || All.SnapFormat > 3)
    {
      if(ThisTask == 0)
	printf("Unsupported File-Format\n");
      endrun(0);
    }
#ifndef  HAVE_HDF5
  if(All.SnapFormat == 3)
    {
      if(ThisTask == 0)
	printf("Code wasn't compiled with HDF5 support enabled!\n");
      endrun(0);
    }
#endif

#ifdef PERIODIC
  if(All.PeriodicBoundariesOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched on.\n");
	  printf("You must set `PeriodicBoundariesOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.PeriodicBoundariesOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched off.\n");
	  printf("You must set `PeriodicBoundariesOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif

#ifdef BG_COOLING
  if(All.CoolingOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with cooling switched on.\n");
	  printf("You must set `CoolingOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.CoolingOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with cooling switched off.\n");
	  printf("You must set `CoolingOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif

  if(All.TypeOfTimestepCriterion >= 3)
    {
      if(ThisTask == 0)
	{
	  printf("The specified timestep criterion\n");
	  printf("is not valid\n");
	}
      endrun(0);
    }

#if defined(LONG_X) ||  defined(LONG_Y) || defined(LONG_Z)
#ifndef NOGRAVITY
  if(ThisTask == 0)
    {
      printf("Code was compiled with LONG_X/Y/Z, but not with NOGRAVITY.\n");
      printf("Stretched periodic boxes are not implemented for gravity yet.\n");
    }
  endrun(0);
#endif
#endif

#ifdef BG_SFR
  if(All.StarformationOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with star formation switched on.\n");
	  printf("You must set `StarformationOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
  /*
  if(All.CoolingOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("You try to use the code with star formation enabled,\n");
	  printf("but you did not switch on cooling.\nThis mode is not supported.\n");
	}
      endrun(0);
    }
  */
#else
  if(All.StarformationOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with star formation switched off.\n");
	  printf("You must set `StarformationOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif

#ifdef BG_SFR
#ifdef BG_STELLAR_EVOLUTION
  if(All.StellarEvolutionOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with stellar evolution switched on.\n");
	  printf("You must set `StellarEvolutionOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.StellarEvolutionOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with stellar evolution switched off.\n");
	  printf("You must set `StellarEvolutionOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif
#endif

#ifdef TIMEDEPDE
#ifndef DARKENERGY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIMEDEPDE, but not with DARKENERGY.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif


#ifndef BG_SFR
#ifdef BG_STELLAR_EVOLUTION
  BG_STELLAR_EVOLUTION compilation flag can be switched on only when BG_SFR is!
#endif
#ifdef BG_SNII_KINETIC_FEEDBACK
  BG_SNII_KINETIC_FEEDBACK compilation falg can be switched on only when BG_SFR is!
#endif
#ifdef BG_SNII_THERMAL_FEEDBACK
  BG_SNII_THERMAL_FEEDBACK compilation falg can be switched on only when BG_SF is!
#endif
#endif
#ifdef BG_SFR
#if defined(BG_SNII_KINETIC_FEEDBACK) && defined(BG_SNII_THERMAL_FEEDBACK)
  BG_SNII_KINETIC_FEEDBACK and BG_SNII_THERMAL_FEEDBACK cannot be on at the same time!
#endif
#ifndef BG_SNII_KINETIC_FEEDBACK
#ifdef BG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING
  BG_SNII_THERMAL_FEEDBACK_SF_DECOUPLING compilation falg can be switched on only when BG_SNII_KINETIC_FEEDBACK is!
#endif
#endif
#endif
#if defined(BG_SFR) && !defined(BG_STELLAR_EVOLUTION)
#ifdef BG_SNIA_IRON
  BG_SNIA_IRON compilation falg can be switched on only when BG_STELLAR_EVOLUTION is!
#endif
#ifdef BG_DOUBLE_IMF
  BG_DOUBLE_IMF compilation falg can be switched on only when BG_STELLAR_EVOLUTION is!
#endif
#ifdef BG_Z_WEIGHTED_REDSHIFT
  BG_Z_WEIGHTED_REDSHIFT compilation falg can be switched on only when BG_STELLAR_EVOLUTION is!
#endif
#ifdef BG_METALSMOOTHING
  BG_METALSMOOTHING compilation falg can be switched on only when BG_STELLAR_EVOLUTION is!
#endif
#ifdef BG_OUTPUT_GRID
  BG_OUTPUT_GRID compilation falg can be switched on only when BG_STELLAR_EVOLUTION is!
#endif
#endif


#undef REAL
#undef STRING
#undef INT
#undef LABEL
#undef MAXTAGS
}


/*! this function reads a table with a list of desired output times. The table
 *  does not have to be ordered in any way, but may not contain more than
 *  MAXLEN_OUTPUTLIST entries.
 */
int read_outputlist(char *fname)
{
  FILE *fd;

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      return 1;
    }

  All.OutputListLength = 0;
  do
    {
      if(fscanf(fd, " %lg ", &All.OutputListTimes[All.OutputListLength]) == 1)
	All.OutputListLength++;
      else
	break;
    }
  while(All.OutputListLength < MAXLEN_OUTPUTLIST);

  fclose(fd);

  printf("\nfound %d times in output-list.\n", All.OutputListLength);

  return 0;
}


/*! If a restart from restart-files is carried out where the TimeMax variable
 * is increased, then the integer timeline needs to be adjusted. The approach
 * taken here is to reduce the resolution of the integer timeline by factors
 * of 2 until the new final time can be reached within TIMEBASE.
 */
void readjust_timebase(double TimeMax_old, double TimeMax_new)
{
  int i;

  long long ti_end;

  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType 'long long' is not 64 bit on this platform\n\n");
      endrun(555);
    }

  if(ThisTask == 0)
    {
      printf("\nAll.TimeMax has been changed in the parameterfile\n");
      printf("Need to adjust integer timeline\n\n\n");
    }

  if(TimeMax_new < TimeMax_old)
    {
      if(ThisTask == 0)
	printf("\nIt is not allowed to reduce All.TimeMax\n\n");
      endrun(556);
    }

  if(All.ComovingIntegrationOn)
    ti_end = (long long) (log(TimeMax_new / All.TimeBegin) / All.Timebase_interval);
  else
    ti_end = (long long) ((TimeMax_new - All.TimeBegin) / All.Timebase_interval);

  while(ti_end > TIMEBASE)
    {
      All.Timebase_interval *= 2.0;

      ti_end /= 2;
      All.Ti_Current /= 2;

#ifdef PMGRID
      All.PM_Ti_begstep /= 2;
      All.PM_Ti_endstep /= 2;
#endif

      for(i = 0; i < NumPart; i++)
	{
	  P[i].Ti_begstep /= 2;
	  P[i].TimeBin--;
	  if(P[i].TimeBin == 0)
	    {
	      printf("Error in readjust_timebase(). Minimum Timebin for particle %d reached.\n", i);
	      endrun(8765);
	    }
	}

      All.Ti_nextlineofsight /= 2;
    }

  All.TimeMax = TimeMax_new;
}
