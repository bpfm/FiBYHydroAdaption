#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#include "forcetree.h"


#if defined(BG_SFR) || defined(BG_COOLING)

#include "bg_cooling.h"
#include "bg_proto.h"
#include "bg_vars.h"

/*
 * This routine does cooling and star formation for
 * the STELLA project on BlueGene.
 */

void cooling_and_starformation(void)
/* cooling routine when star formation is enabled */
{
  int i;

  double dt, dtime, dtime_in_s, ascale = 1, hubble_a = 0, a3inv;

  double time_hubble_a, tstart, tend, timecooling = 0, timesfr = 0;

  double redshift;

#if defined(BG_SFR) && defined(BG_COOLING)
  double egycurrent;
#endif

#ifdef BG_COOLING
  int j, z_index;

  float d_z;

  double dmax1, dmax2;

  double rho, dz, unew, new_entropy, element_metallicity[BG_NELEMENTS];
#endif

#ifdef BG_SFR
  int stars_spawned, tot_spawned, stars_converted, tot_converted, number_of_stars_generated, bin;

  unsigned int bits;

  double sfrrate, totsfrrate, rate_in_msunperyear, pressure;

  double sm, sum_sm, total_sm, rate, sum_mass_stars, total_sum_mass_stars;
#ifdef BG_POPIII
  double popiii_sfrrate, popiii_totsfrrate, popiii_sum_sm, popiii_total_sm;
  double popiii_sum_mass_stars, popiii_total_sum_mass_stars;
#endif
  double prob, mass_of_star, frac;
#if defined(BG_COOLING) && defined(BG_SFR)
 double effective_entropy;
#endif

#if defined(BG_SNII_KINETIC_FEEDBACK) && defined(BG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING)
  double time_in_wind_in_Gyr;
#endif

#ifdef BG_STELLAR_EVOLUTION
  int k;
#endif
#endif /* BG_SFR */

#if defined(BH_THERMALFEEDBACK)
  double temp, u_to_temp_fac;

  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1
    * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#endif

/* #if defined(FLTROUNDOFFREDUCTION) && ((defined(BLACK_HOLES) && defined(BH_THERMALFEEDBACK)) || defined(BG_SNII_THERMAL_FEEDBACK)) */
/*   for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) */
/*     if(P[i].Type == 0) */
/*       SphP[i].i.Injected_Energy = FLT(SphP[i].i.dInjected_Energy); */
/* #endif */

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time)
#ifdef DARKENERGY
	+ DarkEnergy_a(All.Time);
#else
	+ All.OmegaLambda;
#endif
      hubble_a = All.Hubble * sqrt(hubble_a);

      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;
      redshift = 1 / ascale - 1;
    }
  else
    {
      a3inv = ascale = time_hubble_a = 1;
      redshift = 0;
    }

#ifdef BG_SFR
  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinActive[bin])
      TimeBinSfr[bin] = 0;
#ifdef BG_POPIII
  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinActive[bin])
      TimeBinPOPIIISfr[bin] = 0;
#endif

  for(bits = 0; All.Generations > (1 << bits); bits++);

  stars_spawned = stars_converted = 0;
  sum_sm = sum_mass_stars = 0;
#ifdef BG_POPIII
  popiii_sum_sm = popiii_sum_mass_stars = 0;
#endif
#endif

#ifdef BG_COOLING
  get_redshift_index(redshift, &z_index, &d_z);
  LoadCoolingTables(z_index);
#endif

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
	dt = TISTEP(P[i].TimeBin) * All.Timebase_interval;

	/* the actual time-step */
	if(All.ComovingIntegrationOn)
	  dtime = All.Time * dt / time_hubble_a;
	else
	  dtime = dt;


	tstart = second();

#ifdef BG_COOLING


	/* ======================================================================
	 *  COOLING
	 * ====================================================================== */


	/* compute element metallicity as element mass fractions */
	if(All.MetDepCoolingOn == 1)
	  {
#ifdef BG_STELLAR_EVOLUTION
#ifdef BG_METALSMOOTHING
	    for(j = 0; j < BG_NELEMENTS; j++)
	      element_metallicity[j] = SphP[i].MetalsSmoothed[j] / P[i].Mass;
#else
	    for(j = 0; j < BG_NELEMENTS; j++)
	      element_metallicity[j] = SphP[i].Metals[j] / P[i].Mass;
#endif

#ifdef BG_DUST_METAL_COOLING_CORRECTION
#ifdef BG_METALSMOOTHING
            element_metallicity[element_index("Silicon")] -= 0.2426 * SphP[i].DusticitySmoothed;
            element_metallicity[element_index("Oxygen")] -= 0.5978 * SphP[i].DusticitySmoothed;
            element_metallicity[element_index("Magnesium")] -= 0.1028 * SphP[i].DusticitySmoothed;
            element_metallicity[element_index("Iron")] -= 0.056 * SphP[i].DusticitySmoothed;
#else
            element_metallicity[element_index("Silicon")] -= 0.2426 * SphP[i].Dusticity;
            element_metallicity[element_index("Oxygen")] -= 0.5978 * SphP[i].Dusticity;
            element_metallicity[element_index("Magnesium")] -= 0.1028 * SphP[i].Dusticity;
            element_metallicity[element_index("Iron")] -= 0.056 * SphP[i].Dusticity;
#endif

            if(element_metallicity[element_index("Silicon")] < 0)
              element_metallicity[element_index("Silicon")] = 0;
            if(element_metallicity[element_index("Oxygen")] < 0)
              element_metallicity[element_index("Oxygen")] = 0;
            if(element_metallicity[element_index("Magnesium")] < 0)
              element_metallicity[element_index("Magnesium")] = 0;
            if(element_metallicity[element_index("Iron")] < 0)
              element_metallicity[element_index("Iron")] = 0;
#endif

#else /* BG_STELLAR_EVOLUTION */
	    element_metallicity[element_index("Hydrogen")] = All.InitAbundance_Hydrogen;
	    element_metallicity[element_index("Helium")] = All.InitAbundance_Helium;
	    element_metallicity[element_index("Carbon")] = All.InitAbundance_Carbon;
	    element_metallicity[element_index("Nitrogen")] = All.InitAbundance_Nitrogen;
	    element_metallicity[element_index("Oxygen")] = All.InitAbundance_Oxygen;
	    element_metallicity[element_index("Neon")] = All.InitAbundance_Neon;
	    element_metallicity[element_index("Magnesium")] = All.InitAbundance_Magnesium;
	    element_metallicity[element_index("Silicon")] = All.InitAbundance_Silicon;
	    element_metallicity[element_index("Iron")] = All.InitAbundance_Iron;
#endif /* BG_STELLAR_EVOLUTION */
	  }
	else
	  {
	    for(j = 0; j < BG_NELEMENTS; j++)
	      if(j == element_index("Hydrogen"))
		element_metallicity[j] = All.InitAbundance_Hydrogen;
	      else if(j == element_index("Helium"))
		element_metallicity[j] = All.InitAbundance_Helium;
	      else
		element_metallicity[j] = 0;
	  }

	/* get the internal energy per unit mass */
	unew = (SphP[i].Entropy + SphP[i].e.DtEntropy * dt) /
	  GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);

	/* check against minimum allowed energy */
	unew = DMAX(All.MinEgySpec, unew);

	/*                                                                   */
	/*                                                                   */
	/*                                                                   */
	/* SHALL WE HAVE A CHECK AGAINST SOME MAXIMUM ALLOWED ENERGY HERE??? */
	/*                                                                   */
	/*                                                                   */
	/*                                                                   */

	/* convert to CGS units for the cooling routine */
	unew *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

	rho = SphP[i].d.Density * a3inv;
	rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

	dtime_in_s = dtime * All.UnitTime_in_s / All.HubbleParam;

	dz = 1 / All.Time - 1 / exp(log(All.Time) - dt);

	/*
        if(bg_get_temperature(i) > 1e6)
          printf("[TASK %2d] before cooling: PID=%d, temp=%g, entr=%g, dentr=%g\n",
		 ThisTask, P[i].ID, bg_get_temperature(i), SphP[i].Entropy, SphP[i].e.DtEntropy);
	*/

	/* do the cooling */
#ifdef BG_MOL_COOLING
	double x_H2, x_HD;

#ifdef BG_METALSMOOTHING
	x_H2 = SphP[i].x_H2_Smoothed + SphP[i].x_H2p_Smoothed;
	x_HD = SphP[i].x_HD_Smoothed;
#else
	x_H2 = SphP[i].x_H2 + SphP[i].x_H2p;
	x_HD = SphP[i].x_HD;
#endif
	unew = DoCooling(unew, rho, dtime_in_s, dz, redshift, d_z, element_metallicity, x_H2, x_HD);
#else
	unew = DoCooling(unew, rho, dtime_in_s, dz, redshift, d_z, element_metallicity);
#endif

	/* convert energy back to code units */
	unew *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs;

	/* check against minimum allowed energy */
	unew = DMAX(All.MinEgySpec, unew);

/* #if defined(BH_THERMALFEEDBACK) || defined(BG_SNII_THERMAL_FEEDBACK) */
/* 	/\* add thermal feedback internal energy if needed *\/ */
/* 	if(SphP[i].i.Injected_Energy) */
/* 	  { */
/* 	    if(P[i].Mass == 0) */
/* 	      SphP[i].i.Injected_Energy = 0; */
/* 	    else */
/* 	      { */
/* 		unew += SphP[i].i.Injected_Energy / P[i].Mass; */

/* 		/\* NEW *\/ */
/* 		SphP[i].Entropy = unew * GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1); */
/* 		SphP[i].e.DtEntropy = 0; */
/* 		/\* NEW *\/ */
/* 	      } */
/* #ifdef FLTROUNDOFFREDUCTION */
/* 	    SphP[i].i.dInjected_Energy = 0; */
/* #else */
/* 	    SphP[i].i.Injected_Energy = 0; */
/* #endif */
/* 	  } */
/* 	else */
/* 	  { */
/* #endif */

	if(P[i].TimeBin)	/* upon start-up, we need to protect against dt==0 */
	  {
	    /* note: the adiabatic rate has been already added in! */
	    if(dt > 0)
	      {
		new_entropy = unew * GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);
#if defined(BG_SFR)
		/* prevent particles to cool below the effective EoS */
/* 		if(SphP[i].d.Density * a3inv >= All.EOSPhysDensThresh && */
/* 		   (SphP[i].d.Density > All.OverDensThresh || All.ComovingIntegrationOn == 0)) */
/* 		  { */
 		if(SphP[i].d.Density * a3inv >= All.EOSPhysDensThresh &&
		   (SphP[i].d.Density >= All.EOSOverDensThresh || All.ComovingIntegrationOn == 0))
		  {
		    /* physical energy in cgs units */
		    egycurrent = All.SF_EOSEnergyAtThreshold_ERG * All.UnitMass_in_g / All.UnitEnergy_in_cgs;

		    /* physical pressure */
		    pressure = GAMMA_MINUS1 * All.EOSPhysDensThresh * egycurrent *
		      pow(SphP[i].d.Density * a3inv / All.EOSPhysDensThresh, All.SF_EOSGammaEffective);

		    /* convert pressure to comoving pressure */
		    pressure /= pow(a3inv, GAMMA);

		    effective_entropy = pressure / pow(SphP[i].d.Density, GAMMA);

		    new_entropy = max(new_entropy, effective_entropy);
		  }
/* 		  } */
#endif /* BG_SFR */
		SphP[i].e.DtEntropy = (new_entropy - SphP[i].Entropy) / dt;
	      }
	  }

/* #if defined(BH_THERMALFEEDBACK) || defined(BG_SNII_THERMAL_FEEDBACK) */
/* 	  } */
/* #endif */


	/*
        if(bg_get_temperature(i) > 1e6)
          printf("[TASK %2d]  after cooling: PID=%d, temp=%g, entr=%g, dentr=%g\n",
		 ThisTask, P[i].ID, bg_get_temperature(i), SphP[i].Entropy, SphP[i].e.DtEntropy);
	*/

#endif /* BG_COOLING */

	tend = second();
	timecooling += timediff(tstart, tend);


	tstart = second();

#ifdef BG_SFR


	/* ======================================================================
	 *  CHECK WIND AND STAR FORMING PARTICLES
	 * ====================================================================== */

	int oneos;
	double temperature;

	temperature = bg_get_temperature(i);

	/* this flag tells whether a particle fullfills all the star formation criteria:
	   (a) its density is above the threshold for star formation AND
	   (b) (its density is above the over density threshold OR no comoving integration) AND
	   ( (c) ( its temperature is below the EoS threshold AND its density is above the EoS threshold ) OR
	   (d) its temperature is below the fixed temperature threshold )
	   that is  ======= if ( a AND b AND (c OR d) ) ======= */
	oneos = SphP[i].d.Density * a3inv >= All.PhysDensThresh &&
	  ( SphP[i].d.Density > All.OverDensThresh || All.ComovingIntegrationOn == 0 ) &&
	  ( ( log10(SphP[i].Entropy) < log10(bg_get_entropy(i)) + All.SF_EOSTempThreshMargin_DEX &&
	      SphP[i].d.Density * a3inv >= All.EOSPhysDensThresh ) ||
	    temperature < All.SF_THRESH_MaxTemp_K );

#ifdef BG_SNII_KINETIC_FEEDBACK
#ifdef BG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING
	/*
	 * WindFlag = 0     never been in the wind
	 * WindFlag = -aexp when it was first put in the wind
	 *
	 * OnEOS = 0     never been on the EoS
	 * OnEOS = 1     currently on the EoS
	 * OnEOS = -aexp last on the EoS
	 */

	if(SphP[i].WindFlag < 0)
	  {
	    /* elapsed time since the particle was last put in the wind */
	    time_in_wind_in_Gyr = bg_get_elapsed_time(-(double) SphP[i].WindFlag, All.Time, 1); /* note: this is valid for a flat universe! */
	  }
	else
	  time_in_wind_in_Gyr = 0;

	/* if the particle is not in the wind (WindFlag == 0) OR it was in the wind longer than
	   the SF decoupling time, check if it should be put on the effective EoS */
	if(SphP[i].WindFlag == 0 || time_in_wind_in_Gyr > All.SNII_WindDecouplingTime_YR * 1.0e-9)
	  if(SphP[i].OnEOS != 1 && oneos)
	    SphP[i].OnEOS = 1; /* flag the particle and put it on the effective EoS */

        /* check if a particle should leave the effective EoS */
        if(SphP[i].OnEOS == 1 && !oneos)
          SphP[i].OnEOS = -(MyFloat) All.Time; /* flag the particle with the current time */

#else /* BG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING */

	/* check if a particle should be put on the effective EoS */
	if(SphP[i].OnEOS != 1 && oneos)
	  SphP[i].OnEOS = 1; /* flag the particle and put it on the effective EoS */

	/* check if a particle should leave the effective EoS */
	if(SphP[i].OnEOS == 1 && !oneos)
	  SphP[i].OnEOS = -(MyFloat) All.Time; /* flag the particle with the current time */

#endif /* BG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING */

#else /* BG_SNII_KINETIC_FEEDBACK */

	/* check if a particle should be put on the effective EoS */
	if(SphP[i].OnEOS != 1 && oneos)
	  SphP[i].OnEOS = 1; /* flag the particle and put it on the effective EoS */

	/* check if a particle should leave the effective EoS */
	if(SphP[i].OnEOS == 1 && !oneos)
	  SphP[i].OnEOS = -(MyFloat) All.Time; /* flag the particle with the current time */

#endif /* BG_SNII_KINETIC_FEEDBACK */


#ifdef BLACK_HOLES
	if(P[i].Mass == 0)
	  SphP[i].OnEOS = 0;	/* this particle has been swallowed and can't make stars for that reason */
#endif


	/* ======================================================================
	 *  STAR FORMATION
	 * ====================================================================== */


#ifdef TIMESTEP_LIMITER
	if(SphP[i].FlagFeedbackDtLimiter & (1 << BITFLAG_DT_LIMITER))
	  {
	    /* dt from last active time to now */
            dt = (All.Ti_Current - SphP[i].Ti_begstep) * All.Timebase_interval;

	    /* unflag particle */
	    SphP[i].FlagFeedbackDtLimiter ^= (1 << BITFLAG_DT_LIMITER);

	    if(dt < 0)
	      endrun(535353);

            /* the actual time-step */
            if(All.ComovingIntegrationOn)
              dtime = All.Time * dt / time_hubble_a;
            else
              dtime = dt;

	    /*
            printf("[LIMITER, TASK %d]   unflagged particle, ID=%d, flag=%d\n",
                   ThisTask, P[i].ID, SphP[i].FlagFeedbackDtLimiter);
	    */
	  }
#endif


	SphP[i].Sfr = 0;


	if(SphP[i].OnEOS == 1)	/* active star forming particle */
	  {
	    /* the upper bits of the gas particle ID store how many stars this gas
	       particle already generated */
	    if(bits == 0)
	      number_of_stars_generated = 0;
	    else
	      number_of_stars_generated = (P[i].ID >> (32 - bits));

	    mass_of_star = P[i].Mass / (All.Generations - number_of_stars_generated);

	    /* compute the pressure */
	    pressure = SphP[i].Entropy * pow(SphP[i].d.Density, GAMMA);

	    /* convert to physical pressure */
	    pressure *= pow(a3inv, GAMMA);

	    /* convert to cgs units */
	    pressure *= All.UnitPressure_in_cgs;

	    dtime_in_s = dtime * All.UnitTime_in_s / All.HubbleParam;

	    /* mass of stars produced in this time step (code units) */
	    if(SphP[i].d.Density * a3inv >= All.MaxPhysDensThresh && All.SF_THRESH_MaxPhysDensOn == 1)
	      sm = P[i].Mass / (All.Generations - number_of_stars_generated);
	    else
	      {
#ifdef BG_DOUBLE_IMF
		if(SphP[i].d.Density * a3inv > All.IMF_PhysDensThresh)
		  sm = All.SF_SchmidtLawCoeff1_GpSpCM2 * P[i].Mass *
		    pow(pressure * All.HubbleParam * All.HubbleParam * GAMMA / GRAVITY,
			(All.SF_SchmidtLawExponent - 1) / 2) * dtime_in_s;
		else
		  sm = All.SF_SchmidtLawCoeff_GpSpCM2 * P[i].Mass *
		    pow(pressure * All.HubbleParam * All.HubbleParam * GAMMA / GRAVITY,
			(All.SF_SchmidtLawExponent - 1) / 2) * dtime_in_s;
#else
		sm = All.SF_SchmidtLawCoeff_GpSpCM2 * P[i].Mass *
		  pow(pressure * All.HubbleParam * All.HubbleParam * GAMMA / GRAVITY,
		      (All.SF_SchmidtLawExponent - 1) / 2) * dtime_in_s;
#endif
	      }

	    if(sm > P[i].Mass)
	      {
		printf("Mass in new stars greater than particle mass!!! sm = %g\n", sm);
		printf("ID = %lld, entr = %e, dens = %e, pres = %e, dt = %e\n", P[i].ID, SphP[i].Entropy,
		       SphP[i].d.Density, pressure, dtime_in_s);
		sm = P[i].Mass;
	      }

	    sum_sm += sm;
#ifdef BG_POPIII
#ifdef BG_METALSMOOTHING
            if(SphP[i].MetallicitySmoothed < All.POPIII_MetallicityThreshold)
	      popiii_sum_sm += sm;
#else
            if(SphP[i].Metallicity < All.POPIII_MetallicityThreshold)
	      popiii_sum_sm += sm;
#endif
#endif

	    /* store SFR in [solar_mass / year] */
	    if(SphP[i].d.Density * a3inv >= All.MaxPhysDensThresh && All.SF_THRESH_MaxPhysDensOn == 1)
	      SphP[i].Sfr =
		P[i].Mass / (All.Generations - number_of_stars_generated) *
		(All.UnitMass_in_g / SOLAR_MASS) / (dtime_in_s / SEC_PER_YEAR);
	    else
	      {
#ifdef BG_DOUBLE_IMF
		if(SphP[i].d.Density * a3inv > All.IMF_PhysDensThresh)
		  SphP[i].Sfr =
		    All.SF_SchmidtLawCoeff1_GpSpCM2 * P[i].Mass * All.UnitMass_in_g / All.HubbleParam *
		    pow(pressure * All.HubbleParam * All.HubbleParam * GAMMA / GRAVITY,
			(All.SF_SchmidtLawExponent - 1) / 2) * (SEC_PER_YEAR / SOLAR_MASS);
		else
		  SphP[i].Sfr =
		    All.SF_SchmidtLawCoeff_GpSpCM2 * P[i].Mass * All.UnitMass_in_g / All.HubbleParam *
		    pow(pressure * All.HubbleParam * All.HubbleParam * GAMMA / GRAVITY,
			(All.SF_SchmidtLawExponent - 1) / 2) * (SEC_PER_YEAR / SOLAR_MASS);
#else
		SphP[i].Sfr =
		  All.SF_SchmidtLawCoeff_GpSpCM2 * P[i].Mass * All.UnitMass_in_g / All.HubbleParam *
		  pow(pressure * All.HubbleParam * All.HubbleParam * GAMMA / GRAVITY,
		      (All.SF_SchmidtLawExponent - 1) / 2) * (SEC_PER_YEAR / SOLAR_MASS);
#endif
	      }

	    TimeBinSfr[P[i].TimeBin] += SphP[i].Sfr;
#ifdef BG_POPIII
#ifdef BG_METALSMOOTHING
	    if(SphP[i].MetallicitySmoothed < All.POPIII_MetallicityThreshold)
	      TimeBinPOPIIISfr[P[i].TimeBin] += SphP[i].Sfr;
#else
	    if(SphP[i].Metallicity < All.POPIII_MetallicityThreshold)
	      TimeBinPOPIIISfr[P[i].TimeBin] += SphP[i].Sfr;
#endif
#endif

	    if(SphP[i].d.Density * a3inv >= All.MaxPhysDensThresh && All.SF_THRESH_MaxPhysDensOn == 1)
	      {
		printf("[SFR] Task = %d, SphP[%d] above max dens threshold\n", ThisTask, i);
		fflush(stdout);
		prob = 1;
	      }
	    else
	      prob = sm / mass_of_star;

	    if(get_random_number(P[i].ID + 1) < prob)	/* ok, make a star */
	      {
		if(number_of_stars_generated == (All.Generations - 1))
		  {
		    /* here we turn the gas particle itself into a star */
		    Stars_converted++;
		    stars_converted++;

		    sum_mass_stars += P[i].Mass;
#ifdef BG_POPIII
#ifdef BG_METALSMOOTHING
                    if(SphP[i].MetallicitySmoothed < All.POPIII_MetallicityThreshold)
		      popiii_sum_mass_stars += P[i].Mass;
#else
                    if(SphP[i].Metallicity < All.POPIII_MetallicityThreshold)
		      popiii_sum_mass_stars += P[i].Mass;
#endif
#endif
		    P[i].Type = 4;
		    TimeBinCountSph[P[i].TimeBin]--;
		    TimeBinSfr[P[i].TimeBin] -= SphP[i].Sfr;
#ifdef BG_POPIII
#ifdef BG_METALSMOOTHING
		    if(SphP[i].MetallicitySmoothed < All.POPIII_MetallicityThreshold)
		      TimeBinPOPIIISfr[P[i].TimeBin] -= SphP[i].Sfr;
#else
		    if(SphP[i].Metallicity < All.POPIII_MetallicityThreshold)
		      TimeBinPOPIIISfr[P[i].TimeBin] -= SphP[i].Sfr;
#endif
#endif
		    P[i].StarID = N_star;
		    StarP[N_star].PID = i;

		    StarP[N_star].StarBirthTime = All.Time;
		    StarP[N_star].GasDensity = SphP[i].d.Density;
		    StarP[N_star].InitialMass = P[i].Mass;

#ifdef BG_EXTRA_ARRAYS
		    StarP[N_star].MaximumEntropy = SphP[i].MaximumEntropy;
		    StarP[N_star].MaximumTemperature = SphP[i].MaximumTemperature;
		    StarP[N_star].TimeMaximumEntropy = SphP[i].TimeMaximumEntropy;
		    StarP[N_star].TimeMaximumTemperature = SphP[i].TimeMaximumTemperature;
#endif

#ifdef BG_STELLAR_EVOLUTION
		    StarP[N_star].SNIaRate = 0;

#ifdef BG_Z_WEIGHTED_REDSHIFT
		    StarP[N_star].MetallicityWeightedRedshift = SphP[i].MetallicityWeightedRedshift;
#endif

#ifdef BG_SNIA_IRON
		    StarP[N_star].IronFromSNIa = SphP[i].IronFromSNIa;
#endif

#ifdef BG_METALSMOOTHING
		    for(k = 0; k < BG_NELEMENTS; k++)
		      StarP[N_star].MetalsSmoothed[k] = SphP[i].MetalsSmoothed[k];

		    StarP[N_star].MetallicitySmoothed = SphP[i].MetallicitySmoothed;	/* get_particle_smoothed_metallicity(i); */
#ifdef BG_SNIA_IRON
		    StarP[N_star].IronFromSNIaSmoothed = SphP[i].IronFromSNIaSmoothed;
#endif
#endif
		    for(k = 0; k < BG_NELEMENTS; k++)
		      StarP[N_star].Metals[k] = SphP[i].Metals[k];

		    StarP[N_star].Metallicity = SphP[i].Metallicity;
#endif
		    N_star++;
		  }
		else
		  {
		    /* here we spawn a new star particle */

		    if(NumPart + stars_spawned >= All.MaxPart)
		      {
			printf
			  ("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPart=%d)\n",
			   ThisTask, NumPart, stars_spawned, All.MaxPart);
			fflush(stdout);
			endrun(8888);
		      }

		    P[NumPart + stars_spawned] = P[i];
		    P[NumPart + stars_spawned].Type = 4;

		    NextActiveParticle[NumPart + stars_spawned] = FirstActiveParticle;
		    FirstActiveParticle = NumPart + stars_spawned;
		    NumForceUpdate++;

		    TimeBinCount[P[NumPart + stars_spawned].TimeBin]++;

		    PrevInTimeBin[NumPart + stars_spawned] = i;
		    NextInTimeBin[NumPart + stars_spawned] = NextInTimeBin[i];
		    if(NextInTimeBin[i] >= 0)
		      PrevInTimeBin[NextInTimeBin[i]] = NumPart + stars_spawned;
		    NextInTimeBin[i] = NumPart + stars_spawned;
		    if(LastInTimeBin[P[i].TimeBin] == i)
		      LastInTimeBin[P[i].TimeBin] = NumPart + stars_spawned;

		    P[i].ID += (1 << (32 - bits));

		    /* mass of the star particle */
		    P[NumPart + stars_spawned].Mass = mass_of_star;

		    /* fraction of mass removed from the gas particle */
		    frac = mass_of_star / P[i].Mass;

		    /* new mass of the gas particle */
		    P[i].Mass -= P[NumPart + stars_spawned].Mass;

		    sum_mass_stars += P[NumPart + stars_spawned].Mass;
#ifdef BG_POPIII
#ifdef BG_METALSMOOTHING
		    if(SphP[i].MetallicitySmoothed < All.POPIII_MetallicityThreshold)
		      popiii_sum_mass_stars += P[NumPart + stars_spawned].Mass;
#else
		    if(SphP[i].Metallicity < All.POPIII_MetallicityThreshold)
		      popiii_sum_mass_stars += P[NumPart + stars_spawned].Mass;
#endif
#endif
		    P[NumPart + stars_spawned].StarID = N_star;
		    StarP[N_star].PID = NumPart + stars_spawned;

		    StarP[N_star].StarBirthTime = All.Time;
		    StarP[N_star].GasDensity = SphP[i].d.Density;
		    StarP[N_star].InitialMass = mass_of_star;

#ifdef BG_EXTRA_ARRAYS
		    StarP[N_star].MaximumEntropy = SphP[i].MaximumEntropy;
		    StarP[N_star].MaximumTemperature = SphP[i].MaximumTemperature;
		    StarP[N_star].TimeMaximumEntropy = SphP[i].TimeMaximumEntropy;
		    StarP[N_star].TimeMaximumTemperature = SphP[i].TimeMaximumTemperature;
#endif

#ifdef BG_STELLAR_EVOLUTION
		    StarP[N_star].SNIaRate = 0;

#ifdef BG_SNIA_IRON
		    StarP[N_star].IronFromSNIa = frac * SphP[i].IronFromSNIa;
#endif

#ifdef BG_Z_WEIGHTED_REDSHIFT
		    StarP[N_star].MetallicityWeightedRedshift = SphP[i].MetallicityWeightedRedshift;
#endif

#ifdef BG_METALSMOOTHING
		    for(k = 0; k < BG_NELEMENTS; k++)
		      StarP[N_star].MetalsSmoothed[k] = frac * SphP[i].MetalsSmoothed[k];

		    StarP[N_star].MetallicitySmoothed = SphP[i].MetallicitySmoothed;
#ifdef BG_SNIA_IRON
		    StarP[N_star].IronFromSNIaSmoothed = SphP[i].IronFromSNIaSmoothed;
#endif
#endif

		    for(k = 0; k < BG_NELEMENTS; k++)
		      {
			/* metal mass of the star (fraction of the mass removed from the gas particle) */
			StarP[N_star].Metals[k] = frac * SphP[i].Metals[k];
			/* metal mass left in the gas particle */
			SphP[i].Metals[k] -= StarP[N_star].Metals[k];
		      }

		    StarP[N_star].Metallicity = SphP[i].Metallicity;
#endif

		    force_add_star_to_tree(i, NumPart + stars_spawned);

		    stars_spawned++;
		    N_star++;
		  }
	      }
	  }
#endif /* BG_SFR */

	tend = second();
	timesfr += timediff(tstart, tend);

      }				/* end of main loop over active particles */

    CPU_Step[CPU_COOLING] += timecooling;
    CPU_Step[CPU_SFR] += timesfr;


#ifdef BG_SFR
  MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&stars_converted, &tot_converted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(tot_spawned > 0 || tot_converted > 0)
    {
      if(ThisTask == 0)
	{
	  printf("\nSFR: spawned %d stars, converted %d gas particles into stars\n\n",
		 tot_spawned, tot_converted);
	  fflush(stdout);
	}

      All.TotNumPart += tot_spawned;
      All.TotN_gas -= tot_converted;
      All.TotN_star += tot_spawned + tot_converted;
      NumPart += stars_spawned;

      /* Note: N_gas is only reduced once rearrange_particle_sequence is called */

      /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
    }

  for(bin = 0, sfrrate = 0; bin < TIMEBINS; bin++)
    if(TimeBinCount[bin])
      sfrrate += TimeBinSfr[bin];
#ifdef BG_POPIII
  for(bin = 0, popiii_sfrrate = 0; bin < TIMEBINS; bin++)
    if(TimeBinCount[bin])
      popiii_sfrrate += TimeBinPOPIIISfr[bin];
#endif

  MPI_Allreduce(&sfrrate, &totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Reduce(&sum_sm, &total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum_mass_stars, &total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef BG_POPIII
  MPI_Allreduce(&popiii_sfrrate, &popiii_totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Reduce(&popiii_sum_sm, &popiii_total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&popiii_sum_mass_stars, &popiii_total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
  if(ThisTask == 0)
    {
      if(All.TimeStep > 0)
	rate = total_sm / (All.TimeStep / time_hubble_a);
      else
	rate = 0;

      /* convert to solar masses per yr */
      rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      fprintf(FdSfr, "%e  %e  %e  %e  %e\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear,
	      total_sum_mass_stars);
      fflush(FdSfr);
#ifdef BG_POPIII
      if(All.TimeStep > 0)
        rate = popiii_total_sm / (All.TimeStep / time_hubble_a);
      else
        rate = 0;

      /* convert to solar masses per yr */
      rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      fprintf(FdPOPIIISfr, "%e  %e  %e  %e  %e\n", All.Time, popiii_total_sm, popiii_totsfrrate, rate_in_msunperyear,
              popiii_total_sum_mass_stars);
      fflush(FdPOPIIISfr);
#endif
    }
#endif /* BG_SFR */
}


#ifdef BG_SFR
void set_units_sfr(void)
{
  All.OverDensThresh = All.SF_THRESH_MinOverDens *
    All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  All.PhysDensThresh = All.SF_THRESH_MinPhysDens_HpCM3 *
    PROTONMASS / HYDROGEN_MASSFRAC / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);

  All.EOSOverDensThresh = All.SF_EOSMinOverDens *
    All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  All.EOSPhysDensThresh = All.SF_EOSMinPhysDens_HpCM3 *
    PROTONMASS / HYDROGEN_MASSFRAC / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);

  All.MaxPhysDensThresh = All.SF_THRESH_MaxPhysDens_HpCM3 *
    PROTONMASS / HYDROGEN_MASSFRAC / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);

  All.SF_SchmidtLawCoeff_GpSpCM2 = All.SF_SchmidtLawCoeff_MSUNpYRpKPC2 *
    pow(pow(CM_PER_MPC / 1.e6, 2) / SOLAR_MASS, All.SF_SchmidtLawExponent - 1) / (1.e6 * SEC_PER_YEAR) *
    pow(All.GasFraction, (All.SF_SchmidtLawExponent - 1) / 2);
#ifdef BG_DOUBLE_IMF
  All.SF_SchmidtLawCoeff1_GpSpCM2 = All.SF_SchmidtLawCoeff1_MSUNpYRpKPC2 *
    pow(pow(CM_PER_MPC / 1.e6, 2) / SOLAR_MASS, All.SF_SchmidtLawExponent - 1) / (1.e6 * SEC_PER_YEAR) *
    pow(All.GasFraction, (All.SF_SchmidtLawExponent - 1) / 2);

  All.IMF_PhysDensThresh = All.IMF_PhysDensThresh_HpCM3 *
    PROTONMASS / HYDROGEN_MASSFRAC / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);
#endif
#ifdef BG_POPIII
  All.POPIII_MetallicityThreshold = All.POPIII_MetallicityThreshold_SOLAR * ZSOLAR;
#endif
}
#endif


void rearrange_particle_sequence(void)
{
  int flag = 0, flag_sum;

#ifdef BG_SFR
  int i, j;

  struct particle_data psave;
#endif

#ifdef BLACK_HOLES
  int count_elim, count_gaselim, tot_elim, tot_gaselim;
#endif

#ifdef BG_SFR
  if(Stars_converted)
    {
      N_gas -= Stars_converted;
      Stars_converted = 0;

      for(i = 0; i < N_gas; i++)
	if(P[i].Type != 0)
	  {
	    for(j = N_gas; j < NumPart; j++)
	      if(P[j].Type == 0)
		break;

	    if(j >= NumPart)
	      endrun(181170);

	    if(P[i].Type == 4)
	      StarP[P[i].StarID].PID = j;

	    psave = P[i];
	    P[i] = P[j];
	    SphP[i] = SphP[j];
	    P[j] = psave;
	  }
      flag = 1;
    }
#endif

#ifdef BLACK_HOLES
  count_elim = 0;
  count_gaselim = 0;

  for(i = 0; i < NumPart; i++)
    if(P[i].Mass == 0)
      {
	TimeBinCount[P[i].TimeBin]--;

	if(TimeBinActive[P[i].TimeBin])
	  NumForceUpdate--;

	if(P[i].Type == 0)
	  {
	    TimeBinCountSph[P[i].TimeBin]--;

	    P[i] = P[N_gas - 1];
	    SphP[i] = SphP[N_gas - 1];

	    P[N_gas - 1] = P[NumPart - 1];
#ifdef BG_SFR
	    if(P[N_gas - 1].Type == 4)
	      StarP[P[N_gas - 1].StarID].PID = N_gas - 1;
#endif
	    N_gas--;

	    count_gaselim++;
	  }
	else
	  {
	    P[i] = P[NumPart - 1];
#ifdef BG_SFR
	    if(P[i].Type == 4)
	      StarP[P[i].StarID].PID = i;
#endif
	  }

	NumPart--;
	i--;

	count_elim++;
      }

  MPI_Allreduce(&count_elim, &tot_elim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&count_gaselim, &tot_gaselim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(count_elim)
    flag = 1;

  if(ThisTask == 0)
    {
      printf("Blackholes: Eliminated %d gas particles and merged away %d black holes.\n",
	     tot_gaselim, tot_elim - tot_gaselim);
      fflush(stdout);
    }

  All.TotNumPart -= tot_elim;
  All.TotN_gas -= tot_gaselim;
  All.TotBHs -= tot_elim - tot_gaselim;
#endif

  MPI_Allreduce(&flag, &flag_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(flag_sum)
    reconstruct_timebins();

#ifdef BG_SFR
  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 4)
      {
	if((i != (int) StarP[P[i].StarID].PID) || ((int) P[i].StarID >= N_star))
	  {
	    printf
	      ("************\n  [rearrange_particle_sequence] Warning in rearrange on Task %i at particle %i\n************\n"
	       "P.StarID = %u  StarP.ID = %lld N_star = %i\n", ThisTask, i, P[i].StarID, StarP[P[i].StarID].PID,
	       N_star);
	    fflush(stdout);
	    endrun(1212);
	  }
      }
#endif
}

#ifdef BG_COOLING
double bg_get_temperaturep(struct sph_particle_data *sphp, struct particle_data *pp)
{
  int z_index;

  double rho, egycurrent, redshift = 0, a3inv = 1;

  double He_frac, n_H, XH;

  float d_z;

  static int HydrogenIndex, HeliumIndex, first_call = 0;

  if(first_call == 0)
    {
      HydrogenIndex = element_index("Hydrogen");
      HeliumIndex = element_index("Helium");
      first_call = 1;
    }

  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      redshift = 1 / All.Time - 1;
    }
  else
    {
      a3inv = 1;
      redshift = 0;
    }

  /* get d_z for u to temp conversion */
  get_redshift_index(redshift, &z_index, &d_z);

  XH = sphp->Metals[HydrogenIndex] / pp->Mass;
  He_frac = sphp->Metals[HeliumIndex] / pp->Mass;

  /* convert internal energy to physiscal units */
  egycurrent = sphp->Entropy / GAMMA_MINUS1 * pow(sphp->d.Density * a3inv, GAMMA_MINUS1);
  egycurrent *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  /* convert density to physical units */
  rho = sphp->d.Density * a3inv;
  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

  /* convert Hydrogen mass fraction in Hydrogen number density */
  n_H = rho * XH / PROTONMASS;

  /* return the temperature */
  return convert_u_to_temp(d_z, egycurrent, n_H, He_frac);
}

double bg_get_temperature(int i)
{
  return bg_get_temperaturep(&SphP[i], &P[i]);
}
#endif

double bg_get_entropyp(struct sph_particle_data *sphp)
{
#ifdef BG_SFR
  double a3inv;

  double pressure, energy;

  /* NEW */
  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;
  /* NEW */

  if(sphp->d.Density * a3inv >= All.EOSPhysDensThresh)
    {
      /* physical energy in cgs units */
      energy = All.SF_EOSEnergyAtThreshold_ERG * All.UnitMass_in_g / All.UnitEnergy_in_cgs;

      /* physical pressure */
      pressure = GAMMA_MINUS1 * All.EOSPhysDensThresh * energy *
	pow(sphp->d.Density * a3inv / All.EOSPhysDensThresh, All.SF_EOSGammaEffective);

      /* convert pressure to comoving pressure */
      pressure /= pow(a3inv, GAMMA);

      return pressure / pow(sphp->d.Density, GAMMA);
    }
  else
#endif
    return sphp->Entropy;
}

double bg_get_entropy(int i)
{
  return bg_get_entropyp(&SphP[i]);
}

#endif /* closing of BG_SFR-conditional */
