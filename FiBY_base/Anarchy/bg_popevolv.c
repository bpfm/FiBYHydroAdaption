#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <hdf5.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "allvars.h"
#include "proto.h"
#include "bg_cooling.h"
#include "bg_proto.h"
#include "bg_vars.h"
#include "bg_yields.h"

/* global variables for SN feedback/rate (definition in bg_yields.h) */
double NumberOfSNIa, NumberOfSNII;

#ifdef BG_STELLAR_EVOLUTION

/*#define TEST_YIELDS */

/*
 * ----------------------------------------------------------------------
 * This routine applies changes to metal content and mass of star
 * Note that the metallicity, P[i].Metallicity, does not change
 * ----------------------------------------------------------------------
 */

void set_particle_metal_content(int iPart, double mass_released)
{
  int iel;

  double frac;

  if(P[iPart].Type != 4)
    {
      printf("Wrong particle type in set_particle_metal_content!");
      endrun(777);
    }

  if(mass_released > P[iPart].Mass)
    {
      printf
	("[set_particle_metal_content()] Released star mass > star mass: id=%lld, m_released=%e, p_mass=%e\n",
	 P[iPart].ID, mass_released, P[iPart].Mass);
      printf("InitialMass %e\n", StarP[P[iPart].StarID].InitialMass);
      printf("Metallicity %e\n", StarP[P[iPart].StarID].Metallicity);
      for(iel = 0; iel < BG_NELEMENTS; iel++)
	printf("Mass fraction of %12s=%f\n", ElementNames[iel],
	       StarP[P[iPart].StarID].Metals[iel] / P[iPart].Mass);
#ifdef BG_METALSMOOTHING
      printf("Smoothed metallicity %e\n", StarP[P[iPart].StarID].MetallicitySmoothed);
      for(iel = 0; iel < BG_NELEMENTS; iel++)
	printf("Smoothed mass fraction of %12s=%f\n", ElementNames[iel],
	       StarP[P[iPart].StarID].MetalsSmoothed[iel] / StarP[P[iPart].StarID].InitialMass);
#endif
      printf("StarBirthTime %e\n", StarP[P[iPart].StarID].StarBirthTime);
      printf("GasDensity %e\n", StarP[P[iPart].StarID].GasDensity);
#ifdef BG_SNII_KINETIC_FEEDBACK
      printf("WindFlag %e\n", StarP[P[iPart].StarID].WindFlag);
#endif
      if(mass_released > P[iPart].Mass)
	endrun(666);
    }

  frac = 1.0 - mass_released / P[iPart].Mass;

  /*
    NOTE: only particle metallicities are scaled. If smoothed metallicities
    are used, their initial values are taken into account. See few lines
    below for more useful comments.
  */
  for(iel = 0; iel < BG_NELEMENTS; iel++)
    StarP[P[iPart].StarID].Metals[iel] *= frac;

  P[iPart].Mass *= frac;
}

/*
 * ----------------------------------------------------------------------
 * This routine computes the initial element fractions of stars
 * ----------------------------------------------------------------------
 */

void get_particle_metal_content(int iPart, double *initial_metals)
{
  int i;

  double mass_m1;

  if(P[iPart].Type != 4)
    {
      printf("Wrong particle type in get_particle_metal_content!");
      endrun(779);
    }

#ifdef BG_METALSMOOTHING
  /*
     if smoothing is on, the initial element fractions are computed from the
     initial element masses divided by the initial star mass
   */
  mass_m1 = 1 / StarP[P[iPart].StarID].InitialMass;

  for(i = 0; i < BG_NELEMENTS; i++)
    initial_metals[i] = StarP[P[iPart].StarID].MetalsSmoothed[i] * mass_m1;
#else
  /*
     if smoothing is off, the initial element fractions are computed from the
     current element masses divided by the star mass (element and particle masses
     are changing as the star releases mass)
   */
  mass_m1 = 1 / P[iPart].Mass;

  for(i = 0; i < BG_NELEMENTS; i++)
    initial_metals[i] = StarP[P[iPart].StarID].Metals[i] * mass_m1;
#endif
}


/*
 * ----------------------------------------------------------------------
 * This routine does the stellar evolution, and yields the metal release 
 * ----------------------------------------------------------------------
 */
#ifdef BG_SNIA_IRON
void bg_stellar_evolution(double age_of_star_in_Gyr, double dtime_in_Gyr,
			  int iPart, MyFloat * metals_released, MyFloat * metal_mass_released, MyFloat *dust_mass_released,
			  MyFloat * iron_from_snia, MyFloat * energy_released) /* , double *Number_of_SNIa) */
#else
void bg_stellar_evolution(double age_of_star_in_Gyr, double dtime_in_Gyr,
			  int iPart, MyFloat * metals_released, MyFloat * metal_mass_released, MyFloat *dust_mass_released,
			  MyFloat * energy_released) /* , double *Number_of_SNIa) */
#endif
{
  int iel;

  /*  double frac; */
  double initial_metals[BG_NELEMENTS];	/* original element abundances */

  double Star_Z, egy;

  double mass, initial_mass;

  /* double Number_of_SNII = 0; */

  double mass_released;

#ifdef BG_DOUBLE_IMF
  double phys_dens;
#endif

  static int first_call = 1, HydrogenIndex, HeliumIndex;

  if(first_call == 1)
    {
      HydrogenIndex = element_index("Hydrogen");
      HeliumIndex = element_index("Helium");
      first_call = 0;
    }

  double a3inv;

  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  /* *Number_of_SNIa = 0; */

  /* CURRENT mass and INITIAL mass (internal units) */
  mass = P[iPart].Mass;
  initial_mass = StarP[P[iPart].StarID].InitialMass;

  /* get net metal content Z=1-X-Y, and initial element fractions */
#ifdef BG_METALSMOOTHING
  Star_Z = StarP[P[iPart].StarID].MetallicitySmoothed;
#else
  Star_Z = StarP[P[iPart].StarID].Metallicity;
#endif

  get_particle_metal_content(iPart, initial_metals);

  /* zero amount of elements released */
  for(iel = 0; iel < BG_NELEMENTS; iel++)
    metals_released[iel] = 0;

  *metal_mass_released = 0;
#ifdef BG_DUST
  *dust_mass_released = 0;
#endif

  /* physical density */
#ifdef BG_DOUBLE_IMF
  phys_dens = StarP[P[iPart].StarID].GasDensity * a3inv;
#endif

  /* evolve stars. metals_released is in solar masses, for an IMF normalised
     to 1 solar mass, i.e. it is the metal fraction. */
#ifdef BG_SNIA_IRON
  *iron_from_snia = 0;

#ifdef BG_DOUBLE_IMF
  stellar_evolution(age_of_star_in_Gyr, dtime_in_Gyr, Star_Z, initial_metals, metals_released,
		    metal_mass_released, dust_mass_released, iron_from_snia, phys_dens);
  /* metal_mass_released, iron_from_snia, Number_of_SNIa, &Number_of_SNII, phys_dens); */
#else
  stellar_evolution(age_of_star_in_Gyr, dtime_in_Gyr, Star_Z, initial_metals, metals_released,
		    metal_mass_released, dust_mass_released, iron_from_snia);
  /* metal_mass_released, iron_from_snia, Number_of_SNIa, &Number_of_SNII); */
#endif
#else /* BG_SNIA_IRON */
#ifdef BG_DOUBLE_IMF
  stellar_evolution(age_of_star_in_Gyr, dtime_in_Gyr, Star_Z, initial_metals, metals_released,
		    metal_mass_released, dust_mass_released, phys_dens);
  /* metal_mass_released, Number_of_SNIa, &Number_of_SNII, phys_dens); */
#else
  stellar_evolution(age_of_star_in_Gyr, dtime_in_Gyr, Star_Z, initial_metals, metals_released,
		    metal_mass_released, dust_mass_released);
  /* metal_mass_released, Number_of_SNIa, &Number_of_SNII); */
#endif
#endif /* BG_SNIA_IRON */

  /* convert metals_released to Gadget mass units */
  for(iel = 0; iel < BG_NELEMENTS; iel++)
    metals_released[iel] *= initial_mass;

  /* convert metal_mass_released to Gadget mass units */
  *metal_mass_released *= initial_mass;

#ifdef BG_DUST
  /* convert dust_mass_released to Gadget mass units */
  *dust_mass_released *= initial_mass;
#endif

#ifdef BG_SNIA_IRON
  *iron_from_snia *= initial_mass;
#endif

  /* do same for number of SN of each type */
  /* *Number_of_SNIa *= initial_mass * All.UnitMass_in_g / SOLAR_MASS; */
  /* Number_of_SNII *= initial_mass * All.UnitMass_in_g / SOLAR_MASS; */

/*   /\* DEBUG *\/ */
/*   printf("[TASK %d] NumberOfSNIa = %f, NumberOfSNII = %f\n", */
/* 	 ThisTask, NumberOfSNIa, NumberOfSNII); */
/*   /\* DEBUG *\/ */

  NumberOfSNIa *= initial_mass * All.UnitMass_in_g / SOLAR_MASS;
  NumberOfSNII *= initial_mass * All.UnitMass_in_g / SOLAR_MASS;

  for(iel = 0; iel < BG_NELEMENTS; iel++)	/* H, He and metals */
    if(metals_released[iel] < 0)
      {
	printf("Metal mass released of %s for particle %d is negative! %f\n",
	       ElementNames[iel], iPart, *metal_mass_released);
	*metal_mass_released -= metals_released[iel];
	metals_released[iel] = 0;
      }

  if(*metal_mass_released < 0)
    {
      printf("Total metal mass released is negative! %f\n", *metal_mass_released);
      *metal_mass_released = 0;
    }

#ifdef BG_DUST
  if(*dust_mass_released < 0)
    {
      printf("Total dust mass released is negative! %f\n", *dust_mass_released);
      *dust_mass_released = 0;
    }
#endif

  /* compute total mass relased */
  mass_released = *metal_mass_released + metals_released[HydrogenIndex] + metals_released[HeliumIndex];

  /* apply change in mass to star, without changing its metallicity */
  set_particle_metal_content(iPart, mass_released);

  /* compute physical energy released in code units */
  if(All.SNIa_EnergyTransferOn == 1)
    {
      /* egy = *Number_of_SNIa * All.SNIa_Energy_ERG; */
      egy = NumberOfSNIa * All.SNIa_Energy_ERG;
      egy /= All.UnitEnergy_in_cgs;
      *energy_released = egy;
    }
  else
    *energy_released = 0;
}


/*
 * ----------------------------------------------------------------------
 * This routine sets the IMF, normalized such that total mass == 1
 * ----------------------------------------------------------------------
 */

void init_imf(void)
{
  int i;

  double lm_min, lm_max, dlm, lmass, norm, mass;

  imf_by_number = (double *) mymalloc(N_MASS_BIN * sizeof(double));
  imf_mass_bin = (double *) mymalloc(N_MASS_BIN * sizeof(double));
  imf_mass_bin_log10 = (double *) mymalloc(N_MASS_BIN * sizeof(double));
#ifdef BG_DOUBLE_IMF
  imf_by_number1 = (double *) mymalloc(N_MASS_BIN * sizeof(double));
#endif
  stellar_yield = (double *) mymalloc(N_MASS_BIN * sizeof(double));
  integrand = (double *) mymalloc(N_MASS_BIN * sizeof(double));

  if(ThisTask == 0)
    printf("Setting lifetime model to %s\n", All.IMF_LifetimeModel);

  /* set lifetime model */
  if(strcmp(All.IMF_LifetimeModel, "PM93") == 0)
    All.IMF_LifetimeMode = 0;
  else if(strcmp(All.IMF_LifetimeModel, "MM89") == 0)
    All.IMF_LifetimeMode = 1;
  else if(strcmp(All.IMF_LifetimeModel, "P98") == 0)
    All.IMF_LifetimeMode = 2;
  else
    {
      printf("[init]: wrong lifetime model... check your parameter file!\n");
      printf("        valid models are: PM93, MM89, P98\n");
      endrun(-1);
    }

  lm_min = log10(All.IMF_MinMass_MSUN);
  lm_max = log10(All.IMF_MaxMass_MSUN);

  dlm = (lm_max - lm_min) / (double) (N_MASS_BIN - 1);
  norm = 0;

  /* Power-law IMF (Salpeter for IMF_EXPONENT = 2.35) */
  if(strcmp(All.IMF_Model, "PowerLaw") == 0)
    {
      if(All.IMF_Exponent < 0)
	{
	  printf("imf_exponent is supposed to be > 0\n");
	  endrun(-1);
	}

      for(i = 0; i < N_MASS_BIN; i++)
	{
	  lmass = lm_min + i * dlm;
	  mass = pow(10, lmass);

	  imf_by_number[i] = 1.0 / pow(mass, All.IMF_Exponent);

#ifdef BG_DOUBLE_IMF
	  imf_by_number1[i] = 1.0 / pow(mass, All.IMF_Exponent1);
#endif

	  imf_mass_bin[i] = mass;
	  imf_mass_bin_log10[i] = lmass;
	}
    }
  /* Chabrier 2003 */
  else if(strcmp(All.IMF_Model, "Chabrier") == 0)
    {
      for(i = 0; i < N_MASS_BIN; i++)
	{
	  lmass = lm_min + i * dlm;
	  mass = pow(10, lmass);

	  if(mass > 1.0)
	    /* IMF bug fix
	    imf_by_number[i] = 0.241169 * pow(mass, -2.3);
	    */
	    imf_by_number[i] = 0.237912 * pow(mass, -2.3);
	  else
	    /* IMF bug fix
	    imf_by_number[i] =
	      0.864135 * exp(pow((log10(mass) - log10(0.079)), 2.0) / (-2.0 * pow(0.69, 2))) / mass;
	    */
	    imf_by_number[i] =
	      0.852464 * exp(pow((log10(mass) - log10(0.079)), 2.0) / (-2.0 * pow(0.69, 2))) / mass;

#ifdef BG_DOUBLE_IMF
	  imf_by_number1[i] = 1.0 / pow(mass, All.IMF_Exponent1);
#endif

	  imf_mass_bin[i] = mass;
	  imf_mass_bin_log10[i] = lmass;
	}
    }
  else
    {
      printf("[init_imf]: wrong IMF model... check your parameter file!\n");
      printf("            valid models are: PowerLaw and Chabrier\n");
      endrun(-1);
    }

  norm = integrate_imf(log10(All.IMF_MinMass_MSUN), log10(All.IMF_MaxMass_MSUN), 0.0, 1);

  for(i = 0; i < N_MASS_BIN; i++)
    imf_by_number[i] /= norm;

#ifdef BG_DOUBLE_IMF
  norm = integrate_imf(log10(All.IMF_MinMass_MSUN), log10(All.IMF_MaxMass_MSUN), 0.0, 5);

  for(i = 0; i < N_MASS_BIN; i++)
    imf_by_number1[i] /= norm;
#endif
}

/*
 * ----------------------------------------------------------------------
 * This routine integrates the tabulated IMF by number (mode=0)
 * or by mass (mode=1), or integrate yields (mode=2), 
 * ----------------------------------------------------------------------
 */

double integrate_imf(double log_min_dying_mass, double log_max_dying_mass, double m2, int mode)
{
  double result, u, f;

  int ilow, ihigh, index;

  static double dlm, dm;

  static int first_call = 0;

  if(first_call == 0)
    {
      first_call = 1;

      dlm = imf_mass_bin_log10[1] - imf_mass_bin_log10[0];	/* dlog(m) */
    }

  determine_imf_bins(log_min_dying_mass, log_max_dying_mass, &ilow, &ihigh);

  if(mode == 0)
    for(index = ilow; index < ihigh + 1; index++)
      integrand[index] = imf_by_number[index] * imf_mass_bin[index];	/* integrate by number */
  else if(mode == 1)
    for(index = ilow; index < ihigh + 1; index++)
      integrand[index] = imf_by_number[index] * imf_mass_bin[index] * imf_mass_bin[index];	/* integrate by mass */
  else if(mode == 2)
    for(index = ilow; index < ihigh + 1; index++)
      integrand[index] = stellar_yield[index] * imf_by_number[index] * imf_mass_bin[index];	/* integrate number * yield weighted */
  else if(mode == 3)
    for(index = ilow; index < ihigh + 1; index++)
      {
	u = m2 / imf_mass_bin[index];
	f = pow(2.0, gamma_SNIa + 1) * (gamma_SNIa + 1) * pow(u, gamma_SNIa);
	integrand[index] = f * imf_by_number[index] / imf_mass_bin[index];	/* integrate number * f(u) / M  ... type Ia SN */
      }
#ifdef BG_DOUBLE_IMF
  else if(mode == 4)
    for(index = ilow; index < ihigh + 1; index++)
      integrand[index] = imf_by_number1[index] * imf_mass_bin[index];	/* integrate by number */
  else if(mode == 5)
    for(index = ilow; index < ihigh + 1; index++)
      integrand[index] = imf_by_number1[index] * imf_mass_bin[index] * imf_mass_bin[index];	/* integrate by mass */
  else if(mode == 6)
    for(index = ilow; index < ihigh + 1; index++)
      integrand[index] = stellar_yield[index] * imf_by_number1[index] * imf_mass_bin[index];	/* integrate number * yield weighted */
  else if(mode == 7)
    for(index = ilow; index < ihigh + 1; index++)
      {
	u = m2 / imf_mass_bin[index];
	f = pow(2.0, gamma_SNIa + 1) * (gamma_SNIa + 1) * pow(u, gamma_SNIa);
	integrand[index] = f * imf_by_number1[index] / imf_mass_bin[index];	/* integrate number * f(u) / M  ... type Ia SN */
      }
#endif
  else
    {
      printf("invalid mode in integrate_imf = %d\n", mode);
      endrun(-1);
    }

  /* integrate using trapezoidal rule */
  result = 0;
  for(index = ilow; index < ihigh + 1; index++)
    result += integrand[index];

  result = result - 0.5 * integrand[ilow] - 0.5 * integrand[ihigh];

  /* correct first bin */
  dm = (log_min_dying_mass - imf_mass_bin_log10[ilow]) / dlm;

  if(dm < 0.5)
    result -= dm * integrand[ilow];
  else
    {
      result -= 0.5 * integrand[ilow];
      result -= (dm - 0.5) * integrand[ilow + 1];
    }

  /* correct last bin */
  dm = (log_max_dying_mass - imf_mass_bin_log10[ihigh - 1]) / dlm;

  if(dm < 0.5)
    {
      result -= 0.5 * integrand[ihigh];
      result -= (0.5 - dm) * integrand[ihigh - 1];
    }
  else
    result -= (1 - dm) * integrand[ihigh];

  result *= dlm * log(10.0);	/* log(10) since mass function tabulated as function of log_10(mass) */

  return result;
}


void determine_imf_bins(double log_min_dying_mass, double log_max_dying_mass, int *ilow, int *ihigh)
{
  int i1, i2;

  if(log_min_dying_mass < imf_mass_bin_log10[0])
    log_min_dying_mass = imf_mass_bin_log10[0];

  if(log_min_dying_mass > imf_mass_bin_log10[N_MASS_BIN - 1])
    log_min_dying_mass = imf_mass_bin_log10[N_MASS_BIN - 1];

  if(log_max_dying_mass < imf_mass_bin_log10[0])
    log_max_dying_mass = imf_mass_bin_log10[0];

  if(log_max_dying_mass > imf_mass_bin_log10[N_MASS_BIN - 1])
    log_max_dying_mass = imf_mass_bin_log10[N_MASS_BIN - 1];

  for(i1 = 0; i1 < N_MASS_BIN - 2 && imf_mass_bin_log10[i1 + 1] < log_min_dying_mass; i1++);

  for(i2 = 1; i2 < N_MASS_BIN - 1 && imf_mass_bin_log10[i2] < log_max_dying_mass; i2++);

  *ilow = i1;
  *ihigh = i2;
}


/*
 * ----------------------------------------------------------------------
 * This routine determines mass released per element, and energy feedback,
 * from TypeIa and TypeII SNe, and AGB stars, over the time-interval
 *
 * [age_of_star_in_Gyr, dtime_in_Gyr],
 *
 * For a stellar population with given metallicity metallicity
 * (defined as M_metal/M_total).
 *
 * Results are for a given IMF, normalsed to total mass = 1 solar mass.
 * ----------------------------------------------------------------------
 */

#ifdef BG_SNIA_IRON
#ifdef BG_DOUBLE_IMF
void stellar_evolution(double age_of_star_in_Gyr, double dtime_in_Gyr, double metallicity,
		       double *initial_metals, MyFloat * metals_released, MyFloat * metal_mass_released,
		       MyFloat * dust_mass_released, MyFloat * iron_from_snia, double phys_dens)
     /* MyFloat * iron_from_snia, double *Number_of_SNIa, double *Number_of_SNII, double phys_dens) */
#else
void stellar_evolution(double age_of_star_in_Gyr, double dtime_in_Gyr, double metallicity,
		       double *initial_metals, MyFloat * metals_released, MyFloat * metal_mass_released,
		       MyFloat * dust_mass_released, MyFloat * iron_from_snia)
     /* MyFloat * iron_from_snia, double *Number_of_SNIa, double *Number_of_SNII) */
#endif
#else /* BG_SNIA_IRON */
#ifdef BG_DOUBLE_IMF
void stellar_evolution(double age_of_star_in_Gyr, double dtime_in_Gyr, double metallicity,
		       double *initial_metals, MyFloat * metals_released, MyFloat * metal_mass_released,
		       MyFloat * dust_mass_released, double phys_dens)
     /* double *Number_of_SNIa, double *Number_of_SNII, double phys_dens) */
#else
void stellar_evolution(double age_of_star_in_Gyr, double dtime_in_Gyr, double metallicity,
		       double *initial_metals, MyFloat * metals_released, MyFloat * metal_mass_released,
		       MyFloat * dust_mass_released)
     /* double *Number_of_SNIa, double *Number_of_SNII) */
#endif
#endif				/* BG_SNIA_IRON */
{
  double log_min_dying_mass, log_max_dying_mass, tmp, log_metallicity;

  double Ntot, Mtot;

  static int first_call = -1;

  /* test normalisation */
  if(first_call == -1)
    {
      log_min_dying_mass = log10(All.IMF_MinMass_MSUN);
      log_max_dying_mass = log10(All.IMF_MaxMass_MSUN);

      Ntot = integrate_imf(log_min_dying_mass, log_max_dying_mass, 0.0, 0);
      Mtot = integrate_imf(log_min_dying_mass, log_max_dying_mass, 0.0, 1);

      if(ThisTask == 0)
	printf("[stellar_evolution] Test of normalisation of mass function. Ntot = %g\t Mtot = %g\n", Ntot,
	       Mtot);

      first_call = 0;
    }

  /* minimum and maximum mass of stars that will die during this time-step */
#if BG_DEBUG >= 1
  if(dying_mass_msun(age_of_star_in_Gyr, metallicity) <= 0
     || dying_mass_msun(age_of_star_in_Gyr + dtime_in_Gyr, metallicity) <= 0)
    {
      printf("[stellar_evolution] error: dying mass <= 0 %f, %f\n",
	     dying_mass_msun(age_of_star_in_Gyr, metallicity),
	     dying_mass_msun(age_of_star_in_Gyr + dtime_in_Gyr, metallicity));
      endrun(801);
    }
#endif
  log_max_dying_mass = log10(dying_mass_msun(age_of_star_in_Gyr, metallicity));
  log_min_dying_mass = log10(dying_mass_msun(age_of_star_in_Gyr + dtime_in_Gyr, metallicity));

  if(log_min_dying_mass > log_max_dying_mass)
    {
      tmp = log_max_dying_mass;
      log_max_dying_mass = log_min_dying_mass;
      log_min_dying_mass = tmp;

      printf("[stellar_evolution] well, this should not happen, really %g, %g, %g, %g\n",
	     log_min_dying_mass, log_max_dying_mass, age_of_star_in_Gyr, age_of_star_in_Gyr + dtime_in_Gyr);

      endrun(911);
    }

  /* integration interval is zero - this can happen if minimum and maximum
   * dying masses are above All.IMF_MaxMass_MSUN */
  if(log_min_dying_mass == log_max_dying_mass)
    return;

  if(metallicity > 0)
    log_metallicity = log10(metallicity);
  else
    log_metallicity = MIN_METAL;

#ifdef BG_POPIII
  if(metallicity >= All.POPIII_MetallicityThreshold)
    {
#endif
      /* contribution from TypeIa SNe */
#ifdef BG_SNIA_IRON
#ifdef BG_DOUBLE_IMF
      evolve_SNIa(log_min_dying_mass, log_max_dying_mass, dtime_in_Gyr, age_of_star_in_Gyr,
		  metallicity, metals_released, metal_mass_released, iron_from_snia, 2, phys_dens);
      /* metallicity, metals_released, metal_mass_released, iron_from_snia, Number_of_SNIa, 2, phys_dens); */
#else
      evolve_SNIa(log_min_dying_mass, log_max_dying_mass, dtime_in_Gyr, age_of_star_in_Gyr,
		  metallicity, metals_released, metal_mass_released, iron_from_snia, 2);
      /* metallicity, metals_released, metal_mass_released, iron_from_snia, Number_of_SNIa, 2); */
#endif
#else /* BG_SNIA_IRON */
#ifdef BG_DOUBLE_IMF
      evolve_SNIa(log_min_dying_mass, log_max_dying_mass, dtime_in_Gyr, age_of_star_in_Gyr,
		  metallicity, metals_released, metal_mass_released, 2, phys_dens);
      /* metallicity, metals_released, metal_mass_released, Number_of_SNIa, 2, phys_dens); */
#else
      evolve_SNIa(log_min_dying_mass, log_max_dying_mass, dtime_in_Gyr, age_of_star_in_Gyr,
		  metallicity, metals_released, metal_mass_released, 2);
      /* metallicity, metals_released, metal_mass_released, Number_of_SNIa, 2); */
#endif
#endif /* BG_SNIA_IRON */

  /* contribution from TypeII SNe */
#ifdef BG_DOUBLE_IMF
      evolve_SNII(log_min_dying_mass, log_max_dying_mass, log_metallicity,
		  initial_metals, metals_released, metal_mass_released, dust_mass_released, phys_dens);
      /* initial_metals, metals_released, metal_mass_released, Number_of_SNII, phys_dens); */
#else
      evolve_SNII(log_min_dying_mass, log_max_dying_mass, log_metallicity,
		  initial_metals, metals_released, metal_mass_released, dust_mass_released);
      /* initial_metals, metals_released, metal_mass_released, Number_of_SNII); */
#endif

  /* contribution from AGB stars */
#ifdef BG_DOUBLE_IMF
      evolve_AGB(log_min_dying_mass, log_max_dying_mass, log_metallicity, initial_metals, metals_released,
		 metal_mass_released, dust_mass_released, phys_dens);
#else
      evolve_AGB(log_min_dying_mass, log_max_dying_mass, log_metallicity, initial_metals, metals_released,
		 metal_mass_released, dust_mass_released);
#endif

  /* contribution from massive popIII stars */
#ifdef BG_POPIII
    }
  else
    {
      evolve_POPIII(log_min_dying_mass, log_max_dying_mass, log_metallicity, initial_metals, metals_released,
		    metal_mass_released, dust_mass_released);
    }
#endif
}

/*
 * ----------------------------------------------------------------------
 * This routine computes yields (per unit solar mass) of AGB stars
 * ----------------------------------------------------------------------
 */

#ifdef BG_DOUBLE_IMF
void evolve_AGB(double log_min_mass, double log_max_mass, double log_metallicity,
		double *initial_metals, MyFloat * metals_released,
		MyFloat * metal_mass_released, MyFloat * dust_mass_released,
		double phys_dens)
#else
void evolve_AGB(double log_min_mass, double log_max_mass, double log_metallicity,
		double *initial_metals, MyFloat * metals_released,
		MyFloat * metal_mass_released, MyFloat * dust_mass_released)
#endif
{
  int ilow, ihigh, imass;

  int iz_low, iz_high, i;

  double dz, deltaz, metallicity;

#ifdef BG_DUST
  MyFloat dust_mass;
#endif

  MyFloat metals[BG_NELEMENTS], mass;

  if(All.AGB_MassTransferOn == 0)
    return;

  metallicity = pow(10.0, log_metallicity);

  /* determine integration range, limiting to stars that become AGB stars */
  if(log_max_mass > log10(All.SNII_MinMass_MSUN))
    log_max_mass = log10(All.SNII_MinMass_MSUN);

  if(log_min_mass >= log_max_mass)
    return;

  /* determine which mass bins will contribute */
  determine_imf_bins(log_min_mass, log_max_mass, &ilow, &ihigh);

  /* determine yield of these bins (not equally spaced bins) */
  if(log_metallicity > MIN_METAL)
    {
      for(iz_low = 0; iz_low < yieldsAGB.N_Z - 1 &&
	  log_metallicity > yieldsAGB.Metallicity[iz_low + 1]; iz_low++);

      iz_high = iz_low + 1;

      if(iz_high >= yieldsAGB.N_Z)
	iz_high = yieldsAGB.N_Z - 1;

      if(log_metallicity >= yieldsAGB.Metallicity[0] &&
	 log_metallicity <= yieldsAGB.Metallicity[yieldsAGB.N_Z - 1])
	dz = log_metallicity - yieldsAGB.Metallicity[iz_low];
      else
	dz = 0;

      deltaz = yieldsAGB.Metallicity[iz_high] - yieldsAGB.Metallicity[iz_low];

      if(deltaz > 0)
	dz = dz / deltaz;
      else
	dz = 0;
    }
  else
    {
      iz_low = 0;
      iz_high = 0;
      dz = 0;
    }

  /* compute stellar_yield as function of mass */
  for(i = 0; i < BG_NELEMENTS; i++)
    {
      for(imass = ilow; imass < ihigh + 1; imass++)
	/* yieldsAGB.SPH refers to elements produced, yieldsAGB.Ejecta_SPH refers to elements already in star */
	stellar_yield[imass] =
	  (1 - dz) * (yieldsAGB.SPH[iz_low][i][imass] +
		      initial_metals[i] * yieldsAGB.Ejecta_SPH[iz_low][imass]) +
	  dz * (yieldsAGB.SPH[iz_high][i][imass] + initial_metals[i] * yieldsAGB.Ejecta_SPH[iz_high][imass]);

#ifdef BG_DOUBLE_IMF
      if(phys_dens > All.IMF_PhysDensThresh)
	metals[i] = integrate_imf(log_min_mass, log_max_mass, 0.0, 6);
      else
	metals[i] = integrate_imf(log_min_mass, log_max_mass, 0.0, 2);
#else
      metals[i] = integrate_imf(log_min_mass, log_max_mass, 0.0, 2);
#endif
    }

  for(imass = ilow; imass < ihigh + 1; imass++)
    stellar_yield[imass] = (1 - dz) * (yieldsAGB.TotalMetals_SPH[iz_low][imass] +
				       metallicity * yieldsAGB.Ejecta_SPH[iz_low][imass]) +
      dz * (yieldsAGB.TotalMetals_SPH[iz_high][imass] +
	    metallicity * yieldsAGB.Ejecta_SPH[iz_high][imass]);

#ifdef BG_DOUBLE_IMF
  if(phys_dens > All.IMF_PhysDensThresh)
    mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 6);
  else
    mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 2);
#else
  mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 2);
#endif

#ifdef BG_DUST
  /* determine yield of these bins (not equally spaced bins) */
  if(log_metallicity > MIN_METAL)
    {
      for(iz_low = 0; iz_low < dustAGB.N_Z - 1 &&
          log_metallicity > dustAGB.Metallicity[iz_low + 1]; iz_low++);

      iz_high = iz_low + 1;

      if(iz_high >= dustAGB.N_Z)
	iz_high = dustAGB.N_Z - 1;

      if(log_metallicity >= dustAGB.Metallicity[0] &&
         log_metallicity <= dustAGB.Metallicity[dustAGB.N_Z - 1])
	dz = log_metallicity - dustAGB.Metallicity[iz_low];
      else
        dz = 0;

      deltaz = dustAGB.Metallicity[iz_high] - dustAGB.Metallicity[iz_low];

      if(deltaz > 0)
        dz = dz / deltaz;
      else
        dz = 0;
    }
  else
    {
      iz_low = 0;
      iz_high = 0;
      dz = 0;
    }

  /* this may be wrong */
  for(imass = ilow; imass < ihigh + 1; imass++)
    stellar_yield[imass] = (1 - dz) * dustAGB.TotalDust_SPH[iz_low][imass] +
      dz * dustAGB.TotalDust_SPH[iz_high][imass];

#ifdef BG_DOUBLE_IMF
  if(phys_dens > All.IMF_PhysDensThresh)
    dust_mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 6);
  else
    dust_mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 2);
#else
  dust_mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 2);
#endif
#endif

/* yield normalization */
  int HydrogenIndex, HeliumIndex;

  double norm0, norm1;

  HydrogenIndex = element_index("Hydrogen");
  HeliumIndex = element_index("Helium");

  /* zero all negative values */
  for(i = 0; i < BG_NELEMENTS; i++)
    if(metals[i] < 0)
      metals[i] = 0;

  if(mass < 0)
    mass = 0;

  /* get the total mass ejected from the table */
  for(imass = ilow; imass < ihigh + 1; imass++)
    stellar_yield[imass] = (1 - dz) * yieldsAGB.Ejecta_SPH[iz_low][imass] +
      dz * yieldsAGB.Ejecta_SPH[iz_high][imass];

  norm0 = integrate_imf(log_min_mass, log_max_mass, 0.0, 2);

  /* compute the total mass ejected */
  norm1 = mass + metals[HydrogenIndex] + metals[HeliumIndex];
#ifdef BG_DUST
  norm1 += dust_mass;
#endif

  /* normalize the yields */
  if(norm1 > 0)
    {
      for(i = 0; i < BG_NELEMENTS; i++)
	metals_released[i] += metals[i] * (norm0 / norm1);

      *metal_mass_released += mass * (norm0 / norm1);
#ifdef BG_DUST
      *dust_mass_released += dust_mass * (norm0 / norm1);
#endif
    }
  else
    {
      printf("[evolve_AGB] wrong normalization!!!! norm1 = %e\n", norm1);
      endrun(666);
    }
}

/*
 * ----------------------------------------------------------------------
 * This routine computes yields from and number of SNIa stars
 * ----------------------------------------------------------------------
 */
#ifdef BG_SNIA_IRON
#ifdef BG_DOUBLE_IMF
void evolve_SNIa(double log_min_dying_mass, double log_max_dying_mass,
		 double dt_in_Gyr, double age_of_star_in_Gyr, double metallicity,
		 MyFloat * metals_released, MyFloat * metal_mass_released,
		 MyFloat * iron_from_snia, int mode, double phys_dens)
     /* MyFloat * iron_from_snia, double *Number_of_SNIa, int mode, double phys_dens) */
#else
void evolve_SNIa(double log_min_dying_mass, double log_max_dying_mass,
		 double dt_in_Gyr, double age_of_star_in_Gyr, double metallicity,
		 MyFloat * metals_released, MyFloat * metal_mass_released,
		 MyFloat * iron_from_snia, int mode)
     /* MyFloat * iron_from_snia, double *Number_of_SNIa, int mode) */
#endif
#else
#ifdef BG_DOUBLE_IMF
void evolve_SNIa(double log_min_dying_mass, double log_max_dying_mass,
		 double dt_in_Gyr, double age_of_star_in_Gyr, double metallicity,
		 MyFloat * metals_released, MyFloat * metal_mass_released,
		 int mode, double phys_dens)
     /* double *Number_of_SNIa, int mode, double phys_dens) */
#else
void evolve_SNIa(double log_min_dying_mass, double log_max_dying_mass,
		 double dt_in_Gyr, double age_of_star_in_Gyr, double metallicity,
		 MyFloat * metals_released, MyFloat * metal_mass_released, int mode)
     /* MyFloat * metals_released, MyFloat * metal_mass_released, double *Number_of_SNIa, int mode) */
#endif
#endif
{
  double dm_dtau, part1, part2, m_bmin, m_bmax, frac_wd_beg, frac_wd_end;

  double d_tau = 0.000001, num_of_SNIa_per_msun = 0.0;

  int i;

  double temp_in_Gyr;

  /* outside the mass range for SNIa */
  if(log_min_dying_mass >= log10(TYPEIA_MAX_MASS))
    return;

  /* check integration limits */
  if(log_max_dying_mass > log10(TYPEIA_MAX_MASS))
    {
      log_max_dying_mass = log10(TYPEIA_MAX_MASS);
      temp_in_Gyr = lifetime_in_Gyr(TYPEIA_MAX_MASS, metallicity);
      dt_in_Gyr = age_of_star_in_Gyr + dt_in_Gyr - temp_in_Gyr;
      age_of_star_in_Gyr = temp_in_Gyr;
    }

  /* compute the fraction of white dwarfs */
  /* (note to rob: since the IMF normalises to 1, and wd number cancels, number_wd is not needed) */
#ifdef BG_DOUBLE_IMF
  if(phys_dens > All.IMF_PhysDensThresh)
    {
      frac_wd_beg = integrate_imf(log_max_dying_mass, log10(TYPEIA_MAX_MASS), 0.0, 4);
      frac_wd_end = integrate_imf(log_min_dying_mass, log10(TYPEIA_MAX_MASS), 0.0, 4);
    }
  else
    {
      frac_wd_beg = integrate_imf(log_max_dying_mass, log10(TYPEIA_MAX_MASS), 0.0, 0);
      frac_wd_end = integrate_imf(log_min_dying_mass, log10(TYPEIA_MAX_MASS), 0.0, 0);
    }
#else
  frac_wd_beg = integrate_imf(log_max_dying_mass, log10(TYPEIA_MAX_MASS), 0.0, 0);
  frac_wd_end = integrate_imf(log_min_dying_mass, log10(TYPEIA_MAX_MASS), 0.0, 0);
#endif

  switch (All.SNIa_Mode)
    {
      /* Greggio and Renzini 1983 (as in Lia et al. 2002) */
    case 0:
      {
	/* part 1 (at t = age_of_star_in_Gyr) */
	dm_dtau = (dying_mass_msun(age_of_star_in_Gyr, metallicity) -
		   dying_mass_msun(age_of_star_in_Gyr + d_tau, metallicity)) / d_tau;

	m_bmin = max(3.0, 2.0 * dying_mass_msun(age_of_star_in_Gyr, metallicity));
	m_bmax = min(2.0 * TYPEIA_MAX_MASS,
		     dying_mass_msun(age_of_star_in_Gyr, metallicity) + TYPEIA_MAX_MASS);

#ifdef BG_DOUBLE_IMF
	if(phys_dens > All.IMF_PhysDensThresh)
	  part1 = dm_dtau * integrate_imf(log10(m_bmin), log10(m_bmax),
					  dying_mass_msun(age_of_star_in_Gyr, metallicity), 7);
	else
	  part1 = dm_dtau * integrate_imf(log10(m_bmin), log10(m_bmax),
					  dying_mass_msun(age_of_star_in_Gyr, metallicity), 3);
#else
	part1 = dm_dtau * integrate_imf(log10(m_bmin), log10(m_bmax),
					dying_mass_msun(age_of_star_in_Gyr, metallicity), 3);
#endif

	/* part 2 (at t = age_of_star_in_Gyr + dt_in_Gyr) */
	dm_dtau = (dying_mass_msun(age_of_star_in_Gyr + dt_in_Gyr, metallicity) -
		   dying_mass_msun(age_of_star_in_Gyr + dt_in_Gyr + d_tau, metallicity)) / d_tau;

	m_bmin = max(3.0, 2.0 * dying_mass_msun(age_of_star_in_Gyr + dt_in_Gyr, metallicity));
	m_bmax = min(2.0 * TYPEIA_MAX_MASS,
		     dying_mass_msun(age_of_star_in_Gyr + dt_in_Gyr, metallicity) + TYPEIA_MAX_MASS);

#ifdef BG_DOUBLE_IMF
	if(phys_dens > All.IMF_PhysDensThresh)
	  part2 = dm_dtau * integrate_imf(log10(m_bmin), log10(m_bmax),
					  dying_mass_msun(age_of_star_in_Gyr + dt_in_Gyr, metallicity), 7);
	else
	  part2 = dm_dtau * integrate_imf(log10(m_bmin), log10(m_bmax),
					  dying_mass_msun(age_of_star_in_Gyr + dt_in_Gyr, metallicity), 3);
#else
	part2 = dm_dtau * integrate_imf(log10(m_bmin), log10(m_bmax),
					dying_mass_msun(age_of_star_in_Gyr + dt_in_Gyr, metallicity), 3);
#endif

	/* num_of_SNIa_per_msun = total_mass * A_SNIa_Lia*(part1+part1)*dt*1.0e9/(2.0); */
	num_of_SNIa_per_msun = A_SNIa_Lia * (part1 + part2) * dt_in_Gyr / 2.0;
      }

      break;

      /* Gaussian (Forster 2006) */
    case 1:
      num_of_SNIa_per_msun = All.SNIa_Efficiency_fracwd * dt_in_Gyr * (frac_wd_beg * exp(pow(age_of_star_in_Gyr - TAU_SNIa_MAN_G, 2.0) / (-2.0 * pow(SIGMA_SNIa_G, 2.0))) + frac_wd_end * exp(pow(age_of_star_in_Gyr - TAU_SNIa_MAN_G + dt_in_Gyr, 2.0) / (-2.0 * pow(SIGMA_SNIa_G, 2.0)))) / (2.0 * sqrt(2 * PI * pow(SIGMA_SNIa_G, 2.0)));	/* gaussian */

      break;

      /* Efolding (Forster 2006) */
    case 2:
      num_of_SNIa_per_msun = All.SNIa_Efficiency_fracwd * dt_in_Gyr * (frac_wd_beg * exp(-1.0 * age_of_star_in_Gyr / TAU_SNIa_MAN_E) + frac_wd_end * exp(-1.0 * (age_of_star_in_Gyr + dt_in_Gyr) / TAU_SNIa_MAN_E)) / (2.0 * TAU_SNIa_MAN_E);	/* efolding */

      break;

      /* GaussianEfolding (Mannucci 2005) */
    case 3:
      {
	num_of_SNIa_per_msun = 0.5 * All.SNIa_Efficiency_fracwd * dt_in_Gyr * (exp(pow(age_of_star_in_Gyr - TAU_SNIa_MAN_GC, 2.0) / (-2.0 * pow(SIGMA_SNIa_GC, 2.0))) + exp(pow(age_of_star_in_Gyr - TAU_SNIa_MAN_GC + dt_in_Gyr, 2.0) / (-2.0 * pow(SIGMA_SNIa_GC, 2.0)))) / (2.0 * sqrt(2 * PI * pow(SIGMA_SNIa_GC, 2.0)));	/* gaussian */

	num_of_SNIa_per_msun += 0.5 * All.SNIa_Efficiency_fracwd * dt_in_Gyr * (frac_wd_beg * exp(-1.0 * age_of_star_in_Gyr / TAU_SNIa_MAN_EC) + frac_wd_end * exp(-1.0 * (age_of_star_in_Gyr + dt_in_Gyr) / TAU_SNIa_MAN_EC)) / (2.0 * TAU_SNIa_MAN_EC);	/*efold (combo) */
      }

      break;

    default:
      puts("SNIa mode not defined yet");
      /* Scannapieco and Bildston 2005 not implemented. */
      /* This formula depends on the actual (NOT initial) mass and star formation rate. */
      /* The code still doesn't pass these values. To be updated if needed. */
      /* num_of_SNIa_per_msun = (A_SNIa_Scan * (total_mass/1.0E10) + B_SNIa_Scan * (SFR/1.0E10)) * dt * 1.0e9/100.0; */
    }

  /* *Number_of_SNIa = num_of_SNIa_per_msun; */
  NumberOfSNIa = num_of_SNIa_per_msun;

  if(All.SNIa_MassTransferOn == 1)
    {
      for(i = 0; i < BG_NELEMENTS; i++)
	metals_released[i] += num_of_SNIa_per_msun * yieldsSNIa.SPH[i];

      *metal_mass_released += num_of_SNIa_per_msun * yieldsSNIa.TotalMetals_SPH;
    }

#ifdef BG_SNIA_IRON
  if(All.SNIa_MassTransferOn == 1)
    *iron_from_snia = metals_released[element_index("Iron")];
#endif
}

/*
 * ----------------------------------------------------------------------
 * This routine computes yields from and number of SNII stars
 * ----------------------------------------------------------------------
 */

#ifdef BG_DOUBLE_IMF
void evolve_SNII(double log_min_mass, double log_max_mass, double log_metallicity,
		 double *initial_metals, MyFloat * metals_released,
		 MyFloat * metal_mass_released, MyFloat * dust_mass_released,
		 double phys_dens)
     /* double *Number_of_SNII, double phys_dens) */
#else
void evolve_SNII(double log_min_mass, double log_max_mass, double log_metallicity,
		 double *initial_metals, MyFloat * metals_released,
		 MyFloat * metal_mass_released, MyFloat * dust_mass_released)
     /* double *Number_of_SNII) */
#endif
{
  int ilow, ihigh, imass;

  int iz_low, iz_high, i;

  double dz, deltaz, metallicity;

#ifdef BG_DUST
  MyFloat dust_mass;
#endif

  MyFloat metals[BG_NELEMENTS], mass;

  /* *Number_of_SNII = 0; */

  if(All.SNII_MassTransferOn == 0)
    return;

  metallicity = pow(10.0, log_metallicity);

  /* determine integration range: make sure all these stars actually become SN of type II */
  if(log_min_mass < log10(All.SNII_MinMass_MSUN))
    log_min_mass = log10(All.SNII_MinMass_MSUN);

  if(log_max_mass > log10(All.SNII_MaxMass_MSUN))
    log_max_mass = log10(All.SNII_MaxMass_MSUN);

  if(log_min_mass >= log_max_mass)
    return;

  /* determine which mass bins will contribute */
  determine_imf_bins(log_min_mass, log_max_mass, &ilow, &ihigh);
  /*
#ifdef BG_DOUBLE_IMF
  if(phys_dens > All.IMF_PhysDensThresh)
    *Number_of_SNII = integrate_imf(log_min_mass, log_max_mass, 0.0, 4);
  else
    *Number_of_SNII = integrate_imf(log_min_mass, log_max_mass, 0.0, 0);
#else
  *Number_of_SNII = integrate_imf(log_min_mass, log_max_mass, 0.0, 0);
#endif
  */

#ifdef BG_DOUBLE_IMF
  if(phys_dens > All.IMF_PhysDensThresh)
    NumberOfSNII = integrate_imf(log_min_mass, log_max_mass, 0.0, 4);
  else
    NumberOfSNII = integrate_imf(log_min_mass, log_max_mass, 0.0, 0);
#else
  NumberOfSNII = integrate_imf(log_min_mass, log_max_mass, 0.0, 0);
#endif

  /* determine yield of these bins */
  if(log_metallicity > MIN_METAL)
    {
      for(iz_low = 0;
	  iz_low < yieldsSNII.N_Z - 1 && log_metallicity > yieldsSNII.Metallicity[iz_low + 1]; iz_low++);

      iz_high = iz_low + 1;

      if(iz_high >= yieldsSNII.N_Z)
	iz_high = yieldsSNII.N_Z - 1;

      if(log_metallicity >= yieldsSNII.Metallicity[0]
	 && log_metallicity <= yieldsSNII.Metallicity[yieldsSNII.N_Z - 1])
	dz = log_metallicity - yieldsSNII.Metallicity[iz_low];
      else
	dz = 0;

      deltaz = yieldsSNII.Metallicity[iz_high] - yieldsSNII.Metallicity[iz_low];

      if(deltaz > 0)
	dz = dz / deltaz;
      else
	dz = 0;
    }
  else
    {
      iz_low = 0;
      iz_high = 0;
      dz = 0;
    }

  /* .... compute stellar_yield as function of mass */
  for(i = 0; i < BG_NELEMENTS; i++)
    {
      for(imass = ilow; imass < ihigh + 1; imass++)
	stellar_yield[imass] =
	  (1 - dz) * (yieldsSNII.SPH[iz_low][i][imass] +
		      (initial_metals[i] * yieldsSNII.Ejecta_SPH[iz_low][imass])) +
	  dz * (yieldsSNII.SPH[iz_high][i][imass] +
		(initial_metals[i] * yieldsSNII.Ejecta_SPH[iz_high][imass]));

#ifdef BG_DOUBLE_IMF
      if(phys_dens > All.IMF_PhysDensThresh)
	metals[i] = integrate_imf(log_min_mass, log_max_mass, 0.0, 6);
      else
	metals[i] = integrate_imf(log_min_mass, log_max_mass, 0.0, 2);
#else
      metals[i] = integrate_imf(log_min_mass, log_max_mass, 0.0, 2);
#endif
    }

  for(imass = ilow; imass < ihigh + 1; imass++)
    stellar_yield[imass] = (1 - dz) * (yieldsSNII.TotalMetals_SPH[iz_low][imass] +
				       metallicity * yieldsSNII.Ejecta_SPH[iz_low][imass]) +
      dz * (yieldsSNII.TotalMetals_SPH[iz_high][imass] +
	    metallicity * yieldsSNII.Ejecta_SPH[iz_high][imass]);


#ifdef BG_DOUBLE_IMF
  if(phys_dens > All.IMF_PhysDensThresh)
    mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 6);
  else
    mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 2);
#else
  mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 2);
#endif

#ifdef BG_DUST
  /* determine yield of these bins (not equally spaced bins) */
  if(log_metallicity > MIN_METAL)
    {
      for(iz_low = 0; iz_low < dustSNII.N_Z - 1 &&
	    log_metallicity > dustSNII.Metallicity[iz_low + 1]; iz_low++);

      iz_high = iz_low + 1;

      if(iz_high >= dustSNII.N_Z)
	iz_high = dustSNII.N_Z - 1;

      if(log_metallicity >= dustSNII.Metallicity[0] &&
         log_metallicity <= dustSNII.Metallicity[dustSNII.N_Z - 1])
        dz = log_metallicity - dustSNII.Metallicity[iz_low];
      else
        dz = 0;

      deltaz = dustSNII.Metallicity[iz_high] - dustSNII.Metallicity[iz_low];

      if(deltaz > 0)
        dz = dz / deltaz;
      else
        dz = 0;
    }
  else
    {
      iz_low = 0;
      iz_high = 0;
      dz = 0;
    }

  /* this may be wrong */
  for(imass = ilow; imass < ihigh + 1; imass++)
    stellar_yield[imass] = (1 - dz) * dustSNII.TotalDust_SPH[iz_low][imass] +
      dz * dustSNII.TotalDust_SPH[iz_high][imass];

#ifdef BG_DOUBLE_IMF
  if(phys_dens > All.IMF_PhysDensThresh)
    dust_mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 6);
  else
    dust_mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 2);
#else
  dust_mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 2);
#endif
#endif

  /* yield normalization */
  int HydrogenIndex, HeliumIndex;

  double norm0, norm1;

  HydrogenIndex = element_index("Hydrogen");
  HeliumIndex = element_index("Helium");

  /* zero all negative values */
  for(i = 0; i < BG_NELEMENTS; i++)
    if(metals[i] < 0)
      metals[i] = 0;

  if(mass < 0)
    mass = 0;

  /* get the total mass ejected from the table */
  for(imass = ilow; imass < ihigh + 1; imass++)
    stellar_yield[imass] = (1 - dz) * yieldsSNII.Ejecta_SPH[iz_low][imass] +
      dz * yieldsSNII.Ejecta_SPH[iz_high][imass];

  norm0 = integrate_imf(log_min_mass, log_max_mass, 0.0, 2);

  /* compute the total mass ejected */
  norm1 = mass + metals[HydrogenIndex] + metals[HeliumIndex];
#ifdef BG_DUST
  norm1 += dust_mass;
#endif

  /* normalize the yields */
  if(norm1 > 0)
    {
      for(i = 0; i < BG_NELEMENTS; i++)
	metals_released[i] += metals[i] * (norm0 / norm1);

      *metal_mass_released += mass * (norm0 / norm1);
#ifdef BG_DUST
      *dust_mass_released += dust_mass * (norm0 / norm1);
#endif
    }
  else
    {
      printf("[evolve_SNII] wrong normalization!!!! norm1 = %e\n", norm1);
      endrun(666);
    }
}

/*
 * ----------------------------------------------------------------------
 * This routine computes yields (per unit solar mass) of POPIII stars
 * ----------------------------------------------------------------------
 */

#ifdef BG_POPIII
void evolve_POPIII(double log_min_mass, double log_max_mass, double log_metallicity,
		   double *initial_metals, MyFloat * metals_released,
		   MyFloat * metal_mass_released, MyFloat * dust_mass_released)
{
  int ilow, ihigh, imass;

  int i;

  double metallicity;

#ifdef BG_DUST
  MyFloat dust_mass;
#endif

  MyFloat metals[BG_NELEMENTS], mass;

  if(All.POPIII_MassTransferOn == 0)
    return;

  metallicity = pow(10.0, log_metallicity);

  /* determine integration range, limiting to stars that become POPIII stars */
  if(log_max_mass > log10(All.POPIII_MaxMass_MSUN))
    log_max_mass = log10(All.POPIII_MaxMass_MSUN);

  if(log_min_mass < log10(All.POPIII_MinMass_MSUN))
    log_min_mass = log10(All.POPIII_MinMass_MSUN);

  if(log_min_mass >= log_max_mass)
    return;

  /* determine which mass bins will contribute */
  determine_popiii_imf_bins(log_min_mass, log_max_mass, &ilow, &ihigh);

  /* compute stellar_yield as function of mass */
  for(i = 0; i < BG_NELEMENTS; i++)
    {
      for(imass = ilow; imass < ihigh + 1; imass++)
	popiii_stellar_yield[imass] = yieldsPOPIII.SPH[i][imass] +
	  initial_metals[i] * yieldsPOPIII.Ejecta_SPH[imass];

      metals[i] = integrate_popiii_imf(log_min_mass, log_max_mass, 0.0, 2);
    }

  for(imass = ilow; imass < ihigh + 1; imass++)
    popiii_stellar_yield[imass] = yieldsPOPIII.TotalMetals_SPH[imass] +
      metallicity * yieldsPOPIII.Ejecta_SPH[imass];

  mass = integrate_popiii_imf(log_min_mass, log_max_mass, 0.0, 2);

#ifdef BG_DUST
  for(imass = ilow; imass < ihigh + 1; imass++)
    popiii_stellar_yield[imass] = dustPOPIII.TotalDust_SPH[imass];

  dust_mass = integrate_popiii_imf(log_min_mass, log_max_mass, 0.0, 2);
#endif

  /* yield normalization */
  int HydrogenIndex, HeliumIndex;

  double norm0, norm1;

  HydrogenIndex = element_index("Hydrogen");
  HeliumIndex = element_index("Helium");

  /* zero all negative values */
  for(i = 0; i < BG_NELEMENTS; i++)
    if(metals[i] < 0)
      metals[i] = 0;

  if(mass < 0)
    mass = 0;

  /* get the total mass ejected from the table */
  for(imass = ilow; imass < ihigh + 1; imass++)
    popiii_stellar_yield[imass] = yieldsPOPIII.Ejecta_SPH[imass];

  norm0 = integrate_popiii_imf(log_min_mass, log_max_mass, 0.0, 2);

  /* compute the total mass ejected */
  norm1 = mass + metals[HydrogenIndex] + metals[HeliumIndex];
#ifdef BG_DUST
  norm1 += dust_mass;
#endif

  /* normalize the yields */
  if(norm0 > 0 && norm1 > 0)
    {
      for(i = 0; i < BG_NELEMENTS; i++)
	metals_released[i] += metals[i] * (norm0 / norm1);

      *metal_mass_released += mass * (norm0 / norm1);
#ifdef BG_DUST
      *dust_mass_released += dust_mass * (norm0 / norm1);
#endif
    }
  /*
  else
    {
      printf("[evolve_POPIII] wrong normalization!!!! norm1 = %e\n", norm1);
      printf("[evolve_POPIII] wrong normalization!!!! ilow = %d, ihigh = %d\n", ilow, ihigh);
      printf("[evolve_POPIII] wrong normalization!!!! log(min_mass) = %f, log(max_mass) = %f\n", log_min_mass, log_max_mass);
      endrun(666);
    }
  */

  /* compute the number of SNII from POPIII stellar population */
  if(log_max_mass > log10(All.POPIII_SNII_MaxMass_MSUN))
    log_max_mass = log10(All.POPIII_SNII_MaxMass_MSUN);

  if(log_min_mass < log10(All.POPIII_SNII_MinMass_MSUN))
    log_min_mass = log10(All.POPIII_SNII_MinMass_MSUN);

  if(log_min_mass < log_max_mass)
    NumberOfSNII += integrate_imf(log_min_mass, log_max_mass, 0.0, 0);
}
#endif /* BG_POPIII */

/*
 * ----------------------------------------------------------------------
 * This routine computes the mass of stars dying at some age_of_star_in_Gyr
 * ----------------------------------------------------------------------
 */

double dying_mass_msun(double age_of_star_in_Gyr, double metallicity)
{
  double mass = 0, d_metal, d_time1 = 0, d_time2 = 0, logage, mass1, mass2;

  int i_metal, i, i_time1 = -1, i_time2 = -1;

#ifdef BG_POPIII
  double d_time = 0;

  int i_time = -1;
#endif

  static int first_call = 0;

#ifdef BG_POPIII
  if(metallicity >= All.POPIII_MetallicityThreshold)
    {
#endif
      switch (All.IMF_LifetimeMode)
	{
	case 0:
	  /* padovani & matteucci 1993 */
	  if(first_call == 0 && ThisTask == 0)
	    {
	      if(ThisTask == 0)
		printf("Using Padovani & Matteucci 1993 life-times\n");
	      first_call = 1;
	    }

	  if(age_of_star_in_Gyr > 0.039765318659064693)
	    mass = pow(10, 7.764 - (1.79 - pow(1.338 - 0.1116 * (9 + log10(age_of_star_in_Gyr)), 2)) / 0.2232);
	  else if(age_of_star_in_Gyr > 0.003)
	    mass = pow((age_of_star_in_Gyr - 0.003) / 1.2, -1 / 1.85);
	  else
	    mass = All.IMF_MaxMass_MSUN;
	  break;

	case 1:
	  /* maeder & meynet 1989 */
	  if(first_call == 0 && ThisTask == 0)
	    {
	      if(ThisTask == 0)
		printf("Using Maeder & Meynet 1989 life-times\n");
	      first_call = 1;
	    }
	  
	  if(age_of_star_in_Gyr >= 8.4097378)
	    mass = pow(10, (1 - log10(age_of_star_in_Gyr)) / 0.6545);
	  else if(age_of_star_in_Gyr >= 0.35207776)
	    mass = pow(10, (1.35 - log10(age_of_star_in_Gyr)) / 3.7);
	  else if(age_of_star_in_Gyr >= 0.050931493)
	    mass = pow(10, (0.77 - log10(age_of_star_in_Gyr)) / 2.51);
	  else if(age_of_star_in_Gyr >= 0.010529099)
	    mass = pow(10, (0.17 - log10(age_of_star_in_Gyr)) / 1.78);
	  else if(age_of_star_in_Gyr >= 0.0037734787)
	    mass = pow(10, -(0.94 + log10(age_of_star_in_Gyr)) / 0.86);
	  else if(age_of_star_in_Gyr > 0.003)
	    mass = pow((age_of_star_in_Gyr - 0.003) / 1.2, -0.54054053);
	  else
	    mass = All.IMF_MaxMass_MSUN;
	  break;

	case 2:
	  /* portinari et al. 1998 */
	  if(first_call == 0 && ThisTask == 0)
	    {
	      if(ThisTask == 0)
		printf("Using Portinari et. al 1998 life-times\n");
	      first_call = 1;
	    }

	  if(age_of_star_in_Gyr <= 0)
	    {
	      mass = All.IMF_MaxMass_MSUN;
	      break;
	    }
	  
	  logage = log10(age_of_star_in_Gyr * 1.E9);
	  
	  if(metallicity <= Lifetimes.Metallicity[0])
	    {
	      i_metal = 0;
	      d_metal = 0.0;
	    }
	  else if(metallicity >= Lifetimes.Metallicity[Lifetimes.N_Z - 1])
	    {
	      i_metal = Lifetimes.N_Z - 2;
	      d_metal = 1.0;
	    }
	  else
	    {
	      for(i_metal = 0; i_metal < Lifetimes.N_Z - 1; i_metal++)
		if(Lifetimes.Metallicity[i_metal + 1] > metallicity)
		  break;

	      d_metal = (metallicity - Lifetimes.Metallicity[i_metal]) /
		(Lifetimes.Metallicity[i_metal + 1] - Lifetimes.Metallicity[i_metal]);
	    }

	  if(logage >= Lifetimes.Dyingtime[i_metal][0])
	    {
	      i_time1 = 0;
	      d_time1 = 0.0;
	    }
	  else if(logage <= Lifetimes.Dyingtime[i_metal][Lifetimes.N_MASS - 1])
	    {
	      i_time1 = Lifetimes.N_MASS - 2;
	      d_time1 = 1.0;
	    }
	  
	  if(logage >= Lifetimes.Dyingtime[i_metal + 1][0])
	    {
	      i_time2 = 0;
	      d_time2 = 0.0;
	    }
	  else if(logage <= Lifetimes.Dyingtime[i_metal + 1][Lifetimes.N_MASS - 1])
	    {
	      i_time2 = Lifetimes.N_MASS - 2;
	      d_time2 = 1.0;
	    }
	  
	  i = Lifetimes.N_MASS;
	  while(i >= 0 && (i_time1 == -1 || i_time2 == -1))
	    {
	      i--;
	      if(Lifetimes.Dyingtime[i_metal][i] >= logage && i_time1 == -1)
		{
		  i_time1 = i;
		  d_time1 = (logage - Lifetimes.Dyingtime[i_metal][i_time1]) /
		    (Lifetimes.Dyingtime[i_metal][i_time1 + 1] - Lifetimes.Dyingtime[i_metal][i_time1]);
		}
	      if(Lifetimes.Dyingtime[i_metal + 1][i] >= logage && i_time2 == -1)
		{
		  i_time2 = i;
		  d_time2 = (logage - Lifetimes.Dyingtime[i_metal + 1][i_time2]) /
		    (Lifetimes.Dyingtime[i_metal + 1][i_time2 + 1] - Lifetimes.Dyingtime[i_metal + 1][i_time2]);
		}
	    }
	  
	  mass1 = interpol_1d_dbl(Lifetimes.Mass, i_time1, d_time1);
	  mass2 = interpol_1d_dbl(Lifetimes.Mass, i_time2, d_time2);
	  
	  mass = (1.0 - d_metal) * mass1 + d_metal * mass2;
	  break;
	  
	default:
	  printf("stellar lifetimes not defined\n");
	  endrun(-1);
	}

      if(mass > All.IMF_MaxMass_MSUN)
	mass = All.IMF_MaxMass_MSUN;

      return mass;
#ifdef BG_POPIII
    }
  else
    {
      if(age_of_star_in_Gyr <= 0) /* to avoid log of zero */
	{
	  mass = All.POPIII_IMF_MaxMass_MSUN;

	  return mass;
	}
      else
	{
	  logage = log10(age_of_star_in_Gyr * 1.E9);

	  if(logage >= POPIIILifetimes.Dyingtime[0])
	    {
	      mass = POPIIILifetimes.Mass[0];

              /* not sure the following check is necessary */
              if(mass < All.POPIII_IMF_MinMass_MSUN)
                mass = All.POPIII_IMF_MinMass_MSUN;

	      return mass;
	    }
	  else if(logage <= POPIIILifetimes.Dyingtime[POPIIILifetimes.N_MASS - 1])
	    {
	      mass = POPIIILifetimes.Mass[POPIIILifetimes.N_MASS - 1];

	      /* not sure the following check is necessary */
	      if(mass > All.POPIII_IMF_MaxMass_MSUN)
		mass = All.POPIII_IMF_MaxMass_MSUN;

	      return mass;
	    }

	  i = POPIIILifetimes.N_MASS;
	  while(i >= 0 && i_time == -1)
	    {
	      i--;

	      if(POPIIILifetimes.Dyingtime[i] >= logage && i_time == -1)
		{
		  i_time = i;
		  d_time = (logage - POPIIILifetimes.Dyingtime[i_time]) /
		    (POPIIILifetimes.Dyingtime[i_time + 1] - POPIIILifetimes.Dyingtime[i_time]);
		}
	    }

	  mass = interpol_1d_dbl(POPIIILifetimes.Mass, i_time, d_time);

	  /* not sure the following check is necessary */
	  if(mass < All.POPIII_IMF_MinMass_MSUN)
	    mass = All.POPIII_IMF_MinMass_MSUN;

	  /* not sure the following check is necessary */
	  if(mass > All.POPIII_IMF_MaxMass_MSUN)
	    mass = All.POPIII_IMF_MaxMass_MSUN;

	  return mass;
	}
    }
#endif
}


double lifetime_in_Gyr(double mass, double metallicity)
{
  double time = 0, d_mass, d_metal;

  int i_mass, i_metal;

  static int first_call = 0;

  switch (All.IMF_LifetimeMode)
    {
    case 0:
      /* PM93 (Padovani & Matteucci 1993) */
      if(first_call == 0)
	{
	  if(ThisTask == 0)
	    printf("Using Padovani & Matteucci 1993 life-times\n");
	  first_call = 1;
	}

      if(mass <= 0.6)
	time = 160.0;
      else if(mass <= 6.6)
	time = pow(10.0, (0.334 - sqrt(1.790 - 0.2232 * (7.764 - log10(mass)))) / 0.1116);
      else
	time = 1.2 * pow(mass, -1.85) + 0.003;
      break;

    case 1:
      /* MM89 (Maeder & Meynet 1989) */
      if(first_call == 0)
	{
	  if(ThisTask == 0)
	    printf("Using Maeder & Meynet 1989 life-times\n");
	  first_call = 1;
	}

      if(mass <= 1.3)
	time = pow(10.0, -0.6545 * log10(mass) + 1.0);
      else if(mass <= 3.0)
	time = pow(10.0, -3.7 * log10(mass) + 1.35);
      else if(mass <= 7.0)
	time = pow(10.0, -2.51 * log10(mass) + 0.77);
      else if(mass <= 15.0)
	time = pow(10.0, -1.78 * log10(mass) + 0.17);
      else if(mass <= 60.0)
	time = pow(10.0, -0.86 * log10(mass) - 0.94);
      else
	time = 1.2 * pow(mass, -1.85) + 0.003;
      break;

    case 2:
      /* P98 (Portinari et al. 1998) */
      if(first_call == 0)
	{
	  if(ThisTask == 0)
	    printf("Using Portinari et. al 1998 life-times\n");
	  first_call = 1;
	}

      if(mass <= Lifetimes.Mass[0])
	{
	  i_mass = 0;
	  d_mass = 0.0;
	}
      else if(mass >= Lifetimes.Mass[Lifetimes.N_MASS - 1])
	{
	  i_mass = Lifetimes.N_MASS - 2;
	  d_mass = 1.0;
	}
      else
	{
	  for(i_mass = 0; i_mass < Lifetimes.N_MASS - 1; i_mass++)
	    if(Lifetimes.Mass[i_mass + 1] > mass)
	      break;

	  d_mass = (mass - Lifetimes.Mass[i_mass]) / (Lifetimes.Mass[i_mass + 1] - Lifetimes.Mass[i_mass]);
	}

      if(metallicity <= Lifetimes.Metallicity[0])
	{
	  i_metal = 0;
	  d_metal = 0.0;
	}
      else if(metallicity >= Lifetimes.Metallicity[Lifetimes.N_Z - 1])
	{
	  i_metal = Lifetimes.N_Z - 2;
	  d_metal = 1.0;
	}
      else
	{
	  for(i_metal = 0; i_metal < Lifetimes.N_Z - 1; i_metal++)
	    if(Lifetimes.Metallicity[i_metal + 1] > metallicity)
	      break;

	  d_metal = (metallicity - Lifetimes.Metallicity[i_metal]) /
	    (Lifetimes.Metallicity[i_metal + 1] - Lifetimes.Metallicity[i_metal]);
	}

      /* time in years */
      time = pow(10.0, interpol_2d_dbl(Lifetimes.Dyingtime, i_metal, i_mass, d_metal, d_mass));
      time /= 1.0e9;		/* now in Gyr */

      break;

    default:
      printf("stellar lifetimes not defined\n");
      endrun(-1);
    }

  return time;
}

#endif /* BG_STELLAR_EVOLUTION */


/*
 * ----------------------------------------------------------------------
 * This function returns the elapsed time between a0, a1
 * ----------------------------------------------------------------------
 */

double bg_get_elapsed_time(double a0, double a1, int mode)
{
#define WORKSIZE 100000

  gsl_function F;

  gsl_integration_workspace *workspace;

  double result, abserr, age = 0, t0, t1, factor1, factor2, term1, term2;

  if(All.ComovingIntegrationOn)
    {
      switch (mode)
	{
	case (0):
	  F.function = &bg_time_integ;
	  workspace = gsl_integration_workspace_alloc(WORKSIZE);
	  gsl_integration_qag(&F, a0, a1, 1 / All.Hubble,	/* note: absolute error just a dummy */
			      1.0e-7, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
	  gsl_integration_workspace_free(workspace);

	  age = result * All.UnitTime_in_s / All.HubbleParam;	/* now in seconds */
	  age /= SEC_PER_MEGAYEAR * 1000;	/* now in gigayears */
	  break;

	case (1):
	  if(All.OmegaLambda + All.Omega0 != 1)
	    {
	      printf("[bg_get_elapsed_time] wrong mode = %d, universe is not flat!!!", mode);
	      endrun(666);
	    }

	  factor1 = 2.0 / (3.0 * sqrt(All.OmegaLambda));

	  term1 = sqrt(All.OmegaLambda / All.Omega0) * pow(a0, 1.5);
	  term2 = sqrt(1 + All.OmegaLambda / All.Omega0 * pow(a0, 3));
	  factor2 = log(term1 + term2);

	  t0 = factor1 * factor2;

	  term1 = sqrt(All.OmegaLambda / All.Omega0) * pow(a1, 1.5);
	  term2 = sqrt(1 + All.OmegaLambda / All.Omega0 * pow(a1, 3));
	  factor2 = log(term1 + term2);

	  t1 = factor1 * factor2;

	  result = t1 - t0;

	  age = result / (HUBBLE * All.HubbleParam);	/* now in seconds */
	  age /= SEC_PER_MEGAYEAR * 1000;	/* now in gigayears */
	  break;
	}
    }
  else
    {
      age = (a1 - a0) * All.UnitTime_in_s / All.HubbleParam;	/* now in seconds */
      age /= SEC_PER_MEGAYEAR * 1000;	/* now in gigayears */
    }

  return age;

}


double bg_time_integ(double a, void *param)
{
  double hubble_a;

  hubble_a = All.Omega0 / (a * a * a) + (1 - All.Omega0 - All.OmegaLambda) / (a * a)
#ifdef DARKENERGY
    + DarkEnergy_a(a);
#else
    + All.OmegaLambda;
#endif
  hubble_a = All.Hubble * sqrt(hubble_a);

  return 1 / (a * hubble_a);
}
