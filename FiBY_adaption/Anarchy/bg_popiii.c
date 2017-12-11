#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "allvars.h"
#include "proto.h"
#include "bg_proto.h"
#include "bg_yields.h"

#ifdef BG_POPIII

void init_popiii_imf(void)
{
  int i;
  double lm_min, lm_max, dlm, norm, lmass, mass;

  popiii_imf_by_number = (double *) mymalloc(POPIII_N_MASS_BIN * sizeof(double));
  popiii_imf_mass_bin = (double *) mymalloc(POPIII_N_MASS_BIN * sizeof(double));
  popiii_imf_mass_bin_log10 = (double *) mymalloc(POPIII_N_MASS_BIN * sizeof(double));
  popiii_integrand = (double *) mymalloc(POPIII_N_MASS_BIN * sizeof(double));
  popiii_stellar_yield = (double *) mymalloc(POPIII_N_MASS_BIN * sizeof(double));

  /* -------------------------------------------------------------------
     define the popIII IMF and normalize it
     ------------------------------------------------------------------- */

  lm_min = log10(All.POPIII_IMF_MinMass_MSUN);
  lm_max = log10(All.POPIII_IMF_MaxMass_MSUN);

  dlm = (lm_max - lm_min) / (double) (POPIII_N_MASS_BIN - 1);

  for(i = 0; i < POPIII_N_MASS_BIN; i++)
    {
      lmass = lm_min + i * dlm;
      mass = pow(10, lmass);

      popiii_imf_by_number[i] = 1.0 / pow(mass, All.POPIII_IMF_Exponent);
      popiii_imf_mass_bin[i] = mass;
      popiii_imf_mass_bin_log10[i] = lmass;
    }

  for(i = 0; i < POPIII_N_MASS_BIN; i++)
    popiii_integrand[i] = popiii_imf_by_number[i] * popiii_imf_mass_bin[i] * popiii_imf_mass_bin[i];

  /* integrate using some rule */
  for(i = 0, norm = 0; i < POPIII_N_MASS_BIN; i++)
    norm += popiii_integrand[i];

  /* correct first and last bins */
  norm -= 0.5 * integrand[0];
  norm -= 0.5 * integrand[POPIII_N_MASS_BIN - 1];

  /* multiply by the bin size */
  norm *= dlm * log(10.0);    /* log(10) since mass function tabulated as function of log_10(mass) */

  /* normalize the imf */
  for(i = 0; i < POPIII_N_MASS_BIN; i++)
    popiii_imf_by_number[i] /= norm;
}


double integrate_popiii_imf(double log_min_dying_mass, double log_max_dying_mass, double m2, int mode)
{
  double result;
  int ilow, ihigh, index;
  double dlm, dm;

  dlm = popiii_imf_mass_bin_log10[1] - popiii_imf_mass_bin_log10[0];	/* dlog(m) */

  determine_popiii_imf_bins(log_min_dying_mass, log_max_dying_mass, &ilow, &ihigh);

  if(mode == 0)
    /* integrate by number */
    for(index = ilow; index < ihigh + 1; index++)
      popiii_integrand[index] = popiii_imf_by_number[index] * popiii_imf_mass_bin[index];
  else if(mode == 1)
    /* integrate by mass */
    for(index = ilow; index < ihigh + 1; index++)
      popiii_integrand[index] = popiii_imf_by_number[index] * popiii_imf_mass_bin[index] * popiii_imf_mass_bin[index];
  else if(mode == 2)
    /* integrate number * yield weighted */
    for(index = ilow; index < ihigh + 1; index++)
      popiii_integrand[index] = popiii_imf_by_number[index] * popiii_stellar_yield[index] * popiii_imf_mass_bin[index];
  else
    {
      printf("invalid mode in integrate_imf = %d\n", mode);
      endrun(-1);
    }

  /* integrate using some rule */
  result = 0;
  for(index = ilow; index < ihigh + 1; index++)
    result += popiii_integrand[index];

  result = result - 0.5 * popiii_integrand[ilow] - 0.5 * popiii_integrand[ihigh];

  /* correct first bin */
  dm = (log_min_dying_mass - popiii_imf_mass_bin_log10[ilow]) / dlm;

  if(dm < 0.5)
    result -= dm * popiii_integrand[ilow];
  else
    {
      result -= 0.5 * popiii_integrand[ilow];
      result -= (dm - 0.5) * popiii_integrand[ilow + 1];
    }

  /* correct last bin */
  dm = (log_max_dying_mass - popiii_imf_mass_bin_log10[ihigh - 1]) / dlm;
  if(dm < 0.5)
    {
      result -= 0.5 * popiii_integrand[ihigh];
      result -= (0.5 - dm) * popiii_integrand[ihigh - 1];
    }
  else
    result -= (1 - dm) * popiii_integrand[ihigh];

  /* multiply by the bin size */
  result *= dlm * log(10.0);	/* log(10) since mass function tabulated as function of log_10(mass) */

  return result;
}


void determine_popiii_imf_bins(double log_min_dying_mass, double log_max_dying_mass, int *ilow, int *ihigh)
{
  int i1, i2;

  if(log_min_dying_mass < popiii_imf_mass_bin_log10[0])
    log_min_dying_mass = popiii_imf_mass_bin_log10[0];

  if(log_min_dying_mass > popiii_imf_mass_bin_log10[POPIII_N_MASS_BIN - 1])
    log_min_dying_mass = popiii_imf_mass_bin_log10[POPIII_N_MASS_BIN - 1];

  if(log_max_dying_mass < popiii_imf_mass_bin_log10[0])
    log_max_dying_mass = popiii_imf_mass_bin_log10[0];

  if(log_max_dying_mass > popiii_imf_mass_bin_log10[POPIII_N_MASS_BIN - 1])
    log_max_dying_mass = popiii_imf_mass_bin_log10[POPIII_N_MASS_BIN - 1];

  for(i1 = 0; i1 < POPIII_N_MASS_BIN - 2 && popiii_imf_mass_bin_log10[i1 + 1] < log_min_dying_mass; i1++);

  for(i2 = 1; i2 < POPIII_N_MASS_BIN - 1 && popiii_imf_mass_bin_log10[i2] < log_max_dying_mass; i2++);

  *ilow = i1;
  *ihigh = i2;
}


#ifdef BG_POPIII_BH_SEEDS
void popiii_black_holes()
{
  int i, converted, tot_converted;

  double metallicity, age_of_star_in_Gyr;


  converted = 0;

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 4)
	{
#ifdef BG_METALSMOOTHING
	  metallicity = StarP[P[i].StarID].MetallicitySmoothed;
#else
	  metallicity = StarP[P[i].StarID].Metallicity;
#endif

	  if(metallicity < All.POPIII_MetallicityThreshold) /* we have a POPIII particle */
	    {
	      age_of_star_in_Gyr = bg_get_elapsed_time(StarP[P[i].StarID].StarBirthTime, All.Time, 1); /* note: this is valid for a flat universe! */

	      if(dying_mass_msun(age_of_star_in_Gyr, metallicity) <= All.POPIII_MinMass_MSUN) /* end of life of POPIII star, turn it into BH */
		{
		  /* printf("[TASK %d] ID = %d, converted into BH\n", ThisTask, P[i].ID); */

		  P[i].Type = 5;
#ifdef BLACK_HOLES
		  P[i].BH_BirthTime = All.Time;
#ifdef BH_THERMALFEEDBACK
		  P[i].BH_Energy = 0;
#endif
		  P[i].b1.BH_Density = 1e-10;
		  P[i].b2.BH_Entropy = 0.0;
		  P[i].b3.BH_SurroundingGasVel[0] = 0.0;
		  P[i].b3.BH_SurroundingGasVel[1] = 0.0;
		  P[i].b3.BH_SurroundingGasVel[2] = 0.0;
		  P[i].BH_Mass = P[i].Mass; //All.SeedBHMassOverGasMass * All.OrigGasMass;
		  P[i].BH_Mdot = 0;
#endif
		  /* rearrange stars */
		  if(P[i].StarID < N_star - 1)
		    {
		      StarP[P[i].StarID] = StarP[N_star - 1];
		      P[StarP[N_star - 1].PID].StarID = P[i].StarID;
		    }

		  converted++;
		  N_star--;
		}
	    }
	}
    }

  MPI_Allreduce(&converted, &tot_converted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(tot_converted)
    All.TotN_star -= tot_converted;
}
#endif

#endif /* BG_POPIII */
