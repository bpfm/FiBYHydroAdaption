#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "allvars.h"
#include "bg_dust.h"
#include "bg_proto.h"

struct tdust_params
{
  double T_gas;
  double n_gas;
  double T_rad;
};


#ifdef BG_DUST

/* -------------------------------------------------------------------- */
/*  Dust destruction routine. Two mechanisms are actually implemented:  */
/*  thermal sputtering and sublimation.                                 */
/* -------------------------------------------------------------------- */

void dust_destruction(int j)
{
#if defined(BG_DUST_DESTRUCTION_SUBLIMATION) || defined(BG_DUST_DESTRUCTION_SPUTTERING)
  double T_gas, rho, dest_fraction = 0;
  double x_H_tot, x_He_tot, n_H_tot, n_He_tot;
  double a3inv, hubble_a, time_hubble_a, dt, dtime;

#ifdef BG_DUST_DESTRUCTION_SUBLIMATION
  double T_dust, n_gas, T_rad;
  double da, new_a, old_a;
#endif

  if(All.ComovingIntegrationOn)	/* Factors for comoving integration of hydro */
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = All.Hubble * sqrt(All.Omega0 /
				   (All.Time * All.Time * All.Time)
				   + (1 - All.Omega0 - All.OmegaLambda) /
				   (All.Time * All.Time) + All.OmegaLambda);
      time_hubble_a = All.Time * hubble_a;
    }
  else
    a3inv = time_hubble_a = 1;

  if(All.ComovingIntegrationOn)
    dtime = All.Time * dt / time_hubble_a;
  else
    dtime = dt;

  dtime *= All.UnitTime_in_s / All.HubbleParam;	/* convert to  seconds */

  /* get gas temperature from tables */
  T_gas = bg_get_temperature(j);

  /* convert to physical density */
  rho = SphP[j].d.Density * a3inv;
  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

  /* get total mass fractions of hydrogen and helium */
  x_H_tot = SphP[j].Metals[element_index("Hydrogen")] / P[j].Mass;
  x_He_tot = SphP[j].Metals[element_index("Helium")] / P[j].Mass;

  /* total hydrogen and helium number density */
  n_H_tot = x_H_tot * rho / PROTONMASS;
  n_He_tot = x_He_tot * rho / (4 * PROTONMASS);


  /* ------------------------------ */
  /* thermal sputtering (T > 1e6 K) */
  /* ------------------------------ */

#ifdef BG_DUST_DESTRUCTION_SPUTTERING
  if(T_gas >= THERMAL_SPUTTERING_TEMP_THRESH)
#ifdef BG_DUST_DESTRUCTION_SUBLIMATION
    dest_fraction += THERMAL_SPUTTERING_TIME_SCALE / n_H_tot * SphP[j].DustGrainSize / dtime;
#else
    dest_fraction += THERMAL_SPUTTERING_TIME_SCALE / n_H_tot * All.InitDustGrainSize_nM / dtime;
#endif
#endif

  /* ---------------- */
  /* dust sublimation */
  /* ---------------- */

#ifdef BG_DUST_DESTRUCTION_SUBLIMATION
  /* gas number density */
  n_gas = n_H_tot + n_He_tot;

  /* radiation temperature (only CMB) */
  if(All.ComovingIntegrationOn)
    T_rad = T_CMB0 / All.Time;
  else
    T_rad = T_CMB0;

  /* get dust temperature */
  T_dust = get_dust_temperature(T_gas, n_gas, T_rad);

  if(SUBLIMATION_C2 / T_dust < 100)
    {
      da = -SUBLIMATION_C0 * SUBLIMATION_C1 * exp(-SUBLIMATION_C2 / T_dust) * dtime; /* cm */
      da *= 1e7; /* now in nm */

      old_a = SphP[j].DustGrainSize;
      new_a = old_a + da;

      if(new_a < 0)
	new_a = 0;

      /* variation of average grain size */
      SphP[j].DustGrainSize = new_a;

      /* variation of grain volume at constant density */
      dest_fraction += 1 - (new_a * new_a * new_a) / (old_a * old_a * old_a);
    }
#endif

  if(dest_fraction < 0)
    dest_fraction = 0;

  if(dest_fraction > 1)
    dest_fraction = 1;

  SphP[j].Dusticity *= 1 - dest_fraction;
#endif
}


/* -------------------------------------------------------------------- */
/*  Returns the dust temperature following the formalism of             */
/*  Omukai et al., 2005. Dust temperature is given by equilibrium       */
/*  between collisional cooling and CMB heating.                        */
/* -------------------------------------------------------------------- */

double get_dust_temperature(double T_gas, double n_gas, double T_rad)
{
  int status;
  int iter = 0, max_iter = 100;

  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;

  double T_dust = 0;
  double log_x_lo = 0, log_x_hi = 3;

  gsl_function F;

  struct tdust_params params;

  params.T_gas = T_gas;
  params.n_gas = n_gas;
  params.T_rad = T_rad;

  F.function = &dust_temperature_evaluate;
  F.params = &params;
     
  T = gsl_root_fsolver_brent;
  //T = gsl_root_fsolver_bisection;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, log_x_lo, log_x_hi);

  /*
  printf("using %s method\n",
	  gsl_root_fsolver_name(s));
     
  printf("%5s [%9s, %9s] %9s %12s %9s\n",
	  "iter", "lower", "upper", "root",
	  "temp", "err");
  */

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      T_dust = gsl_root_fsolver_root(s);

      log_x_lo = gsl_root_fsolver_x_lower(s);
      log_x_hi = gsl_root_fsolver_x_upper(s);

      status = gsl_root_test_interval(log_x_lo, log_x_hi,
				      0, 0.02);
     
      if(status == GSL_SUCCESS)
	{
	  printf ("[TASK %2d] %5d [%.7f, %.7f] %.7f %e %.7f\n",
		  ThisTask, iter, log_x_lo, log_x_hi,
		  T_dust, pow(10, T_dust),
		  log_x_hi - log_x_lo);
	}
    }
  while(status == GSL_CONTINUE && iter < max_iter);

  if(iter == max_iter)
    printf("Solver did not converge!!!");

  gsl_root_fsolver_free(s);
     
  return status;
}


double dust_temperature_evaluate(double x, void *params)
{
  struct tdust_params *p = (struct tdust_params *) params;

  double T_gas = p->T_gas;
  double n_gas = p->n_gas;
  double T_rad = p->T_rad;

  double T_dust;

  double v = 5.3e-3;    /* volume fraction: fiducial value v = 5.3e-3.
			  Changing between ~1.e-2 and 1.e-4 implies
			  a factor of ~2 in the dust grain temperature.
			  Changing between ~1.e-1 and 1.e-5 implies
			  a factor of ~5 in the dust grain temperature.
			  (UM, 5/2010)
		       */

  double K0gr = 4.e-4;   /* coefficient of the grain opacity: fiducial value
			   K0gr = 4.e-4. The full law of the grain opacity is
			   Kgr = K0gr*Tgr^2.
			   Changing K0gr between ~1.e-3 and 1.e-5 implies
			   a factor of ~2 in the dust grain temperature.
			  (UM, 5/2010)
			*/

  double Cf = 1.1e-5 * n_gas * v * sqrt(T_gas / 1000) *
    ( 1 - 0.8 * exp(-75 / T_gas) ) / (4 * SIGMA * K0gr);
                       /* cooling prefactor: corresponds to the gas-dust cooling
			  divided by ((T-Tgr)/Tgr^2).
			  (UM, 5/2010)
			*/

  /*--------------------------------------------------------------------------*/

  /* Start calculations (see Omukai et al., 2005). */

  T_dust = pow(10, x);

  return T_dust * T_dust * T_dust * T_dust -
    (T_rad * T_rad * T_rad * T_rad + Cf * (T_gas - T_dust) / (T_dust * T_dust));
}

#endif /* BG_DUST */
