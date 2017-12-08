#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#include "bg_proto.h"
#include "bg_molecules.h"

#ifdef BG_MOL_NETWORK

#ifdef BG_MOL_NETWORK_TABULATED_RATES
/* rate coefficient tables declaration */
double *ttab;

double *k1tab;
double *k2tab;
double *k3tab;
double *k4tab;
double *k5tab;
double *k6tab;
double *k7tab;
double *k8tab;
double *k9tab;

double *k10tab;
double *k11tab;
double *k12tab;
double *k13tab;
double *k14tab;
double *k15tab;
double *k16tab;
double *k17tab;
double *k18tab;
double *k19tab;

double *k20tab;
double *k21tab;
double *k22tab;
double *k23tab;
double *k24tab;

double *d1tab;
double *d2tab;
double *d3tab;
double *d4tab;
double *d5tab;
double *d6tab;
#endif /* BG_MOL_NETWORK_TABULATED_RATES */

/* interpolated rate coefficients */
double kk1;
double kk2;
double kk3;
double kk4;
double kk5;
double kk6;
double kk7;
double kk8;
double kk9;

double kk10;
double kk11;
double kk12;
double kk13;
double kk14;
double kk15;
double kk16;
double kk17;
double kk18;
double kk19;

double kk20;
double kk21;
double kk22;
double kk23;
double kk24;

double dd1;
double dd2;
double dd3;
double dd4;
double dd5;
double dd6;

double k_diss_hm;
double k_diss_h2;

/* debug */
double min_chem_tstep, min_chem_tstep_dens, min_chem_tstep_temp;
double min_chem_tstep_H, min_chem_tstep_H2, min_chem_tstep_e;
double min_chem_tstep_dt, min_chem_tstep_nH, min_chem_tstep_nH2, min_chem_tstep_ne;
/* debug */


/* ============================================================================ */
/*  THE MOLECULAR NETWORK MAIN FUNCTION                                         */
/* ============================================================================ */

void bg_molecules()
{
  int i;

  double H2_mass_old = 0, tot_H2_mass_old;
  double H2_mass_new = 0, tot_H2_mass_new;

/* debug */
  min_chem_tstep = 1e30;
  min_chem_tstep_dens = 0;
  min_chem_tstep_temp = 0;
  min_chem_tstep_dt = 0;
  min_chem_tstep_H = 0;
  min_chem_tstep_H2 = 0;
  min_chem_tstep_e = 0;
/* debug */

  if(ThisTask == 0)
    {
      printf("Start molecular network...\n");
      fflush(stdout);
    }

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0 && SphP[i].OnEOS == 1)
	SphP[i].Chem_dt = 1e30;

      if(P[i].Type == 0 && SphP[i].OnEOS != 1)
	{
	  k_diss_hm = k_diss_h2 = 0;

#if defined(LW_BACKGROUND) || defined(LW_LOCAL)
	  k_diss_hm = SphP[i].Kdiss_Hm;
	  k_diss_h2 = SphP[i].Kdiss_H2;
#endif

	  H2_mass_old += SphP[i].x_H2 * P[i].Mass;
	  bg_molecules_evaluate(i);
	  H2_mass_new += SphP[i].x_H2 * P[i].Mass;
	}
    }

  MPI_Allreduce(&H2_mass_old, &tot_H2_mass_old, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&H2_mass_new, &tot_H2_mass_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("produced H2 mass: %e solar masses\n",
	     (tot_H2_mass_new - tot_H2_mass_old) * All.UnitMass_in_g / All.HubbleParam / SOLAR_MASS);
      printf("done with molecular network\n");
      fflush(stdout);
    }

  /* debug */
  for(i = 0; i < NTask; i++)
    {
      if(i == ThisTask && min_chem_tstep < 1e30) /* print only if chemistry sets dt */
	{
	  printf("[%3d] chem_tstep=%e, dt=%e, dt_H=%e, dt_H2=%e, dt_e=%e, dens=%e, temp=%e\n",
		 ThisTask,
		 min_chem_tstep,
		 min_chem_tstep_dt,
		 min_chem_tstep_H,
		 min_chem_tstep_H2,
		 min_chem_tstep_e,
		 min_chem_tstep_dens,
		 min_chem_tstep_temp);

	  printf("[%3d] x_H=%e, x_H2=%e, n_e=%e\n",
		 ThisTask,
		 min_chem_tstep_nH,
		 min_chem_tstep_nH2,
		 min_chem_tstep_ne);
	}

      MPI_Barrier(MPI_COMM_WORLD);
    }
  /* debug */
}


/* ============================================================================ */
/*  THE ACTUAL MOLECULAR NETWORK FUNCTION                                       */
/* ============================================================================ */

void bg_molecules_evaluate(int i)
{
  double n_e, n_e_guess;

  double n_H, n_Hp, n_Hm, n_H2, n_H2p;
  double n_H_guess, n_Hp_guess, n_Hm_guess, n_H2_guess, n_H2p_guess;

  double n_He, n_Hep, n_Hepp;
  double n_He_guess, n_Hep_guess, n_Hepp_guess;

  double n_D, n_Dp, n_HD;
  double n_D_guess, n_Dp_guess, n_HD_guess;

  double rho, temp, norm_factor;

  double x_H, x_H_tot, x_He, x_He_tot, x_D, x_D_tot;

  double n_H_tot, n_He_tot, n_D_tot, min_number_density;

  double dt, dtime, hubble_a, time_hubble_a, a3inv;
  double dt_chem = 1e30, dt_H = 1e30, dt_H2 = 1e30, dt_e = 1e30;


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
    }
  else
    a3inv = time_hubble_a = 1;


  /* ============================================================================ */
  /*  Physical quantities for particle i                                          */
  /* ============================================================================ */

  /* convert to physical density */
  rho = SphP[i].d.Density * a3inv;
  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

  /* get temperature of the gas */
  temp = bg_get_temperature(i);


  /* ============================================================================ */
  /*  Check of temperature threshold                                              */
  /* ============================================================================ */

  /* do it only for cold gas T < MOL_NETWORK_TEMP_THRESH */
  if(temp > MOL_NETWORK_TEMP_THRESH)
    {
      /* get total mass fractions of hydrogen, helium and deuterium */
      x_H_tot = SphP[i].Metals[element_index("Hydrogen")] / P[i].Mass;
      x_He_tot = SphP[i].Metals[element_index("Helium")] / P[i].Mass;
      x_D_tot = All.InitAbundance_D; /* this is constant over the all simulation */

      SphP[i].x_Hp = 0.9 * x_H_tot;
      SphP[i].x_Hm = 3e-9 * x_H_tot;
      SphP[i].x_H2 = 4e-8 * x_H_tot;
      SphP[i].x_H2p = 8e-8 * x_H_tot;

      SphP[i].x_Hep = 0.8 * x_He_tot;
      SphP[i].x_Hepp = 1e-2 * x_He_tot;

      SphP[i].x_Dp = 0.9 * x_D_tot;
      SphP[i].x_HD = 1e-8 * x_D_tot;

      SphP[i].n_e = 0.9 * rho / PROTONMASS;

      SphP[i].Chem_dt = 1e30;

      return;
    }

  if(P[i].TimeBin < 5 && P[i].TimeBin > 0)
    {
      x_H_tot = SphP[i].Metals[element_index("Hydrogen")] / P[i].Mass;
      printf("[MOL TimeBin=%d] x_H=%e, x_Hp=%e, x_H2=%e, n_e=%e, T=%e\n",
	     P[i].TimeBin,
	     x_H_tot,
	     SphP[i].x_Hp / x_H_tot,
	     SphP[i].x_H2 / x_H_tot,
	     SphP[i].n_e, temp);
    }

#ifdef BG_MOL_NETWORK_TABULATED_RATES
  if(temp <= MOL_NETWORK_MIN_TEMP)
    return;
#endif

  /* ============================================================================ */
  /*  Actual time step in seconds                                                 */
  /* ============================================================================ */

  if(P[i].TimeBin == 0) /* skip at start-up (dt = 0) */
    {
      dt = 0;
    }
  else
    {
      dt = TISTEP(P[i].TimeBin) * All.Timebase_interval;

      /* the actual time-step */
      if(All.ComovingIntegrationOn)
	dtime = All.Time * dt / time_hubble_a;
      else
	dtime = dt;

      /* now in seconds */
      dt = dtime * All.UnitTime_in_s / All.HubbleParam;
    }


  /* ============================================================================ */
  /*  Compute rate coefficients                                                   */
  /* ============================================================================ */

#ifdef BG_MOL_NETWORK_TABULATED_RATES
  temp = log10(temp); /* note: for table interpolation log space is used */
#endif

  set_rate_coefficients(temp); /* set the rate coefficients for the given temperature:
				  rate coefficients are computed only once for the given
				  temperature to speed up the network; they are given
				  either by functions or interpolated on the tables. */


  /* ============================================================================ */
  /*  Initial number densities                                                    */
  /* ============================================================================ */

  /* get total mass fractions of hydrogen, helium and deuterium */
  x_H_tot = SphP[i].Metals[element_index("Hydrogen")] / P[i].Mass;
  x_He_tot = SphP[i].Metals[element_index("Helium")] / P[i].Mass;
  x_D_tot = All.InitAbundance_D; /* this is constant over the all simulation */

  /* total hydrogen, helium and deuterium number density */
  n_H_tot = x_H_tot * rho / PROTONMASS;
  n_He_tot = x_He_tot * rho / (4 * PROTONMASS);
  n_D_tot = x_D_tot * rho / (2 * PROTONMASS);

  /* convert mass fractions into number densities - particles store the mass fractions of ions/molecules */
  n_H = (x_H_tot - SphP[i].x_Hp - SphP[i].x_Hm - SphP[i].x_H2 - SphP[i].x_H2p - SphP[i].x_HD / 3) * rho / PROTONMASS;   /* 1/3 of the HD mass fraction is considered here */
  n_Hp = SphP[i].x_Hp * rho / PROTONMASS;
  n_Hm = SphP[i].x_Hm * rho / PROTONMASS;
  n_H2 = SphP[i].x_H2 * rho / (2 * PROTONMASS);
  n_H2p = SphP[i].x_H2p * rho / (2 * PROTONMASS);

  n_He = (x_He_tot - SphP[i].x_Hep - SphP[i].x_Hepp) * rho / (4 * PROTONMASS);
  n_Hep = SphP[i].x_Hep * rho / (4 * PROTONMASS);
  n_Hepp = SphP[i].x_Hepp * rho / (4 * PROTONMASS);

  n_D = (x_D_tot - SphP[i].x_Dp - SphP[i].x_HD * 2 / 3) * rho / (2 * PROTONMASS);  /* 2/3 of the HD mass fraction are considered here */
  n_Dp = SphP[i].x_Dp * rho / (2 * PROTONMASS);
  n_HD = SphP[i].x_HD * rho / (3 * PROTONMASS);

  /* set the guess for the solver: n_X is the current values;
     n_X_guess will be the value at the end of the step */
  n_e = n_e_guess = SphP[i].n_e;

  n_H_guess = n_H;
  n_Hp_guess = n_Hp;
  n_Hm_guess = n_Hm;
  n_H2_guess = n_H2;
  n_H2p_guess = n_H2p;

  n_He_guess = n_He;
  n_Hep_guess = n_Hep;
  n_Hepp_guess = n_Hepp;

  n_D_guess = n_D;
  n_Dp_guess = n_Dp;
  n_HD_guess = n_HD;

  /* define the minimum number density for this particle (fraction of the total H number density) */
  min_number_density = MOL_NETWORK_SMALL_NUMBER * n_H_tot;


  /* ============================================================================ */
  /*  Integration over dt                                                         */
  /* ============================================================================ */

  /* slow evolving reactions */
  n_H_guess = (prod_H(n_H_guess, n_Hp_guess, n_Hm_guess, n_H2_guess, n_H2p_guess, n_He_guess, n_D_guess, n_e_guess) * dt + n_H) /
    (1 + dest_H(n_H_guess, n_Hp_guess, n_Hm_guess, n_H2_guess, n_H2p_guess, n_Hep_guess, n_Dp_guess, n_HD_guess, n_e_guess) * dt);

  if(n_H_guess < min_number_density) /* n_H */
    n_H_guess = min_number_density;

  n_Hp_guess = (prod_Hp(n_H_guess, n_H2_guess, n_H2p_guess, n_Hep_guess, n_Dp_guess, n_e_guess) * dt + n_Hp) /
    (1 + dest_Hp(n_H_guess, n_Hm_guess, n_H2_guess, n_He_guess, n_D_guess, n_HD_guess, n_e_guess) * dt);

  if(n_Hp_guess < min_number_density) /* n_H(+) */
    n_Hp_guess = min_number_density;

  n_He_guess = (prod_He(n_H_guess, n_Hep_guess, n_e_guess) * dt + n_He) /
    (1 + dest_He(n_Hp_guess, n_e_guess) * dt);

  if(n_He_guess < min_number_density) /* n_He */
    n_He_guess = min_number_density;

  n_Hep_guess = (prod_Hep(n_Hp_guess, n_He_guess, n_Hepp_guess, n_e_guess) * dt + n_Hepp) /
    (1 + dest_Hep(n_H_guess, n_e_guess) * dt);

  if(n_Hep_guess < min_number_density) /* n_He(+) */
    n_Hep_guess = min_number_density;

  n_Hepp_guess = (prod_Hepp(n_Hep_guess, n_e_guess) * dt + n_Hepp) /
    (1 + dest_Hepp(n_e_guess) * dt);

  if(n_Hepp_guess < min_number_density) /* n_He(++) */
    n_Hepp_guess = min_number_density;

  n_D_guess = (prod_D(n_H_guess, n_Dp_guess, n_HD_guess) * dt + n_D) /
    (1 + dest_D(n_Hp_guess, n_H2_guess) * dt);

  if(n_D_guess < min_number_density) /* n_D */
    n_D_guess = min_number_density;

  n_Dp_guess = (prod_Dp(n_Hp_guess, n_D_guess, n_HD_guess) * dt + n_Dp) /
    (1 + dest_Dp(n_H_guess, n_H2_guess) * dt);

  if(n_Dp_guess < min_number_density) /* n_D(+) */
    n_Dp_guess = min_number_density;

  n_HD_guess = (prod_HD(n_H2_guess, n_D_guess, n_Dp_guess) * dt + n_HD) /
    (1 + dest_HD(n_H_guess, n_Hp_guess) * dt);

  if(n_HD_guess < min_number_density) /* n_HD */
    n_HD_guess = min_number_density;

  n_e_guess = (prod_e(n_H_guess, n_Hp_guess, n_Hm_guess, n_H2_guess, n_He_guess, n_Hep_guess, n_e_guess) * dt + n_e) / 
    (1 + dest_e(n_H_guess, n_Hp_guess, n_Hm_guess, n_H2_guess, n_H2p_guess, n_He_guess, n_Hep_guess, n_Hepp_guess, n_e_guess) * dt);

  if(n_e_guess < min_number_density) /* n_e */
    n_e_guess = min_number_density;

  /* fast evolving reactions */
  n_Hm_guess = (prod_Hm(n_H_guess, n_e_guess) * dt + n_Hm) /
    (1 + dest_Hm(n_H_guess, n_Hp_guess, n_H2p_guess, n_e_guess) * dt);

  if(n_Hm_guess < min_number_density) /* n_H(-) */
    n_Hm_guess = min_number_density;

  n_H2p_guess = (prod_H2p(n_H_guess, n_Hp_guess, n_Hm_guess, n_H2_guess) * dt + n_H2p) /
    (1 + dest_H2p(n_H_guess, n_Hm_guess, n_e_guess) * dt);

  if(n_H2p_guess < min_number_density) /* n_H2(+) */
    n_H2p_guess = min_number_density;

  /* H2 value */
  n_H2_guess = (prod_H2(n_H_guess, n_Hp_guess, n_Hm_guess, n_H2_guess, n_H2p_guess, n_HD_guess) * dt + n_H2) /
    (1 + dest_H2(n_H_guess, n_Hp_guess, n_H2_guess, n_D_guess, n_Dp_guess, n_e_guess) * dt);

  if(n_H2_guess < min_number_density) /* n_H2 */
    n_H2_guess = min_number_density;


  /* ============================================================================ */
  /*  Normalization of number densities                                           */
  /* ============================================================================ */

  /* normalize the number density fractions */
  norm_factor = n_H_tot / (n_H_guess + n_Hp_guess + n_Hm_guess + 2 * n_H2_guess + 2 * n_H2p_guess + n_HD_guess);

  n_H = n_H_guess* norm_factor;
  n_Hp = n_Hp_guess* norm_factor;
  n_Hm = n_Hm_guess* norm_factor;
  n_H2 = n_H2_guess* norm_factor;
  n_H2p = n_H2p_guess* norm_factor;

  norm_factor = n_He_tot / (n_He_guess + n_Hep_guess + n_Hepp_guess);

  n_He = n_He_guess* norm_factor;
  n_Hep = n_Hep_guess* norm_factor;
  n_Hepp = n_Hepp_guess* norm_factor;

  norm_factor = n_D_tot / (n_D_guess + n_Dp_guess + n_HD_guess);

  n_D = n_D_guess* norm_factor;
  n_Dp = n_Dp_guess* norm_factor;
  n_HD = n_HD_guess * norm_factor;

  /* recompute electron number density after normalization */
  n_e = free_electrons(n_Hp, n_Hm, n_H2p, n_Hep, n_Hepp, n_Dp);


  /* ============================================================================ */
  /*  Once abundances are updated compute the chemical dt                         */
  /* ============================================================================ */

#ifdef BG_MOL_NETWORK_TIMESTEP
#ifdef BG_MOL_NETWORK_H_TIMESTEP
  /* H time step = EPS * dH/(dH/dt) */
  dt_H = EPS * n_H /
    fabs((prod_H(n_H, n_Hp, n_Hm, n_H2, n_H2p, n_He, n_D, n_e) -
	  dest_H(n_H, n_Hp, n_Hm, n_H2, n_H2p, n_Hep, n_Dp, n_HD, n_e) * n_H));

  if(dt_chem > dt_H && dt_H > 0)
    dt_chem = dt_H;
#endif

#ifdef BG_MOL_NETWORK_H2_TIMESTEP
  /* H2 time step = EPS * dH2/(dH2/dt) */
  dt_H2 = EPS * n_H2 /
    fabs(prod_H2(n_H, n_Hp, n_Hm, n_H2, n_H2p, n_HD) -
	 dest_H2(n_H, n_Hp, n_H2, n_D, n_Dp, n_e) * n_H2);

  if(dt_chem > dt_H2 && dt_H2 > 0)
    dt_chem = dt_H2;
#endif

  /* e time step = EPS * de/(de/dt) */
  dt_e = EPS * n_e /
    fabs((prod_e(n_H, n_Hp, n_Hm, n_H2, n_He, n_Hep, n_e) -
	  dest_e(n_H, n_Hp, n_Hm, n_H2, n_H2p, n_He, n_Hep, n_Hepp, n_e) * n_e));

  if(dt_chem > dt_e && dt_e > 0)
    dt_chem = dt_e;
#endif

  if(dt_chem < MIN_CHEM_TSTEP_YR * SEC_PER_YEAR)
    dt_chem = MIN_CHEM_TSTEP_YR * SEC_PER_YEAR;

  /* debug */
  if(min_chem_tstep > dt_chem && dt_chem < dt) /* update only if chemistry sets dt */
    {
      min_chem_tstep = dt_chem;
      min_chem_tstep_dens = rho;
      min_chem_tstep_temp = temp; /* note: this is log10(temp) if tables are used */
      min_chem_tstep_dt = dt;
      min_chem_tstep_H = dt_H;
      min_chem_tstep_H2 = dt_H2;
      min_chem_tstep_e = dt_e;
      min_chem_tstep_nH = n_H / n_H_tot;
      min_chem_tstep_nH2 = n_H2 / n_H_tot;
      min_chem_tstep_ne = n_e;
    }
  /* debug */


  /* ============================================================================ */
  /*  Update particle properties                                                  */
  /* ============================================================================ */

  /* store electron number density */
  SphP[i].n_e = n_e;

  /* convert back to mass fractions */
  x_H = n_H * PROTONMASS / rho;
  SphP[i].x_Hp = n_Hp * PROTONMASS / rho;
  SphP[i].x_Hm = n_Hm * PROTONMASS / rho;
  SphP[i].x_H2 = n_H2 * (2 * PROTONMASS) / rho;
  SphP[i].x_H2p = n_H2p * (2 * PROTONMASS) / rho;

  x_He = n_He * (4 * PROTONMASS) / rho;
  SphP[i].x_Hep = n_Hep * (4 * PROTONMASS) / rho;
  SphP[i].x_Hepp = n_Hepp * (4 * PROTONMASS) / rho;

  x_D = n_D * (2 * PROTONMASS) / rho;
  SphP[i].x_Dp = n_Dp * (2 * PROTONMASS) / rho;
  SphP[i].x_HD = n_HD * (3 * PROTONMASS) / rho;

  /* store chemical timestep */
  SphP[i].Chem_dt = dt_chem;

  /* just a check for debugging with the debugger */
  x_H_tot = x_H + SphP[i].x_Hp + SphP[i].x_Hm + SphP[i].x_H2 + SphP[i].x_H2p + SphP[i].x_HD / 3;
  x_He_tot = x_He + SphP[i].x_Hep + SphP[i].x_Hepp;
  x_D_tot = x_D + SphP[i].x_Dp + SphP[i].x_HD * 2 / 3;
}


/* ============================================================================ */
/*  PRODUCTION/DESTRUCTION FUNCTIONS FOR H2                                     */
/* ============================================================================ */

/* --------------------------------------------------------------- */
/*  H production and destruction                                   */
/* --------------------------------------------------------------- */

double prod_H(double n_H, double n_Hp, double n_Hm, double n_H2, double n_H2p, double n_He, double n_D, double n_e)
{
  return kk2 * n_Hp * n_e +
    3 * kk11 * n_H2 * n_H +
    kk12 * n_H2 * n_Hp +
    2 * kk13 * n_H2 * n_e +
    kk15 * n_He * n_Hp +
    kk16 * n_Hm * n_e +
    2 * kk17 * n_Hm * n_H +
    2 * kk18 * n_Hm * n_Hp +
    2 * kk20 * n_H2p * n_e +
    kk21 * n_H2p * n_Hm +
    kk22 * n_H * n_H * n_H +
    2 * kk24 * n_H2 * n_H2 +
    dd1 * n_H2 * n_D +
    dd5 * n_Hp * n_D;
}

double dest_H(double n_H, double n_Hp, double n_Hm, double n_H2, double n_H2p, double n_Hep, double n_Dp, double n_HD, double n_e)
{
  return kk1 * n_e +
    kk7 * n_e +
    kk8 * n_Hm +
    kk9 * n_Hp +
    kk10 * n_H2p +
    kk11 * n_H2 +
    kk14 * n_Hep +
    kk17 * n_Hm +
    3 * kk22 * n_H * n_H +
    2 * kk23 * n_H * n_H2 +
    dd3 * n_HD +
    dd6 * n_Dp;
}

/* --------------------------------------------------------------- */
/*  H(+) production and destruction                                */
/* --------------------------------------------------------------- */

double prod_Hp(double n_H, double n_H2, double n_H2p, double n_Hep, double n_Dp, double n_e)
{
  return kk1 * n_H * n_e +
    kk10 * n_H2p * n_H +
    kk14 * n_Hep * n_H +
    dd2 * n_Dp * n_H2 +
    dd6 * n_Dp * n_H;
}

double dest_Hp(double n_H, double n_Hm, double n_H2, double n_He, double n_D, double n_HD, double n_e)
{
  return kk2 * n_e +
    kk9 * n_H +
    kk12 * n_H2 +
    kk15 * n_He +
    kk18 * n_Hm +
    kk19 * n_Hm +
    dd4 * n_HD +
    dd5 * n_D;
}

/* --------------------------------------------------------------- */
/*  H(-) production and destruction                                */
/* --------------------------------------------------------------- */

double prod_Hm(double n_H, double n_e)
{
  return kk7 * n_H * n_e;
}

double dest_Hm(double n_H, double n_Hp, double n_H2p, double n_e)
{
  return kk8 * n_H +
    kk16 * n_e +
    kk17 * n_H +
    kk18 * n_Hp +
    kk19 * n_Hp +
    kk21 * n_H2p +
    k_diss_hm;
}

/* --------------------------------------------------------------- */
/*  H2 production and destruction                                  */
/* --------------------------------------------------------------- */

double prod_H2(double n_H, double n_Hp, double n_Hm, double n_H2, double n_H2p, double n_HD)
{
  return kk8 * n_Hm * n_H +
    kk10 * n_H2p * n_H +
    kk21 * n_H2p * n_Hm +
    kk22 * n_H * n_H * n_H +
    2 * kk23 * n_H * n_H * n_H2 +
    kk24 * n_H2 * n_H2 +
    dd3 * n_HD * n_H +
    dd4 * n_HD * n_Hp;
}

double dest_H2(double n_H, double n_Hp, double n_H2, double n_D, double n_Dp, double n_e)
{
  return kk11 * n_H +
    kk12 * n_Hp +
    kk13 * n_e +
    kk23 * n_H * n_H +
    2 * kk24 * n_H2 +
    dd1 * n_D +
    dd2 * n_Dp +
    k_diss_h2;
}

/* --------------------------------------------------------------- */
/*  H2(+) production and destruction                               */
/* --------------------------------------------------------------- */

double prod_H2p(double n_H, double n_Hp, double n_Hm, double n_H2)
{
  return kk9 * n_H * n_Hp +
    kk12 * n_H2 * n_Hp +
    kk19 * n_Hm * n_Hp;
}

double dest_H2p(double n_H, double n_Hm, double n_e)
{
  return kk10 * n_H +
    kk20 * n_e +
    kk21 * n_Hm;
}

/* --------------------------------------------------------------- */
/*  He production and destruction                                  */
/* --------------------------------------------------------------- */

double prod_He(double n_H, double n_Hep, double n_e)
{
  return kk4 * n_Hep * n_e +
    kk14 * n_Hep * n_H;
}

double dest_He(double n_Hp, double n_e)
{
  return kk3 * n_e +
    kk15 * n_Hp;
}

/* --------------------------------------------------------------- */
/*  He(+) production and destruction                               */
/* --------------------------------------------------------------- */

double prod_Hep(double n_Hp, double n_He, double n_Hepp, double n_e)
{
  return kk3 * n_He * n_e +
    kk6 * n_Hepp * n_e +
    kk15 * n_He * n_Hp;
}

double dest_Hep(double n_H, double n_e)
{
  return kk4 * n_e +
    kk5 * n_e +
    kk14 * n_H;
}

/* --------------------------------------------------------------- */
/*  He(++) production and destruction                              */
/* --------------------------------------------------------------- */

double prod_Hepp(double n_Hep, double n_e)
{
  return kk5 * n_Hep * n_e;
}

double dest_Hepp(double n_e)
{
  return kk6 * n_e;
}

/* --------------------------------------------------------------- */
/*  e production and destruction                                   */
/* --------------------------------------------------------------- */

double prod_e(double n_H, double n_Hp, double n_Hm, double n_H2, double n_He, double n_Hep, double n_e)
{
  return 2 * kk1 * n_H * n_e +
    2 * kk3 * n_He * n_e +
    2 * kk5 * n_Hep * n_e +
    kk8 * n_Hm * n_H +
    kk13 * n_H2 * n_e +
    2 * kk16 * n_Hm * n_e +
    kk17 * n_Hm * n_H +
    kk19 * n_Hm * n_Hp;
}

double dest_e(double n_H, double n_Hp, double n_Hm, double n_H2, double n_H2p, double n_He, double n_Hep, double n_Hepp, double n_e)
{
  return kk1 * n_H +
    kk2 * n_Hp +
    kk3 * n_He +
    kk4 * n_Hep +
    kk5 * n_Hep +
    kk6 * n_Hepp +
    kk7 * n_H +
    kk13 * n_H2 +
    kk16 * n_Hm +
    kk20 * n_H2p;
}


/* ============================================================================ */
/*  CHARGE CONSERVATION: RETURNS THE ELECTRON NUMBER DENSITY                    */
/* ============================================================================ */

/* --------------------------------------------------------------- */
/*  charge conservation                                            */
/* --------------------------------------------------------------- */

double free_electrons(double n_Hp, double n_Hm, double n_H2p, double n_Hep, double n_Hepp, double n_Dp)
{
  return n_Hp - n_Hm + n_H2p + n_Hep + 2 * n_Hepp + n_Dp;
}


/* ============================================================================ */
/*  RATE COEFFICIENTS FOR H2 NETWORK                                            */
/* ============================================================================ */

/* --------------------------------------------------------------- */
/*  H + e --> H(+) + 2e                                            */
/* --------------------------------------------------------------- */

double k1(double temp)
{
#define A0 -32.71396786
#define A1 13.536556
#define A2 -5.73932875
#define A3 1.56315498
#define A4 -0.2877056
#define A5 3.48255977e-2
#define A6 -2.6319761e-3
#define A7 1.11954395e-4
#define A8 -2.03914984e-6

  double log_temp_ev1, log_temp_ev2, log_temp_ev3, log_temp_ev4;
  double log_temp_ev5, log_temp_ev6, log_temp_ev7, log_temp_ev8;

  if(temp < 2800)
    return MOL_NETWORK_SMALL_NUMBER;

  log_temp_ev1 = log(temp / EV_IN_K);

  log_temp_ev2 = log_temp_ev1 * log_temp_ev1;
  log_temp_ev3 = log_temp_ev1 * log_temp_ev2;
  log_temp_ev4 = log_temp_ev1 * log_temp_ev3;
  log_temp_ev5 = log_temp_ev1 * log_temp_ev4;
  log_temp_ev6 = log_temp_ev1 * log_temp_ev5;
  log_temp_ev7 = log_temp_ev1 * log_temp_ev6;
  log_temp_ev8 = log_temp_ev1 * log_temp_ev7;

  return exp(A0 +
	     A1 * log_temp_ev1 +
	     A2 * log_temp_ev2 +
	     A3 * log_temp_ev3 +
	     A4 * log_temp_ev4 +
	     A5 * log_temp_ev5 +
	     A6 * log_temp_ev6 +
	     A7 * log_temp_ev7 +
	     A8 * log_temp_ev8);
}

/* --------------------------------------------------------------- */
/*  H(+) + e --> H + photon                                        */
/* --------------------------------------------------------------- */

double k2(double temp)
{
#define B0 -28.6130338
#define B1 -0.72411256
#define B2 -2.02604473e-2
#define B3 -2.38086188e-3
#define B4 -3.21260521e-4
#define B5 -1.42150291e-5
#define B6 4.98910892e-6
#define B7 5.75561414e-7
#define B8 -1.85676704e-8
#define B9 -3.07113524e-9

  double log_temp_ev1, log_temp_ev2, log_temp_ev3, log_temp_ev4;
  double log_temp_ev5, log_temp_ev6, log_temp_ev7, log_temp_ev8;
  double log_temp_ev9;

  log_temp_ev1 = log(temp / EV_IN_K);

  log_temp_ev2 = log_temp_ev1 * log_temp_ev1;
  log_temp_ev3 = log_temp_ev1 * log_temp_ev2;
  log_temp_ev4 = log_temp_ev1 * log_temp_ev3;
  log_temp_ev5 = log_temp_ev1 * log_temp_ev4;
  log_temp_ev6 = log_temp_ev1 * log_temp_ev5;
  log_temp_ev7 = log_temp_ev1 * log_temp_ev6;
  log_temp_ev8 = log_temp_ev1 * log_temp_ev7;
  log_temp_ev9 = log_temp_ev1 * log_temp_ev8;

  return exp(B0 +
	     B1 * log_temp_ev1 +
	     B2 * log_temp_ev2 +
	     B3 * log_temp_ev3 +
	     B4 * log_temp_ev4 +
	     B5 * log_temp_ev5 +
	     B6 * log_temp_ev6 +
	     B7 * log_temp_ev7 +
	     B8 * log_temp_ev8 +
	     B9 * log_temp_ev9);
}

/* --------------------------------------------------------------- */
/*  He + e --> He(+) + 2e                                          */
/* --------------------------------------------------------------- */

double k3(double temp)
{
#define C0 -44.09864886
#define C1 23.91596563
#define C2 -10.7532302
#define C3 3.05803875
#define C4 -5.6851189e-1
#define C5 6.79539123e-2
#define C6 -5.00905610e-3
#define C7 2.06723616e-4
#define C8 -3.64916141e-6

  double log_temp_ev1, log_temp_ev2, log_temp_ev3, log_temp_ev4;
  double log_temp_ev5, log_temp_ev6, log_temp_ev7, log_temp_ev8;

  if(temp < 5600)
    return MOL_NETWORK_SMALL_NUMBER;

  log_temp_ev1 = log(temp / EV_IN_K);

  log_temp_ev2 = log_temp_ev1 * log_temp_ev1;
  log_temp_ev3 = log_temp_ev1 * log_temp_ev2;
  log_temp_ev4 = log_temp_ev1 * log_temp_ev3;
  log_temp_ev5 = log_temp_ev1 * log_temp_ev4;
  log_temp_ev6 = log_temp_ev1 * log_temp_ev5;
  log_temp_ev7 = log_temp_ev1 * log_temp_ev6;
  log_temp_ev8 = log_temp_ev1 * log_temp_ev7;

  return exp(C0 +
	     C1 * log_temp_ev1 +
	     C2 * log_temp_ev2 +
	     C3 * log_temp_ev3 +
	     C4 * log_temp_ev4 +
	     C5 * log_temp_ev5 +
	     C6 * log_temp_ev6 +
	     C7 * log_temp_ev7 +
	     C8 * log_temp_ev8);
}

/* --------------------------------------------------------------- */
/*  He(+) + e --> He + photon                                      */
/* --------------------------------------------------------------- */

double k4(double temp)
{
#define D0 3.925e-13
#define D1 -0.6353

#define E0 1.544e-9
#define E1 -1.5
#define E2 -48.596
#define E3 0.3
#define E4 8.10

  double temp_ev, k4r, k4d;

  temp_ev = temp / EV_IN_K;

  k4r = D0 * pow(temp_ev, D1);

  if(temp > 18230)
    k4d = E0 * pow(temp_ev, E1) * exp(E2 / temp_ev) * (E3 + exp(E4 / temp_ev));
  else
    k4d = 0;

  return k4r + k4d;
}

/* --------------------------------------------------------------- */
/*  He(+) + e --> He(++) + 2e                                      */
/* --------------------------------------------------------------- */

double k5(double temp)
{
#define F0 -68.71040990
#define F1 43.93347633
#define F2 -18.4806699
#define F3 4.70162649
#define F4 -7.6924663e-1
#define F5 8.113042e-2
#define F6 -5.32402063e-3
#define F7 1.97570531e-4
#define F8 -3.16558106e-6

  double log_temp_ev1, log_temp_ev2, log_temp_ev3, log_temp_ev4;
  double log_temp_ev5, log_temp_ev6, log_temp_ev7, log_temp_ev8;

  if(temp < 11500)
    return MOL_NETWORK_SMALL_NUMBER;

  log_temp_ev1 = log(temp / EV_IN_K);

  log_temp_ev2 = log_temp_ev1 * log_temp_ev1;
  log_temp_ev3 = log_temp_ev1 * log_temp_ev2;
  log_temp_ev4 = log_temp_ev1 * log_temp_ev3;
  log_temp_ev5 = log_temp_ev1 * log_temp_ev4;
  log_temp_ev6 = log_temp_ev1 * log_temp_ev5;
  log_temp_ev7 = log_temp_ev1 * log_temp_ev6;
  log_temp_ev8 = log_temp_ev1 * log_temp_ev7;

  return exp(F0 +
	     F1 * log_temp_ev1 +
	     F2 * log_temp_ev2 +
	     F3 * log_temp_ev3 +
	     F4 * log_temp_ev4 +
	     F5 * log_temp_ev5 +
	     F6 * log_temp_ev6 +
	     F7 * log_temp_ev7 +
	     F8 * log_temp_ev8);
}

/* --------------------------------------------------------------- */
/*  He(++) + e --> He(+) + photon                                  */
/* --------------------------------------------------------------- */

double k6(double temp)
{
  return 2 * k2(temp / 4);
}

/* --------------------------------------------------------------- */
/*  H + e --> H(-) + photon                                        */
/* --------------------------------------------------------------- */

double k7(double temp)
{
#define G0 1.4e-18
#define G1 0.928
#define G2 16200.0

    return G0 * pow(temp, G1) * exp(-temp / G2);
}

/* --------------------------------------------------------------- */
/*  H(-) + H --> H2 + e                                            */
/* --------------------------------------------------------------- */

double k8(double temp)
{
#define H0 4.0e-9
#define H1 -0.17

  return H0 * pow(temp, H1);
}

/* --------------------------------------------------------------- */
/*  H + H(+) --> H2(+) + photon                                    */
/* --------------------------------------------------------------- */

double k9(double temp)
{
#define I0 -19.38
#define I1 -1.523
#define I2 1.118
#define I3 -0.1269

  double log_temp1, log_temp2, log_temp3;

  log_temp1 = log10(temp);

  log_temp2 = log_temp1 * log_temp1;
  log_temp3 = log_temp1 * log_temp2;

  return pow(10, I0 +
	     I1 * log_temp1 +
	     I2 * log_temp2 +
	     I3 * log_temp3);
}

/* --------------------------------------------------------------- */
/*  H2(+) + H --> H2(*) + H(+)                                     */
/* --------------------------------------------------------------- */

double k10(void)
{
#define I4 6e-10

  return I4;
}

/* --------------------------------------------------------------- */
/*  H2 + H --> 3H                                                  */
/* --------------------------------------------------------------- */

double k11(double temp)
{
#define J0 1.067e-10
#define J1 2.012
#define J2 -4.463
#define J3 0.2472
#define J4 3.512

  double temp_ev;

  temp_ev = temp / EV_IN_K;

  if(temp > 1240 && temp < 20000)
    return J0 * pow(temp_ev, J1) * exp(J2 / temp_ev) / pow(1 + J3 * temp_ev, J4);
  else
    return MOL_NETWORK_SMALL_NUMBER;
}

/* --------------------------------------------------------------- */
/*  H2 + H(+) --> H2(+) + H                                        */
/* --------------------------------------------------------------- */

double k12(double temp)
{
#define U0 -21237.15
#define U1 -3.3232183e-7
#define U2 3.3735382e-7
#define U3 -1.4491368e-7
#define U4 3.4172805e-8
#define U5 -4.7813720e-9
#define U6 3.9731542e-10
#define U7 -1.8171411e-11
#define U8 3.5311932e-13

  double ln_temp1, ln_temp2, ln_temp3, ln_temp4;
  double ln_temp5, ln_temp6, ln_temp7;

  if(temp < 460)
    return MOL_NETWORK_SMALL_NUMBER;

  ln_temp1 = log(temp);

  ln_temp2 = ln_temp1 * ln_temp1;
  ln_temp3 = ln_temp1 * ln_temp2;
  ln_temp4 = ln_temp1 * ln_temp3;
  ln_temp5 = ln_temp1 * ln_temp4;
  ln_temp6 = ln_temp1 * ln_temp5;
  ln_temp7 = ln_temp1 * ln_temp6;

  return exp(U0 / temp) *
    (U1 +
     U2 * ln_temp1 +
     U3 * ln_temp2 +
     U4 * ln_temp3 +
     U5 * ln_temp4 +
     U6 * ln_temp5 +
     U7 * ln_temp6 +
     U8 * ln_temp7);
}

/* --------------------------------------------------------------- */
/*  H2 + e --> 2H + e                                              */
/* --------------------------------------------------------------- */

double k13(double temp)
{
#define K0 3.73e-9
#define K1 0.1121
#define K2 -99430.0

  if(temp < 1960)
    return MOL_NETWORK_SMALL_NUMBER;
  else
    return K0 * pow(temp, K1) * exp(K2 / temp);
}

/* --------------------------------------------------------------- */
/*  He(+) + H --> He + H(+) + photon                               */
/* --------------------------------------------------------------- */

double k14(double temp)
{
#define L0 1.20e-15
#define L1 300.0
#define L2 0.25

  return L0 * pow(temp / L1, L2);
}

/* --------------------------------------------------------------- */
/*  He + H(+) --> He(+) + H                                        */
/* --------------------------------------------------------------- */

double k15(double temp)
{
#define M0 1.26e-9
#define M1 -0.75
#define M2 -127500.0

#define N0 4.0e-37
#define N1 4.74

  if(temp < 3000)
    return MOL_NETWORK_SMALL_NUMBER;

  if(temp < 10000)
    return M0 * pow(temp, M1) * exp(M2 / temp);
  else
    return N0 * pow(temp, N1);
}

/* --------------------------------------------------------------- */
/*  H(-) + e --> H + 2e                                            */
/* --------------------------------------------------------------- */

double k16(double temp)
{
#define O0 -18.01849334
#define O1 2.3608522
#define O2 -0.28274430
#define O3 1.62331664e-2
#define O4 -3.36501203e-2
#define O5 1.17832978e-2
#define O6 -1.65619470e-3
#define O7 1.06827520e-4
#define O8 -2.63128581e-6

  double log_temp_ev1, log_temp_ev2, log_temp_ev3, log_temp_ev4;
  double log_temp_ev5, log_temp_ev6, log_temp_ev7, log_temp_ev8;

  if(temp < 185)
    return MOL_NETWORK_SMALL_NUMBER;

  log_temp_ev1 = log(temp / EV_IN_K);

  log_temp_ev2 = log_temp_ev1 * log_temp_ev1;
  log_temp_ev3 = log_temp_ev1 * log_temp_ev2;
  log_temp_ev4 = log_temp_ev1 * log_temp_ev3;
  log_temp_ev5 = log_temp_ev1 * log_temp_ev4;
  log_temp_ev6 = log_temp_ev1 * log_temp_ev5;
  log_temp_ev7 = log_temp_ev1 * log_temp_ev6;
  log_temp_ev8 = log_temp_ev1 * log_temp_ev7;

  return exp(O0 +
	     O1 * log_temp_ev1 +
	     O2 * log_temp_ev2 +
	     O3 * log_temp_ev3 +
	     O4 * log_temp_ev4 +
	     O5 * log_temp_ev5 +
	     O6 * log_temp_ev6 +
	     O7 * log_temp_ev7 +
	     O8 * log_temp_ev8);
}

/* --------------------------------------------------------------- */
/*  H(-) + H --> 2H + e                                            */
/* --------------------------------------------------------------- */

double k17(double temp)
{
#define P0 -20.37260896
#define P1 1.13944933
#define P2 -1.4210135e-1
#define P3 8.4644554e-3
#define P4 -1.4327641e-3
#define P5 2.0122503e-4
#define P6 8.6639632e-5
#define P7 -2.5850097e-5
#define P8 2.4555012e-6
#define P9 -8.0683825e-8

#define P10 2.5634e-9
#define P11 1.78186

  double log_temp_ev1, log_temp_ev2, log_temp_ev3, log_temp_ev4;
  double log_temp_ev5, log_temp_ev6, log_temp_ev7, log_temp_ev8;
  double log_temp_ev9;

  log_temp_ev1 = log(temp / EV_IN_K);

  if(temp / EV_IN_K > 0.1)
    {
      log_temp_ev2 = log_temp_ev1 * log_temp_ev1;
      log_temp_ev3 = log_temp_ev1 * log_temp_ev2;
      log_temp_ev4 = log_temp_ev1 * log_temp_ev3;
      log_temp_ev5 = log_temp_ev1 * log_temp_ev4;
      log_temp_ev6 = log_temp_ev1 * log_temp_ev5;
      log_temp_ev7 = log_temp_ev1 * log_temp_ev6;
      log_temp_ev8 = log_temp_ev1 * log_temp_ev7;
      log_temp_ev9 = log_temp_ev1 * log_temp_ev8;

      return exp(P0 +
		 P1 * log_temp_ev1 +
		 P2 * log_temp_ev2 +
		 P3 * log_temp_ev3 +
		 P4 * log_temp_ev4 +
		 P5 * log_temp_ev5 +
		 P6 * log_temp_ev6 +
		 P7 * log_temp_ev7 +
		 P8 * log_temp_ev8 +
		 P9 * log_temp_ev9);
    }
  else
    {
      return P10 * pow(temp / EV_IN_K, P11);
    }
}

/* --------------------------------------------------------------- */
/*  H(-) + H(+) --> 2H                                             */
/* --------------------------------------------------------------- */

double k18(double temp)
{
#define Q0 6.3e-8
#define Q1 5.7e-6
#define Q2 -9.2e-11
#define Q3 4.4e-13

  return Q0 +
    Q1 / sqrt(temp) +
    Q2 * sqrt(temp) +
    Q3 * temp;
}

/* --------------------------------------------------------------- */
/*  H(-) + H(+) --> H2(+) + e                                      */
/* --------------------------------------------------------------- */

double k19(double temp)
{
#define R0 4.0e-4
#define R1 -1.4
#define R2 -15100.0

#define S0 1.0e-8
#define S1 -0.4

#define T_THRESH 1.0e4

  if(temp > T_THRESH)
    return R0 * pow(temp, R1) * exp(R2 / temp);
  else
    return S0 * pow(temp, S1);
}

/* --------------------------------------------------------------- */
/*  H2(+) + e --> 2H                                               */
/* --------------------------------------------------------------- */

double k20(double temp)
{
#define V0 2e-7

  return V0 / sqrt(temp);
}

/* --------------------------------------------------------------- */
/*  H2(+) + H(-) --> H + H2                                        */
/* --------------------------------------------------------------- */

double k21(double temp)
{
#define W0 5e-7
#define W1 100.0

  return W0 / sqrt(temp / W1);
}

/* --------------------------------------------------------------- */
/*  3H --> H2 + H                                                  */
/* --------------------------------------------------------------- */

double k22(double temp)
{
#define X0 5e-29

  return X0 / temp;
}

/* --------------------------------------------------------------- */
/*  2H + H2 --> 2H2                                                */
/* --------------------------------------------------------------- */

double k23(double temp)
{
  return k22(temp) / 8;
}

/* --------------------------------------------------------------- */
/*  H2 + H2 --> 2H + H2                                            */
/* --------------------------------------------------------------- */

double k24(double temp)
{
#define T0 8.125e-8
#define T1 -52000.0
#define T2 -6000.0

  if(temp < 1050)
    return MOL_NETWORK_SMALL_NUMBER;
  else
    return T0 / sqrt(temp) * exp(T1 / temp) * (1 - exp(T2 / temp));
}


/* ============================================================================ */
/*  PRODUCTION/DESTRUCTION FUNCTIONS FOR HD                                     */
/* ============================================================================ */

/* --------------------------------------------------------------- */
/*  D production and destruction                                   */
/* --------------------------------------------------------------- */

double prod_D(double n_H, double n_Dp, double n_HD)
{
  return dd3 * n_HD * n_H +
    dd6 * n_H * n_Dp;
}

double dest_D(double n_Hp, double n_H2)
{
  return dd1 * n_H2 +
    dd5 * n_Hp;
}

/* --------------------------------------------------------------- */
/*  D(+) production and destruction                                */
/* --------------------------------------------------------------- */

double prod_Dp(double n_Hp, double n_D, double n_HD)
{
  return dd4 * n_Hp * n_HD +
    dd5 * n_Hp * n_D;
}

double dest_Dp(double n_H, double n_H2)
{
  return dd2 * n_H2 +
    dd6 * n_H;
}

/* --------------------------------------------------------------- */
/*  HD production and destruction                                  */
/* --------------------------------------------------------------- */

double prod_HD(double n_H2, double n_D, double n_Dp)
{
  return dd1 * n_H2 * n_D +
    dd2 * n_H2 * n_Dp;
}

double dest_HD(double n_H, double n_Hp)
{
  return dd3 * n_H +
    dd4 * n_Hp;
}


/* ============================================================================ */
/*  RATE COEFFICIENTS FOR HD NETWORK                                            */
/* ============================================================================ */

/* --------------------------------------------------------------- */
/*  H2 + D --> HD + H                                              */
/* --------------------------------------------------------------- */

double d1(double temp)
{
#define DA0 9.0e-11
#define DA1 -3876.0

  if(temp > 83)
    return DA0 * exp(DA1 / temp);
  else
    return MOL_NETWORK_SMALL_NUMBER;
}

/* --------------------------------------------------------------- */
/*  H2 + D(+) --> HD + H(+)                                        */
/* --------------------------------------------------------------- */

double d2(double temp)
{
#define DB0 1.6e-9

  return DB0;
}

/* --------------------------------------------------------------- */
/*  HD + H --> D + H2                                              */
/* --------------------------------------------------------------- */

double d3(double temp)
{
#define DC0 3.2e-11
#define DC1 -3624.0

  if(temp > 80)
    return DC0 * exp(DC1 / temp);
  else
    return MOL_NETWORK_SMALL_NUMBER;
}

/* --------------------------------------------------------------- */
/*  HD + H(+) --> D(+) + H2                                        */
/* --------------------------------------------------------------- */

double d4(double temp)
{
#define DD0 1.0e-9
#define DD1 -464.0

  if(temp > 10)
    return DD0 * exp(DD1 / temp);
  else
    return MOL_NETWORK_SMALL_NUMBER;
}

/* --------------------------------------------------------------- */
/*  H(+) + D --> H + D(+)                                          */
/* --------------------------------------------------------------- */

double d5(double temp)
{
#define DE0 2.0e-10
#define DE1 0.402
#define DE2 -37.1
#define DE3 -3.31e-17
#define DE4 1.48

  if(temp > 2.5)
    return DE0 * pow(temp, DE1) * exp(DE2 / temp) +
      DE3 * pow(temp, DE4);
  else
    return MOL_NETWORK_SMALL_NUMBER;
}

/* --------------------------------------------------------------- */
/*  H + D(+) --> H(+) + D                                          */
/* --------------------------------------------------------------- */

double d6(double temp)
{
#define DF0 2.06e-10
#define DF1 0.396
#define DF2 -33.0
#define DF3 2.03e-9
#define DF4 -0.332

  return DF0 * pow(temp, DF1) * exp(DF2 / temp) +
    DF3 * pow(temp, DF4);
}


/* ============================================================================ */
/*  RATE COEFFICIENT TABLES                                                     */
/* ============================================================================ */

#ifdef BG_MOL_NETWORK_TABULATED_RATES
/* --------------------------------------------------------------- */
/*  Initialization of rate coefficient tabels                      */
/* --------------------------------------------------------------- */

void init_rate_coefficient_tables(void)
{
  int i;

  double log_min_temp, log_max_temp, dlog_temp, temp;

  ttab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));

  k1tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k2tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k3tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k4tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k5tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k6tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k7tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k8tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k9tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));

  k10tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k11tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k12tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k13tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k14tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k15tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k16tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k17tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k18tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k19tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));

  k20tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k21tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k22tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k23tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  k24tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));

  d1tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  d2tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  d3tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  d4tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  d5tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));
  d6tab = (double *) mymalloc(MOL_NETWORK_TABLE_SIZE * sizeof(double));

  log_min_temp = log10(MOL_NETWORK_MIN_TEMP);
  log_max_temp = log10(MOL_NETWORK_MAX_TEMP);

  dlog_temp = (log_max_temp - log_min_temp) / (MOL_NETWORK_TABLE_SIZE - 1);

  for(i = 0; i < MOL_NETWORK_TABLE_SIZE; i++)
    {
      temp = pow(10.0, log_min_temp + i * dlog_temp);

      ttab[i] = log_min_temp + i * dlog_temp;

      k1tab[i] = log10(k1(temp));
      k2tab[i] = log10(k2(temp));
      k3tab[i] = log10(k3(temp));
      k4tab[i] = log10(k4(temp));
      k5tab[i] = log10(k5(temp));
      k6tab[i] = log10(k6(temp));
      k7tab[i] = log10(k7(temp));
      k8tab[i] = log10(k8(temp));
      k9tab[i] = log10(k9(temp));

      k10tab[i] = log10(k10());
      k11tab[i] = log10(k11(temp));
      k12tab[i] = log10(k12(temp));
      k13tab[i] = log10(k13(temp));
      k14tab[i] = log10(k14(temp));
      k15tab[i] = log10(k15(temp));
      k16tab[i] = log10(k16(temp));
      k17tab[i] = log10(k17(temp));
      k18tab[i] = log10(k18(temp));
      k19tab[i] = log10(k19(temp));

      k20tab[i] = log10(k20(temp));
      k21tab[i] = log10(k21(temp));
      k22tab[i] = log10(k22(temp));
      k23tab[i] = log10(k23(temp));
      k24tab[i] = log10(k24(temp));

      d1tab[i] = log10(d1(temp));
      d2tab[i] = log10(d2(temp));
      d3tab[i] = log10(d3(temp));
      d4tab[i] = log10(d4(temp));
      d5tab[i] = log10(d5(temp));
      d6tab[i] = log10(d6(temp));
    }

  /* debug */
  if(ThisTask == 0)
    {
      output_rates();
      output_rate_tables();
    }
  /* debug */
}
#endif /* BG_MOL_NETWORK_TABULATED_RATES */


/* ============================================================================ */
/*  FUNCTION FOR INTERPOLATION OF RATE COEFFICIENTS                             */
/* ============================================================================ */

#ifdef BG_MOL_NETWORK_TABULATED_RATES
/* --------------------------------------------------------------- */
/*  Interpolation of rate coefficient tabels                      */
/* --------------------------------------------------------------- */

double k_table_interpol(double *table, double log_temp)
{
  int i = 0;

  double d_temp = 0;

  /*  if(log_temp >= ttab[0] && log_temp < ttab[MOL_NETWORK_TABLE_SIZE - 1]) */
  if(log_temp >= MOL_NETWORK_MIN_TEMP_LOG10 && log_temp < MOL_NETWORK_MAX_TEMP_LOG10)
    {
      /* i = floor((log_temp - ttab[0]) / MOL_NETWORK_DLOGTEMP); */
      i = floor((log_temp - MOL_NETWORK_MIN_TEMP_LOG10) * MOL_NETWORK_DLOGTEMP_M1);

      if(i == MOL_NETWORK_TABLE_SIZE - 1)
	{
	  i--;
	  d_temp = 1;
	}
      else
	d_temp = (log_temp - ttab[i]) * MOL_NETWORK_DLOGTEMP_M1;
    }
  else
    {
      printf("[MOL_NETWORK] Out of table bounds!!! Task=%d, log_temp=%f, min_log_temp=%f, max_log_temp=%f\n",
	     ThisTask, log_temp, MOL_NETWORK_MIN_TEMP_LOG10, MOL_NETWORK_MAX_TEMP_LOG10);

      endrun(30035);
    }

  return pow(10.0, (table[i] + d_temp * (table[i + 1] - table[i])));
}
#endif /* BG_MOL_NETWORK_TABULATED_RATES */


/* ============================================================================ */
/*  FUNCTION FOR COMPUTATION OF RATE COEFFICIENT AT GIVEN TEMPERATURE           */
/* ============================================================================ */

void set_rate_coefficients(double temp)
{
#ifdef BG_MOL_NETWORK_TABULATED_RATES
  int i;

  double d_temp;

  if(temp >= MOL_NETWORK_MIN_TEMP_LOG10 && temp < MOL_NETWORK_MAX_TEMP_LOG10)
    {
      i = floor((temp - MOL_NETWORK_MIN_TEMP_LOG10) * MOL_NETWORK_DLOGTEMP_M1);

      if(i == MOL_NETWORK_TABLE_SIZE - 1)
	{
	  i--;
	  d_temp = 1;
	}
      else
	d_temp = (temp - ttab[i]) * MOL_NETWORK_DLOGTEMP_M1;
    }
  else
    {
      printf("[MOL_NETWORK] Out of table bounds!!! Task=%d, log_temp=%f, min_log_temp=%f, max_log_temp=%f\n",
	     ThisTask, temp, MOL_NETWORK_MIN_TEMP_LOG10, MOL_NETWORK_MAX_TEMP_LOG10);

      endrun(30035);
    }

  /* rate coefficients for given temperature from tables */
  kk1 = pow(10.0, (k1tab[i] + d_temp * (k1tab[i + 1] - k1tab[i])));
  kk2 = pow(10.0, (k2tab[i] + d_temp * (k2tab[i + 1] - k2tab[i])));
  kk3 = pow(10.0, (k3tab[i] + d_temp * (k3tab[i + 1] - k3tab[i])));
  kk4 = pow(10.0, (k4tab[i] + d_temp * (k4tab[i + 1] - k4tab[i])));
  kk5 = pow(10.0, (k5tab[i] + d_temp * (k5tab[i + 1] - k5tab[i])));
  kk6 = pow(10.0, (k6tab[i] + d_temp * (k6tab[i + 1] - k6tab[i])));
  kk7 = pow(10.0, (k7tab[i] + d_temp * (k7tab[i + 1] - k7tab[i])));
  kk8 = pow(10.0, (k8tab[i] + d_temp * (k8tab[i + 1] - k8tab[i])));
  kk9 = pow(10.0, (k9tab[i] + d_temp * (k9tab[i + 1] - k9tab[i])));

  kk10 = k10(); /* note: k10 returns a constant value */
  kk11 = pow(10.0, (k11tab[i] + d_temp * (k11tab[i + 1] - k11tab[i])));
  kk12 = pow(10.0, (k12tab[i] + d_temp * (k12tab[i + 1] - k12tab[i])));
  kk13 = pow(10.0, (k13tab[i] + d_temp * (k13tab[i + 1] - k13tab[i])));
  kk14 = pow(10.0, (k14tab[i] + d_temp * (k14tab[i + 1] - k14tab[i])));
  kk15 = pow(10.0, (k15tab[i] + d_temp * (k15tab[i + 1] - k15tab[i])));
  kk16 = pow(10.0, (k16tab[i] + d_temp * (k16tab[i + 1] - k16tab[i])));
  kk17 = pow(10.0, (k17tab[i] + d_temp * (k17tab[i + 1] - k17tab[i])));
  kk18 = pow(10.0, (k18tab[i] + d_temp * (k18tab[i + 1] - k18tab[i])));
  kk19 = pow(10.0, (k19tab[i] + d_temp * (k19tab[i + 1] - k19tab[i])));

  kk20 = pow(10.0, (k20tab[i] + d_temp * (k20tab[i + 1] - k20tab[i])));
  kk21 = pow(10.0, (k21tab[i] + d_temp * (k21tab[i + 1] - k21tab[i])));
  kk22 = pow(10.0, (k22tab[i] + d_temp * (k22tab[i + 1] - k22tab[i])));
  kk23 = pow(10.0, (k23tab[i] + d_temp * (k23tab[i + 1] - k23tab[i])));
  kk24 = pow(10.0, (k24tab[i] + d_temp * (k24tab[i + 1] - k24tab[i])));

  dd1 = pow(10.0, (d1tab[i] + d_temp * (d1tab[i + 1] - d1tab[i])));
  dd2 = pow(10.0, (d2tab[i] + d_temp * (d2tab[i + 1] - d2tab[i])));
  dd3 = pow(10.0, (d3tab[i] + d_temp * (d3tab[i + 1] - d3tab[i])));
  dd4 = pow(10.0, (d4tab[i] + d_temp * (d4tab[i + 1] - d4tab[i])));
  dd5 = pow(10.0, (d5tab[i] + d_temp * (d5tab[i + 1] - d5tab[i])));
  dd6 = pow(10.0, (d6tab[i] + d_temp * (d6tab[i + 1] - d6tab[i])));

#else /* BG_MOL_NETWORK_TABULATED_RATES */

  /* rate coefficients for given temperature from functions */
  kk1 = k1(temp);
  kk2 = k2(temp);
  kk3 = k3(temp);
  kk4 = k4(temp);
  kk5 = k5(temp);
  kk6 = k6(temp);
  kk7 = k7(temp);
  kk8 = k8(temp);
  kk9 = k9(temp);

  kk10 = k10(); /* note: k10 returns a constant value */
  kk11 = k11(temp);
  kk12 = k12(temp);
  kk13 = k13(temp);
  kk14 = k14(temp);
  kk15 = k15(temp);
  kk16 = k16(temp);
  kk17 = k17(temp);
  kk18 = k18(temp);
  kk19 = k19(temp);

  kk20 = k20(temp);
  kk21 = k21(temp);
  kk22 = k22(temp);
  kk23 = k23(temp);
  kk24 = k24(temp);

  dd1 = d1(temp);
  dd2 = d2(temp);
  dd3 = d3(temp);
  dd4 = d4(temp);
  dd5 = d5(temp);
  dd6 = d6(temp);
#endif /* BG_MOL_NETWORK_TABULATED_RATES */
}

/* debug */
void output_rates(void)
{
#define N_TEMP_BINS 200
#define MIN_TEMP_BIN 1
#define MAX_TEMP_BIN 3e4

  int i;

  double dlogtemp, temp;

  FILE *FdRateCoeff;

  char file[200];


  sprintf(file, "%s/rate_coefficients.data", All.OutputDir);
  if(!(FdRateCoeff = fopen(file, "w")))
    {
      printf("error in opening file '%s'\n", file);
      endrun(1);
    }

  dlogtemp = (log10(MAX_TEMP_BIN) - log10(MIN_TEMP_BIN)) / N_TEMP_BINS;

  for(i = 0; i < N_TEMP_BINS; i++)
    {
      temp = pow(10, log10(MIN_TEMP_BIN) + i * dlogtemp);

      fprintf(FdRateCoeff, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
             temp,
             k1(temp),
             k2(temp),
             k3(temp),
             k4(temp),
             k5(temp),
             k6(temp),
             k7(temp),
             k8(temp),
             k9(temp),
             k10(),
             k11(temp),
             k12(temp),
             k13(temp),
             k14(temp),
             k15(temp),
             k16(temp),
             k17(temp),
             k18(temp),
             k19(temp),
             k20(temp),
             k21(temp),
             k22(temp),
             k23(temp),
             k24(temp),
             d1(temp),
             d2(temp),
             d3(temp),
             d4(temp),
             d5(temp),
             d6(temp));
    }

  fclose(FdRateCoeff);
}
/* debug */

#ifdef BG_MOL_NETWORK_TABULATED_RATES
/*debug */
void output_rate_tables(void)
{
  int i;

  FILE *FdRateCoeff;

  char file[200];

  sprintf(file, "%s/rate_coefficients.table", All.OutputDir);
  if(!(FdRateCoeff = fopen(file, "w")))
    {
      printf("error in opening file '%s'\n", file);
      endrun(1);
    }

  for(i = 0; i < MOL_NETWORK_TABLE_SIZE; i++)
    {
      fprintf(FdRateCoeff, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
             ttab[i],
             k1tab[i],
             k2tab[i],
             k3tab[i],
             k4tab[i],
             k5tab[i],
             k6tab[i],
             k7tab[i],
             k8tab[i],
             k9tab[i],
             k10tab[i],
             k11tab[i],
             k12tab[i],
             k13tab[i],
             k14tab[i],
             k15tab[i],
             k16tab[i],
             k17tab[i],
             k18tab[i],
             k19tab[i],
             k20tab[i],
             k21tab[i],
             k22tab[i],
             k23tab[i],
             k24tab[i],
             d1tab[i],
             d2tab[i],
             d3tab[i],
             d4tab[i],
             d5tab[i],
             d6tab[i]);
    }

  fclose(FdRateCoeff);
}
/* debug */
#endif /* BG_MOL_NETWORK_TABULATED_RATES */

#endif /* BG_MOL_NETWORK */
