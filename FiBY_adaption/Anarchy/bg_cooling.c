/*
 * Updated to match the procedure outlined in Wiersma, Schaye, and
 * Smith 2009  In addition, new cooling tables were prepared
 * (previous tables inherited a bug from Cloudy that incorrectly
 * calculated the iron cooling)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#if defined(BG_COOLING) && !defined(BG_COOLING_OLD_TABLES)

#include "allvars.h"
#include "proto.h"
#include "bg_cooling.h"
#include "bg_proto.h"
#include "bg_vars.h"


float ****cooling_MetalsNetHeating;
float ****cooling_HplusHeNetHeating;
float ****cooling_HplusHeElectronAbundance;
float ****cooling_ThermalToTemp;

float ***cooling_SolarElectronAbundance; /* ***cooling_CollisionalElectronAbundance; */

#ifdef BG_COOLING_SHIELDING
float ***cooling_MetalsNetHeatingCollisional;
float ***cooling_HplusHeNetHeatingCollisional;
float ***cooling_HplusHeElectronAbundanceCollisional;
float ***cooling_ThermalToTempCollisional;

float **cooling_SolarElectronAbundanceCollisional;
#endif

float *cooling_Temp;
float *cooling_nH;
float *cooling_Therm;
float *cooling_HeFrac;
float *cooling_Redshifts;
float *cooling_SolarAbundances;

char **cooling_HeNames;
char **cooling_ElementNames;
char **cooling_SolarAbundanceNames;

int cooling_N_He; //nHes;
int cooling_N_Temp; //ntemps;
int cooling_N_nH; //nhdens;
int cooling_N_Elements; //nCoolHeats
int cooling_N_SolarAbundances; //nSolar_Abundances
int cooling_N_Redshifts; //nZs

double *cooling_ElementAbundance_SOLAR;
double *cooling_ElementAbundance_SOLARM1;

#ifdef BG_MOL_COOLING
int mol_cooling_H2_N_Temp;
float *mol_cooling_H2_Temp;
float *mol_cooling_H2_Rate;

int mol_cooling_HD_N_Temp;
float *mol_cooling_HD_Temp;
float *mol_cooling_HD_Rate;
#endif

int *ElementNamePointers;
int *SolarAbundanceNamePointers;


double Hefrac;

#ifdef BG_COOLING_SHIELDING
double get_shielding_factor(double n_h);

double get_shielding_factor(double n_h)
{
  return 1 - All.UVBG_PhysDensThresh_HpCM3 * All.UVBG_PhysDensThresh_HpCM3 / (n_h * n_h);
}
#endif

/*
 * ----------------------------------------------------------------------
 * This routine returns the temperature in Kelvin for a given internal
 * energy and a mean molecular weight based on the Helium fraction.
 * Temperature is computed by interpolating on the u_to_temp table.
 * ----------------------------------------------------------------------
 */

double convert_u_to_temp(float d_z, double u, double n_h, double inHe)
{
  int u_i, He_i, n_h_i;

  float d_n_h, d_He, d_u, logT, T;

#ifdef BG_COOLING_SHIELDING
  double shielding_factor;
#endif

  get_index_1d(cooling_HeFrac, cooling_N_He, inHe, &He_i, &d_He);
  get_index_1d(cooling_nH, cooling_N_nH, log10(n_h), &n_h_i, &d_n_h);
  get_index_1d(cooling_Therm, cooling_N_Temp, log10(u), &u_i, &d_u);

#ifdef BG_COOLING_SHIELDING
  /* 1 for full shielding, 0 for no shielding */
  shielding_factor = get_shielding_factor(n_h);
#endif

#ifdef BG_COOLING_SHIELDING
  if(n_h <= All.UVBG_PhysDensThresh_HpCM3)
#endif
    logT = interpol_4d(cooling_ThermalToTemp,
		       0, 1, He_i, He_i + 1,
		       n_h_i, n_h_i + 1, u_i, u_i + 1,
		       d_z, d_He, d_n_h, d_u);
#ifdef BG_COOLING_SHIELDING
  else
    /* in case of self-shielding returns a weighted temperature */
    logT = shielding_factor * interpol_3d(cooling_ThermalToTempCollisional,
					  He_i, n_h_i, u_i, d_He, d_n_h, d_u) +
      (1 - shielding_factor) * interpol_4d(cooling_ThermalToTemp,
					   0, 1, He_i, He_i + 1,
					   n_h_i, n_h_i + 1, u_i, u_i + 1,
					   d_z, d_He, d_n_h, d_u);
  /*
  printf("[DEBUG - convert_u_to_temp] shielding_factor = %g (%g, %g)\n", shielding_factor, All.UVBG_PhysDensThresh_HpCM3, n_h);

  printf("[DEBUG - convert_u_to_temp] temperature = (%g, %g)\n",
	 interpol_3d(cooling_ThermalToTempCollisional,
		     He_i, n_h_i, u_i, d_He, d_n_h, d_u),
	 interpol_4d(cooling_ThermalToTemp,
		     0, 1, He_i, He_i + 1,
		     n_h_i, n_h_i + 1, u_i, u_i + 1,
		     d_z, d_He, d_n_h, d_u));
  */
#endif

  T = pow(10.0, logT);

  if(u_i == 0 && d_u == 0)
    {
      T *= u / pow(10.0, cooling_Therm[0]);

      return T;
    }

  return T;
}


/*
 * ----------------------------------------------------------------------
 * This routine calculates (heating rate-cooling rate)/n_h^2 in cgs units
 * Heating - Cooling is computed by adding up all the contributions
 * from the different elements.
 * ----------------------------------------------------------------------
 */

double CoolingRate(double u, float n_h, double z, float d_z, double particle_Z[], double n_H2, double n_HD)
{
  double LambdaNet = 0.0, electron_abundance, solar_electron_abundance;

  int i;

  double temp;

  int n_h_i, He_i, temp_i;

  float d_n_h, d_He, d_temp;

#ifdef BG_MOL_COOLING
  double dlogt;
#endif

#ifdef BG_COOLING_SHIELDING
  double electron_abundance_collisional;
  double solar_electron_abundance_collisional;
  double shielding_factor;
#endif

  temp = convert_u_to_temp(d_z, u, n_h, Hefrac);

  get_index_1d(cooling_Temp, cooling_N_Temp, log10(temp), &temp_i, &d_temp);
  get_index_1d(cooling_HeFrac, cooling_N_He, Hefrac, &He_i, &d_He);
  get_index_1d(cooling_nH, cooling_N_nH, log10(n_h), &n_h_i, &d_n_h);

#ifdef BG_COOLING_SHIELDING
  /* 1 for full shielding, 0 for no shielding */
  shielding_factor = get_shielding_factor(n_h);
#endif

  /* ------------------ */
  /* Metal-free cooling */
  /* ------------------ */

#ifdef BG_COOLING_SHIELDING
  if(n_h <= All.UVBG_PhysDensThresh_HpCM3)
    {
#endif
      LambdaNet = interpol_4d(cooling_HplusHeNetHeating,
			      0, 1, He_i, He_i + 1,
			      n_h_i, n_h_i + 1, temp_i, temp_i + 1,
			      d_z, d_He, d_n_h, d_temp);

      electron_abundance = interpol_4d(cooling_HplusHeElectronAbundance,
                                       0, 1, He_i, He_i + 1,
                                       n_h_i, n_h_i + 1, temp_i, temp_i + 1,
                                       d_z, d_He, d_n_h, d_temp);
#ifdef BG_COOLING_SHIELDING
    }
  else
    {
      /* weighted cooling rate in case of self-shielding */
      LambdaNet = shielding_factor * interpol_3d(cooling_HplusHeNetHeatingCollisional,
						 He_i, n_h_i, temp_i, d_He, d_n_h, d_temp) +
	(1 - shielding_factor) * interpol_4d(cooling_HplusHeNetHeating, 0, 1,
					     He_i, He_i + 1, n_h_i, n_h_i + 1,
					     temp_i, temp_i + 1, d_z, d_He, d_n_h, d_temp);
      /*
      printf("[DEBUG - CoolingRate] shielding_factor = %g (%g, %g)\n", shielding_factor, All.UVBG_PhysDensThresh_HpCM3, n_h);
      printf("[DEBUG - CoolingRate]      temperature = %g\n", temp);

      printf("[DEBUG - CoolingRate]     cooling rate = (%g, %g)\n",
	     interpol_3d(cooling_HplusHeNetHeatingCollisional,
			 He_i, n_h_i, temp_i, d_He, d_n_h, d_temp),
	     interpol_4d(cooling_HplusHeNetHeating, 0, 1,
			 He_i, He_i + 1, n_h_i, n_h_i + 1,
			 temp_i, temp_i + 1, d_z, d_He, d_n_h, d_temp));
      */

      /* electron abundance for self-shielded gas */
      electron_abundance_collisional = interpol_3d(cooling_HplusHeElectronAbundanceCollisional,
						   He_i, n_h_i, temp_i, d_He, d_n_h, d_temp);

      /* electron abundance for non-shielded gas */
      electron_abundance = interpol_4d(cooling_HplusHeElectronAbundance,
				       0, 1, He_i, He_i + 1,
				       n_h_i, n_h_i + 1, temp_i, temp_i + 1,
				       d_z, d_He, d_n_h, d_temp);
    }
#endif

  /* ------------------ */
  /* Compton cooling    */
  /* ------------------ */

  if(z > cooling_Redshifts[cooling_N_Redshifts - 1] || z > All.REION_H_ZCenter)
    {
      /* inverse Compton cooling is not in collisional table
	 before reionisation so add now */
      LambdaNet -= COMPTON_RATE * (temp - T_CMB0 * (1 + z)) *
	pow((1 + z), 4) * electron_abundance / n_h;
    }
#ifdef BG_COOLING_SHIELDING
  else if(n_h > All.UVBG_PhysDensThresh_HpCM3)
    {
      /* inverse Compton cooling is not in collisional table
         needs to be included for self_shielded gas after re-ionization */
      LambdaNet -= shielding_factor * COMPTON_RATE * (temp - T_CMB0 * (1 + z)) *
        pow((1 + z), 4) * electron_abundance_collisional / n_h;
    }
#endif

  /* ------------- */
  /* Metal cooling */
  /* ------------- */

  if(All.MetDepCoolingOn == 1)
    {
      /* for each element, find the abundance and multiply it
	 by the interpolated heating-cooling */

      //printf("[DEBUG - CoolingRate] do metal cooling\n");

      set_cooling_SolarAbundances(particle_Z, cooling_ElementAbundance_SOLAR);

#ifdef BG_COOLING_SHIELDING
      if(n_h <= All.UVBG_PhysDensThresh_HpCM3)
	{
#endif
	  solar_electron_abundance = interpol_3d(cooling_SolarElectronAbundance,
						 0, n_h_i, temp_i, d_z, d_n_h, d_temp);	/* ne/n_h */

	  for(i = 0; i < cooling_N_Elements; i++)
	    LambdaNet += interpol_4d(cooling_MetalsNetHeating,
				     0, 1, i, i, n_h_i, n_h_i + 1,
				     temp_i, temp_i + 1,
				     d_z, 0.0, d_n_h, d_temp) *
	      (electron_abundance / solar_electron_abundance) *
	      cooling_ElementAbundance_SOLAR[i];
#ifdef BG_COOLING_SHIELDING
	}
      else
	{
	  /* solar electron abundance for self-shielded gas (is it necessary?) */
	  solar_electron_abundance_collisional = interpol_2d(cooling_SolarElectronAbundanceCollisional,
							     n_h_i, temp_i, d_n_h, d_temp);

	  /* solar electron abundance for non-shielded gas */
	  solar_electron_abundance = interpol_3d(cooling_SolarElectronAbundance,
						 0, n_h_i, temp_i, d_z, d_n_h, d_temp); /* ne/n_h */

	  /* weighted cooling rate in case of self-shielding */
	  for(i = 0; i < cooling_N_Elements; i++)
	    {
	      LambdaNet += shielding_factor * interpol_3d_metals(cooling_MetalsNetHeatingCollisional,
								 i, n_h_i, temp_i, d_n_h, d_temp) *
		(electron_abundance_collisional / solar_electron_abundance_collisional) * cooling_ElementAbundance_SOLAR[i] +
		(1 - shielding_factor) * interpol_4d(cooling_MetalsNetHeating,
						     0, 1, i, i, n_h_i, n_h_i + 1, temp_i, temp_i + 1,
						     d_z, 0.0, d_n_h, d_temp) *
		(electron_abundance / solar_electron_abundance) * cooling_ElementAbundance_SOLAR[i];
	      /*
	      printf("[DEBUG - CoolingRate] z cooling rate[%d] = (%g, %g), %g\n", i,
		     interpol_3d_metals(cooling_MetalsNetHeatingCollisional,
					i, n_h_i, temp_i, d_n_h, d_temp),
		     interpol_4d(cooling_MetalsNetHeating,
				 0, 1, i, i, n_h_i, n_h_i + 1, temp_i, temp_i + 1,
				 d_z, 0.0, d_n_h, d_temp),
		     cooling_ElementAbundance_SOLAR[i]);

	      printf("[DEBUG - CoolingRate]        e abundance = (%g, %g)\n", electron_abundance_collisional, electron_abundance);
	      printf("[DEBUG - CoolingRate]  solar e abundance = (%g, %g)\n", solar_electron_abundance_collisional, solar_electron_abundance);
	      */
	    }
	}
#endif
    }


#ifdef BG_MOL_COOLING
  /* ---------- */
  /* H2 cooling */
  /* ---------- */

  if(log10(temp) < mol_cooling_H2_Temp[mol_cooling_H2_N_Temp - 1])
    {
      dlogt = mol_cooling_H2_Temp[1] - mol_cooling_H2_Temp[0];
      temp_i = floor((log10(temp) - mol_cooling_H2_Temp[0]) / dlogt);

      if(temp_i == mol_cooling_H2_N_Temp - 1)
	{
	  temp_i--;
	  d_temp = 1;
	}
      else
	d_temp = (log10(temp) - mol_cooling_H2_Temp[temp_i]) / dlogt;

      LambdaNet -= (mol_cooling_H2_Rate[temp_i] + d_temp * (mol_cooling_H2_Rate[temp_i + 1] - mol_cooling_H2_Rate[temp_i])) * n_H2 / n_h;
    }

  /* ---------- */
  /* HD cooling */
  /* ---------- */

  if(log10(temp) < mol_cooling_HD_Temp[mol_cooling_HD_N_Temp - 1])
    {
      dlogt = mol_cooling_HD_Temp[1] - mol_cooling_HD_Temp[0];
      temp_i = floor((log10(temp) - mol_cooling_HD_Temp[0]) / dlogt);

      if(temp_i == mol_cooling_HD_N_Temp - 1)
	{
	  temp_i--;
	  d_temp = 1;
	}
      else
	d_temp = (log10(temp) - mol_cooling_HD_Temp[temp_i]) / dlogt;

      LambdaNet -= (mol_cooling_HD_Rate[temp_i] + d_temp * (mol_cooling_HD_Rate[temp_i + 1] - mol_cooling_HD_Rate[temp_i])) * n_HD / n_h;
    }

  /* ------- */
  /* Heating */
  /* ------- */

  /* maybe, in the future */
#endif

  return LambdaNet;
}


/*
 * ----------------------------------------------------------------------
 * This routine does the cooling for an active gas particle. 
 * ----------------------------------------------------------------------
 */

#ifdef BG_MOL_COOLING
double DoCooling(double u_old, double rho, double dt, double dz,
		 double inz, float d_z, double particle_Z[], double x_H2, double x_HD)
#else
double DoCooling(double u_old, double rho, double dt, double dz,
		 double inz, float d_z, double particle_Z[])
#endif
{
  float inn_h, XH;
 
  double n_H2, n_HD;

  double du, ratefact, u, u_upper, u_lower, LambdaNet, LambdaTune = 0;

  int i;

  static int HydrogenIndex, HeliumIndex, first_call = 0;


  if(first_call == 0)
    {
      HydrogenIndex = element_index("Hydrogen");
      HeliumIndex = element_index("Helium");
      first_call = 1;
    }

  u = u_old;
  u_lower = u;
  u_upper = u;

#ifdef BG_MOL_COOLING
  XH = particle_Z[HydrogenIndex] - x_H2;
#else
  XH = particle_Z[HydrogenIndex];
#endif
  Hefrac = particle_Z[HeliumIndex] / (XH + particle_Z[HeliumIndex]);

  /* convert Hydrogen mass fraction in Hydrogen number density */
  inn_h = rho * XH / PROTONMASS;

#ifdef BG_MOL_COOLING
  n_H2 = rho * x_H2 / (2 * PROTONMASS);
  n_HD = rho * x_HD / (3 * PROTONMASS);
#else
  n_H2 = n_HD = 0;
#endif

  ratefact = inn_h * inn_h / rho;

  /* set helium and hydrogen reheating term */
  LambdaTune = Reionization(inz, dz, inn_h);

  /* iterative, implicit cooling */
  if(dt > 0)
    LambdaNet = (LambdaTune / (dt * ratefact)) + CoolingRate(u, inn_h, inz, d_z, particle_Z, n_H2, n_HD);
  else
    LambdaNet = CoolingRate(u, inn_h, inz, d_z, particle_Z, n_H2, n_HD);


  if(fabs(ratefact * LambdaNet * dt) < 0.05 * u_old)
    {
      /* cooling rate is small, take the explicit solution */
      u = u_old + ratefact * LambdaNet * dt;

      return u;
    }


  i = 0;

  /* bracketing  */
  if(u - u_old - ratefact * LambdaNet * dt < 0)	/* heating  */
    {
      u_upper *= sqrt(1.1);
      u_lower /= sqrt(1.1);

      while(u_upper - u_old - LambdaTune - ratefact * CoolingRate(u_upper, inn_h, inz, d_z, particle_Z, n_H2, n_HD) * dt <
	    0 && i < MAXITER)
	{
	  u_upper *= 1.1;
	  u_lower *= 1.1;
	  i++;
	}
    }

  if(i == MAXITER)
    printf("Problem with cooling finding upper bound\n");

  if(u - u_old - ratefact * LambdaNet * dt > 0)	/* cooling */
    {
      u_lower /= sqrt(1.1);
      u_upper *= sqrt(1.1);

      i = 0;

      while(u_lower - u_old - LambdaTune - ratefact * CoolingRate(u_lower, inn_h, inz, d_z, particle_Z, n_H2, n_HD) * dt >
	    0 && i < MAXITER)
	{
	  u_upper /= 1.1;
	  u_lower /= 1.1;
	  i++;
	}
    }

  if(i == MAXITER)
    printf("Problem with cooling finding lower bound\n");

  i = 0;

  do				/* iterate to convergence */
    {
      u = 0.5 * (u_lower + u_upper);

      LambdaNet = (LambdaTune / (dt * ratefact)) + CoolingRate(u, inn_h, inz, d_z, particle_Z, n_H2, n_HD);

      if(u - u_old - ratefact * LambdaNet * dt > 0)
	u_upper = u;
      else
	u_lower = u;

      du = u_upper - u_lower;

      i++;

      if(i >= (MAXITER - 10))
	printf("u = %g\n", u);
    }
  while(fabs(du / u) > 1.0e-6 && i < MAXITER);

  if(i >= MAXITER)
    printf("failed to converge in DoCooling()\n");

  return u;
}


double Reionization(double z, double dz, float n_h)
{
  double extra_heating = 0.0;

#ifdef BG_COOLING_SHIELDING
  double shielding_factor;
#endif

#ifdef BG_COOLING_SHIELDING
  shielding_factor = get_shielding_factor(n_h);
#endif

#ifdef BG_COOLING_SHIELDING
  if(n_h <= All.UVBG_PhysDensThresh_HpCM3)
    {
#endif
      /* Hydrogen reionization */
      if(z <= All.REION_H_ZCenter && z - dz > All.REION_H_ZCenter)
	extra_heating += All.REION_H_Heating_ERGpG *
	  (erf((0) / (pow(2., 0.5) * All.REION_H_ZSigma)) -
	   erf((z - All.REION_H_ZCenter) / (pow(2., 0.5) * All.REION_H_ZSigma)));
      else if(z <= All.REION_H_ZCenter)
	extra_heating += All.REION_H_Heating_ERGpG *
	  (erf((z - dz - All.REION_H_ZCenter) / (pow(2., 0.5) * All.REION_H_ZSigma)) -
	   erf((z - All.REION_H_ZCenter) / (pow(2., 0.5) * All.REION_H_ZSigma)));

      /* Helium reionization */
      extra_heating += All.REION_He_Heating_ERGpG *
	(erf((z - dz - All.REION_He_ZCenter) / (pow(2., 0.5) * All.REION_He_ZSigma)) -
	 erf((z - All.REION_He_ZCenter) / (pow(2., 0.5) * All.REION_He_ZSigma))) / 2.0;
#ifdef BG_COOLING_SHIELDING
    }
  else
    {
      //printf("[DEBUG] shielding_factor = %g (%g, %g)\n", shielding_factor, All.UVBG_PhysDensThresh_HpCM3, n_h);

      /* Hydrogen reionization */
      if(z <= All.REION_H_ZCenter && z - dz > All.REION_H_ZCenter)
	extra_heating += (1 - shielding_factor) * All.REION_H_Heating_ERGpG *
          (erf((0) / (pow(2., 0.5) * All.REION_H_ZSigma)) -
           erf((z - All.REION_H_ZCenter) / (pow(2., 0.5) * All.REION_H_ZSigma)));
      else if(z <= All.REION_H_ZCenter)
        extra_heating += (1 - shielding_factor) * All.REION_H_Heating_ERGpG *
          (erf((z - dz - All.REION_H_ZCenter) / (pow(2., 0.5) * All.REION_H_ZSigma)) -
           erf((z - All.REION_H_ZCenter) / (pow(2., 0.5) * All.REION_H_ZSigma)));

      /* Helium reionization */
      extra_heating += (1 - shielding_factor) * All.REION_He_Heating_ERGpG *
	(erf((z - dz - All.REION_He_ZCenter) / (pow(2., 0.5) * All.REION_He_ZSigma)) -
         erf((z - All.REION_He_ZCenter) / (pow(2., 0.5) * All.REION_He_ZSigma))) / 2.0;
    }
#endif

  return extra_heating;
}


/*
 * ----------------------------------------------------------------------
 * Given abundance (solar units) of elements carried by SPH particle,
 * set abundance (solar) of all metals in cooling table index of given
 * element in CoolHeat table in SPH particle is CoolHeat_element_index
 * ----------------------------------------------------------------------
 */

int set_cooling_SolarAbundances(double *element_abundance, double *cooling_element_abundance)
{
  int i, index;

  int static Silicon_SPH_Index = -1;
  int static Calcium_SPH_Index = -1;
  int static Sulphur_SPH_Index = -1;

  int static Silicon_CoolHeat_Index = -1;
  int static Calcium_CoolHeat_Index = -1;
  int static Sulphur_CoolHeat_Index = -1;

  int static first_call = 0;


  if(first_call == 0)
    {
      /* determine (inverse of) solar abundance of these elements */
      for(i = 0; i < cooling_N_Elements; i++)
	{
	  index = get_element_index(cooling_SolarAbundanceNames,
				    cooling_N_SolarAbundances,
				    cooling_ElementNames[i]);

	  if(index < 0)
	    endrun(-12345);

	  index = SolarAbundanceNamePointers[i];

	  cooling_ElementAbundance_SOLARM1[i] = 1. / cooling_SolarAbundances[index];

	  index = ElementNamePointers[i];

	  if(index < 0 && ThisTask == 0)
	    printf("[bg_cooling] element not found %s\n", cooling_ElementNames[i]);
	}

      /* Sulphur tracks Silicon: may choose not to follow Sulphur as SPH element */
      /* Same is true for Calcium */
      /* We will assume the code tracks Silicon, and may need to scale Calcium and Sulphur accordingly */

      Silicon_SPH_Index = element_index("Silicon");
      Calcium_SPH_Index = element_index("Calcium");
      Sulphur_SPH_Index = element_index("Sulphur");

      Silicon_CoolHeat_Index =
	get_element_index(cooling_ElementNames, cooling_N_Elements, "Silicon");
      Calcium_CoolHeat_Index =
	get_element_index(cooling_ElementNames, cooling_N_Elements, "Calcium");
      Sulphur_CoolHeat_Index =
	get_element_index(cooling_ElementNames, cooling_N_Elements, "Sulphur");

      if(Silicon_CoolHeat_Index == -1 ||
	 Calcium_CoolHeat_Index == -1 ||
	 Sulphur_CoolHeat_Index == -1)
	{
	  if(ThisTask == 0)
	    printf("[bg_cooling] error: did not find Si or Ca or S??\n");
	  endrun(-1233);
	}

      first_call = 1;
    }

  for(i = 0; i < cooling_N_Elements; i++)
    {
      if(i == Calcium_CoolHeat_Index && Calcium_SPH_Index == -1)
	/* SPH does not track Calcium: use Si abundance */
	if(Silicon_SPH_Index == -1)
	  cooling_element_abundance[i] = 0.0;
	else
	  cooling_element_abundance[i] = element_abundance[Silicon_SPH_Index] *
	    cooling_ElementAbundance_SOLARM1[Silicon_CoolHeat_Index];
      else if(i == Sulphur_CoolHeat_Index && Sulphur_SPH_Index == -1)
	/* SPH does not track Sulphur: use Si abundance */
	if(Silicon_SPH_Index == -1)
	  cooling_element_abundance[i] = 0.0;
	else
	  cooling_element_abundance[i] = element_abundance[Silicon_SPH_Index] *
	    cooling_ElementAbundance_SOLARM1[Silicon_CoolHeat_Index];
      else
	cooling_element_abundance[i] = element_abundance[ElementNamePointers[i]] *
	  cooling_ElementAbundance_SOLARM1[i];
    }

  return 0;
}

#endif
