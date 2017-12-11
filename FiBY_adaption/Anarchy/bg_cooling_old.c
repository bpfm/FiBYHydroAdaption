/*Implementing version */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#if defined(BG_COOLING) && defined(BG_COOLING_OLD_TABLES)

#include "allvars.h"
#include "proto.h"
#include "bg_cooling_old.h"
#include "bg_proto.h"
#include "bg_vars.h"


float ****cooling_MetalsNetHeating;
float ****cooling_HplusHeNetHeating;
float ****cooling_ThermalToTemp;

float ***cooling_CollisionalElectronAbundance;

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

int *ElementNamePointers;
int *SolarAbundanceNamePointers;


float Hefrac;


/*
 * ----------------------------------------------------------------------
 * This routine returns the temperature in Kelvin for a given internal
 * energy and a mean molecular weight based on the Helium fraction.
 * Temperature is computed by interpolating on the u_to_temp table.
 * ----------------------------------------------------------------------
 */

double convert_u_to_temp(float d_z, double u, double inn_h, double inHe)
{
  int u_i, He_i, n_h_i;

  float d_n_h, d_He, d_u, logT, T;

  get_index_1d(cooling_HeFrac, cooling_N_He, inHe, &He_i, &d_He);
  get_index_1d(cooling_nH, cooling_N_nH, log10(inn_h), &n_h_i, &d_n_h);
  get_index_1d(cooling_Therm, cooling_N_Temp, log10(u), &u_i, &d_u);

  logT = interpol_4d(cooling_ThermalToTemp,
		     0, 1, He_i, He_i + 1,
		     n_h_i, n_h_i + 1, u_i, u_i + 1,
		     d_z, d_He, d_n_h, d_u);

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

double CoolingRate(double u, float n_h, double z, float d_z, double particle_Z[])
{
  double T_gam, LambdaNet = 0.0, electron_abundance;

  int i;

  double temp;

  int n_h_i, He_i, temp_i;

  float d_n_h, d_He, d_temp;


  temp = convert_u_to_temp(d_z, u, n_h, Hefrac);

  get_index_1d(cooling_Temp, cooling_N_Temp, log10(temp), &temp_i, &d_temp);
  get_index_1d(cooling_HeFrac, cooling_N_He, Hefrac, &He_i, &d_He);
  get_index_1d(cooling_nH, cooling_N_nH, log10(n_h), &n_h_i, &d_n_h);

  /* ------------------ */
  /* Metal-free cooling */
  /* ------------------ */

  LambdaNet = interpol_4d(cooling_HplusHeNetHeating,
			  0, 1, He_i, He_i + 1,
			  n_h_i, n_h_i + 1, temp_i, temp_i + 1,
			  d_z, d_He, d_n_h, d_temp);

  /* ------------------ */
  /* Compton cooling    */
  /* ------------------ */

  if(z > cooling_Redshifts[cooling_N_Redshifts - 1] || z > All.REION_H_ZCenter)
    {
      /* inverse Compton cooling is not in collisional table before reionisation */
      /* so add now */
      electron_abundance = interpol_3d(cooling_CollisionalElectronAbundance,
				       He_i, n_h_i, temp_i, d_He, d_n_h, d_temp);	/* ne/n_h */

      T_gam = T_CMB0 * (1 + z);

      LambdaNet -= COMPTON_RATE * (temp - T_CMB0 * (1 + z)) *
	pow((1 + z), 4) * electron_abundance / n_h;
    }

  /* ------------- */
  /* Metal cooling */
  /* ------------- */

  if(All.MetDepCoolingOn == 1)
    {
      /* for each element, find the abundance and multiply it by the interpolated heating-cooling */
      set_cooling_SolarAbundances(particle_Z, cooling_ElementAbundance_SOLAR);

      for(i = 0; i < cooling_N_Elements; i++)
	LambdaNet += interpol_4d(cooling_MetalsNetHeating,
				 0, 1, i, i, n_h_i, n_h_i + 1,
				 temp_i, temp_i + 1,
				 d_z, 0.0, d_n_h, d_temp) *
	  cooling_ElementAbundance_SOLAR[i];
    }

  return LambdaNet;
}


/*
 * ----------------------------------------------------------------------
 * This routine does the cooling for an active gas particle. 
 * ----------------------------------------------------------------------
 */

double DoCooling(double u_old, double rho, double dt, double dz,
		 double inz, float d_z, double particle_Z[])
{
  float inn_h, XH;

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

  XH = particle_Z[HydrogenIndex];
  /* cooling interpolation bug fixed here
  Hefrac = particle_Z[HeliumIndex];
  */
  Hefrac = particle_Z[HeliumIndex] / (XH + particle_Z[HeliumIndex]);

  /* convert Hydrogen mass fraction in Hydrogen number density */
  inn_h = rho * XH / PROTONMASS;

  ratefact = inn_h * inn_h / rho;

  /* set helium and hydrogen reheating term */
  LambdaTune = Reionization(inz, dz);

  /* iterative, implicit cooling */
  if(dt > 0)
    LambdaNet = (LambdaTune / (dt * ratefact)) + CoolingRate(u, inn_h, inz, d_z, particle_Z);
  else
    LambdaNet = CoolingRate(u, inn_h, inz, d_z, particle_Z);


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

      while(u_upper - u_old - LambdaTune -ratefact * CoolingRate(u_upper, inn_h, inz, d_z, particle_Z) * dt <
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

      while(u_lower - u_old - LambdaTune - ratefact * CoolingRate(u_lower, inn_h, inz, d_z, particle_Z) * dt >
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

      LambdaNet = (LambdaTune / (dt * ratefact)) + CoolingRate(u, inn_h, inz, d_z, particle_Z);

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


double Reionization(double z, double dz)
{
  double extra_heating = 0.0;

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
    (erf((z - All.REION_He_ZCenter) / (pow(2., 0.5) * All.REION_He_ZSigma)) -
     erf((z + dz - All.REION_He_ZCenter) / (pow(2., 0.5) * All.REION_He_ZSigma))) / 2.0;

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
