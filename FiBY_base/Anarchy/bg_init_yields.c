#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <hdf5.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#ifdef BG_STELLAR_EVOLUTION

#include "allvars.h"
#include "proto.h"
#include "bg_yields.h"
#include "bg_proto.h"
#include "bg_vars.h"


/* #define TEST_YIELDS 1 */

void TestYields(void);


void init_yields(void)
{
  MyFloat lm_min, lm_max, dlm;

  int i, j, k;

  double total_mass;

  /* ---------------------------------------- */
  /* master processor gets tables dimensions, */
  /* allocate tables arrays and read them in  */
  /* ---------------------------------------- */

  if(ThisTask == 0)
    read_yield_tables();

  /* ------------------------------------------------- */
  /* broadcast the table dimensions packed in a buffer */
  /* ------------------------------------------------- */

  bcast_yield_table_dim();

  /* ----------------------------------------------- */
  /* allocate yield tables for non-master processors */
  /* ----------------------------------------------- */

  if(ThisTask != 0)
    allocate_yield_tables();

  /* --------------------------------------- */
  /* broadcast the tables packed in a buffer */
  /* --------------------------------------- */

  bcast_yield_tables();

  /* -------------------------------------------- */
  /* arrays for spline interpolating stellar mass */
  /* -------------------------------------------- */

  lm_min = log10(All.IMF_MinMass_MSUN);	/* min mass in solar masses */
  lm_max = log10(All.IMF_MaxMass_MSUN);	/* max mass in solar masses */

  dlm = (lm_max - lm_min) / (N_MASS_BIN - 1);

  yield_mass_bin = (double *) mymalloc(N_MASS_BIN * sizeof(double));

  for(i = 0; i < N_MASS_BIN; i++)
    yield_mass_bin[i] = dlm * i + lm_min;

#ifdef BG_POPIII
  lm_min = log10(All.POPIII_IMF_MinMass_MSUN); /* min mass in solar masses */
  lm_max = log10(All.POPIII_IMF_MaxMass_MSUN); /* max mass in solar masses */

  dlm = (lm_max - lm_min) / (POPIII_N_MASS_BIN - 1);

  popiii_yield_mass_bin = (double *) mymalloc(POPIII_N_MASS_BIN * sizeof(double));

  for(i = 0; i < POPIII_N_MASS_BIN; i++)
    popiii_yield_mass_bin[i] = dlm * i + lm_min;
#endif

  /* ----------------------- */
  /* convert tables to log10 */
  /* ----------------------- */

  /* SNII tables */
  for(k = 0; k < yieldsSNII.N_MASS; k++)
    yieldsSNII.Mass[k] = log10(yieldsSNII.Mass[k]);

  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      if(yieldsSNII.Metallicity[i] > 0)
	yieldsSNII.Metallicity[i] = log10(yieldsSNII.Metallicity[i]);
      else
	yieldsSNII.Metallicity[i] = MIN_METAL;
    }

#ifdef BG_DUST
  for(i = 0; i < dustSNII.N_MASS; i++)
    dustSNII.Mass[i] = log10(dustSNII.Mass[i]);
#endif

  /* AGB tables */
  for(i = 0; i < yieldsAGB.N_MASS; i++)
    yieldsAGB.Mass[i] = log10(yieldsAGB.Mass[i]);

  for(i = 0; i < yieldsAGB.N_Z; i++)
    {
      if(yieldsAGB.Metallicity[i] > 0)
	yieldsAGB.Metallicity[i] = log10(yieldsAGB.Metallicity[i]);
      else
	yieldsAGB.Metallicity[i] = MIN_METAL;
    }

#ifdef BG_DUST
  for(i = 0; i < dustAGB.N_MASS; i++)
    dustAGB.Mass[i] = log10(dustAGB.Mass[i]);
#endif

   /* POPIII tables */
#ifdef BG_POPIII
  for(i = 0; i < yieldsPOPIII.N_MASS; i++)
    yieldsPOPIII.Mass[i] = log10(yieldsPOPIII.Mass[i]);
#ifdef BG_DUST
  for(i = 0; i < dustPOPIII.N_MASS; i++)
    dustPOPIII.Mass[i] = log10(dustPOPIII.Mass[i]);
#endif
#endif

  /* --------------------- */
  /* allocate table arrays */
  /* --------------------- */

  /* SNIa tables */
  yieldsSNIa.SPH = (double *) mymalloc(BG_NELEMENTS * sizeof(double));

  /* SNII tables */
  yieldsSNII.SPH = (double ***) mymalloc(yieldsSNII.N_Z * sizeof(double **));
  yieldsSNII.Ejecta_SPH = (double **) mymalloc(yieldsSNII.N_Z * sizeof(double *));
  yieldsSNII.TotalMetals_SPH = (double **) mymalloc(yieldsSNII.N_Z * sizeof(double *));
  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      yieldsSNII.SPH[i] = (double **) mymalloc(BG_NELEMENTS * sizeof(double *));
      yieldsSNII.Ejecta_SPH[i] = (double *) mymalloc(N_MASS_BIN * sizeof(double));
      yieldsSNII.TotalMetals_SPH[i] = (double *) mymalloc(N_MASS_BIN * sizeof(double));
      for(j = 0; j < BG_NELEMENTS; j++)
	yieldsSNII.SPH[i][j] = (double *) mymalloc(N_MASS_BIN * sizeof(double));
    }

#ifdef BG_DUST
  dustSNII.TotalDust_SPH = (double **) mymalloc(dustSNII.N_Z * sizeof(double *));
  for(i = 0; i < dustSNII.N_Z; i++)
    dustSNII.TotalDust_SPH[i] = (double *) mymalloc(N_MASS_BIN * sizeof(double));
#endif

  /* AGB tables */
  yieldsAGB.SPH = (double ***) mymalloc(yieldsAGB.N_Z * sizeof(double **));
  yieldsAGB.Ejecta_SPH = (double **) mymalloc(yieldsAGB.N_Z * sizeof(double *));
  yieldsAGB.TotalMetals_SPH = (double **) mymalloc(yieldsAGB.N_Z * sizeof(double *));
  for(i = 0; i < yieldsAGB.N_Z; i++)
    {
      yieldsAGB.Ejecta_SPH[i] = (double *) mymalloc(N_MASS_BIN * sizeof(double));
      yieldsAGB.TotalMetals_SPH[i] = (double *) mymalloc(N_MASS_BIN * sizeof(double));
      yieldsAGB.SPH[i] = (double **) mymalloc(BG_NELEMENTS * sizeof(double *));
      for(j = 0; j < BG_NELEMENTS; j++)
	yieldsAGB.SPH[i][j] = (double *) mymalloc(N_MASS_BIN * sizeof(double));
    }

#ifdef BG_DUST
  dustAGB.TotalDust_SPH = (double **) mymalloc(dustAGB.N_Z * sizeof(double *));
  for(i = 0; i < dustAGB.N_Z; i++)
    dustAGB.TotalDust_SPH[i] = (double *) mymalloc(N_MASS_BIN * sizeof(double));
#endif

  /* POPIII tables */
#ifdef BG_POPIII
  yieldsPOPIII.SPH = (double **) mymalloc(BG_NELEMENTS * sizeof(double *));
  yieldsPOPIII.Ejecta_SPH = (double *) mymalloc(POPIII_N_MASS_BIN * sizeof(double));
  yieldsPOPIII.TotalMetals_SPH = (double *) mymalloc(POPIII_N_MASS_BIN * sizeof(double));
  for(i = 0; i < BG_NELEMENTS; i++)
    yieldsPOPIII.SPH[i] = (double *) mymalloc(POPIII_N_MASS_BIN * sizeof(double));
#ifdef BG_DUST
  dustPOPIII.TotalDust_SPH = (double *) mymalloc(POPIII_N_MASS_BIN * sizeof(double));
#endif
#endif

  /* for each element retable the yields */
  for(i = 0; i < BG_NELEMENTS; i++)
    compute_yields(i);

  /* retable the ejecta and total metal mass released */
  compute_ejecta();

  /* double check scaling */
  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      for(k = 0; k < N_MASS_BIN; k++)
	{

	  total_mass =
	    yieldsSNII.TotalMetals_SPH[i][k] + yieldsSNII.Metallicity[i] * yieldsSNII.Ejecta_SPH[i][k];
	  if(element_index("Hydrogen") > 0)
	    total_mass +=
	      yieldsSNII.SPH[i][element_index("Hydrogen")][k] + (1.0 -
								 yieldsSNII.Metallicity[i]) * 0.75 *
	      yieldsSNII.Ejecta_SPH[i][k];
	  if(element_index("Helium") > 0)
	    total_mass +=
	      yieldsSNII.SPH[i][element_index("Helium")][k] + (1.0 -
							       yieldsSNII.Metallicity[i]) * 0.25 *
	      yieldsSNII.Ejecta_SPH[i][k];


	  if(total_mass > pow(10.0, yield_mass_bin[k]))
	    {
	      printf("Not conserving mass! Your type II scale factors are probably too high");
	      endrun(-2);
	    }
	}
    }

#if TEST_YIELDS == 1
  if(ThisTask == 0)
    printf("Testing yields...");

  TestYields();

  MPI_Barrier(MPI_COMM_WORLD);
  endrun(-10);
#endif
}


void compute_yields(int SPH_element_index)
{
  char element_name[EL_NAME_LENGTH];

  int i, k;

  double *yield, result, typeII_factor = 1;

  int element_index;

  gsl_interp_accel *my_accel_ptr;

  gsl_spline *my_spline_ptr;

  double Mtot, Ntot;

  static int first_call = 0;

#ifdef TEST_YIELDS
  FILE *outfile;

  char filename[100];
#endif

  if(first_call == 0)
    {
      /* first test whether imf is normalised */
      Ntot = integrate_imf(imf_mass_bin_log10[0], imf_mass_bin_log10[N_MASS_BIN - 1], 0.0, 0);
      Mtot = integrate_imf(imf_mass_bin_log10[0], imf_mass_bin_log10[N_MASS_BIN - 1], 0.0, 1);
      if(fabs(Mtot - 1.0) > 1.0e-2)
	{
	  printf("Error: IMF not properly initialized \n");
	  endrun(911);
	}
      first_call = 1;
    }

  /*  if(SPH_element_index < BG_NELEMENTS) */
  sprintf(element_name, ElementNames[SPH_element_index]);	/* metals, including Helium */
  /*  else */
  /*    strcpy(element_name, "Hydrogen"); */

  if(strcmp(element_name, "Hydrogen") == 0)
    typeII_factor = All.SNII_Factor_Hydrogen;
  else if(strcmp(element_name, "Helium") == 0)
    typeII_factor = All.SNII_Factor_Helium;
  else if(strcmp(element_name, "Carbon") == 0)
    typeII_factor = All.SNII_Factor_Carbon;
  else if(strcmp(element_name, "Nitrogen") == 0)
    typeII_factor = All.SNII_Factor_Nitrogen;
  else if(strcmp(element_name, "Oxygen") == 0)
    typeII_factor = All.SNII_Factor_Oxygen;
  else if(strcmp(element_name, "Neon") == 0)
    typeII_factor = All.SNII_Factor_Neon;
  else if(strcmp(element_name, "Magnesium") == 0)
    typeII_factor = All.SNII_Factor_Magnesium;
  else if(strcmp(element_name, "Silicon") == 0)
    typeII_factor = All.SNII_Factor_Silicon;
  else if(strcmp(element_name, "Iron") == 0)
    typeII_factor = All.SNII_Factor_Iron;

  if(ThisTask == 0)
    printf("Computing yield for %s\t index=%d\n", element_name, SPH_element_index);

  /*  SNIa yields */
  element_index = get_element_index(yieldsSNIa.ElementName, yieldsSNIa.N_ELEMENTS, element_name);

  if(element_index < 0)
    {
      if(ThisTask == 0)
	printf("SNIa: element not found %s\n", element_name);

      endrun(-2);
    }

  yieldsSNIa.SPH[SPH_element_index] = yieldsSNIa.Yield[element_index];

  /* SNII yields */
  element_index = get_element_index(yieldsSNII.ElementName, yieldsSNII.N_ELEMENTS, element_name);

  if(element_index < 0)
    {
      if(ThisTask == 0)
	printf("SNII: element not found %s\n", element_name);

      endrun(-2);
    }

  yield = (double *) mymalloc(yieldsSNII.N_MASS * sizeof(double));

  my_accel_ptr = gsl_interp_accel_alloc();
  my_spline_ptr = gsl_spline_alloc(gsl_interp_linear, yieldsSNII.N_MASS);

  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      for(k = 0; k < yieldsSNII.N_MASS; k++)
	yield[k] = yieldsSNII.Yield[i][element_index][k] / pow(10.0, yieldsSNII.Mass[k]);

#ifdef TEST_YIELDS
      sprintf(filename, "SNIIyields.%s.table.%d", element_name, i);
      outfile = fopen(filename, "w");
      for(k = 0; k < yieldsSNII.N_MASS; k++)
	fprintf(outfile, "%g\t%g\n", yieldsSNII.Mass[k], yield[k]);
      fclose(outfile);
#endif

      gsl_spline_init(my_spline_ptr, yieldsSNII.Mass, yield, yieldsSNII.N_MASS);

      for(k = 0; k < N_MASS_BIN; k++)
	{
	  if(yield_mass_bin[k] < yieldsSNII.Mass[0])
	    result = yield[0];
	  else if(yield_mass_bin[k] > yieldsSNII.Mass[yieldsSNII.N_MASS - 1])
	    result = yield[yieldsSNII.N_MASS - 1];
	  else
	    result = gsl_spline_eval(my_spline_ptr, yield_mass_bin[k], my_accel_ptr);

	  yieldsSNII.SPH[i][SPH_element_index][k] = pow(10.0, yield_mass_bin[k]) * result;
	}

#ifdef TEST_YIELDS
      sprintf(filename, "SNIIyields.%s.spline.%d", element_name, i);
      outfile = fopen(filename, "w");
      for(k = 0; k < N_MASS_BIN; k++)
	fprintf(outfile, "%g\t%g\n", yield_mass_bin[k], yieldsSNII.SPH[i][SPH_element_index][k]);
      fclose(outfile);
#endif
    }

  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      for(k = 0; k < N_MASS_BIN; k++)
	{
	  if(strcmp(element_name, "Hydrogen") != 0 || strcmp(element_name, "Helium") != 0)
	    yieldsSNII.TotalMetals_SPH[i][k] += (typeII_factor - 1) * yieldsSNII.SPH[i][SPH_element_index][k];

	  yieldsSNII.SPH[i][SPH_element_index][k] *= typeII_factor;
	}
    }

  myfree(yield);
  gsl_spline_free(my_spline_ptr);
  gsl_interp_accel_free(my_accel_ptr);

  /* AGB yields */
  element_index = get_element_index(yieldsAGB.ElementName, yieldsAGB.N_ELEMENTS, element_name);

  if(element_index < 0)
    {
      if(ThisTask == 0)
	printf("AGB: element not found: %s\n", element_name);

      for(i = 0; i < yieldsAGB.N_Z; i++)
	{
	  for(k = 0; k < N_MASS_BIN; k++)
	    yieldsAGB.SPH[i][SPH_element_index][k] = 0.0;
	}
    }
  else
    {
      yield = (double *) mymalloc(yieldsAGB.N_MASS * sizeof(double));

      my_accel_ptr = gsl_interp_accel_alloc();
      my_spline_ptr = gsl_spline_alloc(gsl_interp_linear, yieldsAGB.N_MASS);

      for(i = 0; i < yieldsAGB.N_Z; i++)
	{
	  for(k = 0; k < yieldsAGB.N_MASS; k++)
	    yield[k] = yieldsAGB.Yield[i][element_index][k] / pow(10.0, yieldsAGB.Mass[k]);

#ifdef TEST_YIELDS
	  sprintf(filename, "AGByields.%s.table.%d", element_name, i);
	  outfile = fopen(filename, "w");
	  for(k = 0; k < yieldsAGB.N_MASS; k++)
	    fprintf(outfile, "%g\t%g\n", yieldsAGB.Mass[k], yield[k]);
	  fclose(outfile);
#endif

	  gsl_spline_init(my_spline_ptr, yieldsAGB.Mass, yield, yieldsAGB.N_MASS);

	  for(k = 0; k < N_MASS_BIN; k++)
	    {
	      if(yield_mass_bin[k] < yieldsAGB.Mass[0])
		result = yield[0];
	      else if(yield_mass_bin[k] > yieldsAGB.Mass[yieldsAGB.N_MASS - 1])
		result = yield[yieldsAGB.N_MASS - 1];
	      else
		result = gsl_spline_eval(my_spline_ptr, yield_mass_bin[k], my_accel_ptr);

	      yieldsAGB.SPH[i][SPH_element_index][k] = pow(10.0, yield_mass_bin[k]) * result;
	    }

#ifdef TEST_YIELDS
	  sprintf(filename, "AGByields.%s.spline.%d", element_name, i);
	  outfile = fopen(filename, "w");
	  for(k = 0; k < N_MASS_BIN; k++)
	    fprintf(outfile, "%g\t%g\n", yield_mass_bin[k], yieldsAGB.SPH[i][SPH_element_index][k]);
	  fclose(outfile);
#endif
	}

      myfree(yield);
      gsl_spline_free(my_spline_ptr);
      gsl_interp_accel_free(my_accel_ptr);
    }

  /* POPIII yields */
#ifdef BG_POPIII
  element_index = get_element_index(yieldsPOPIII.ElementName, yieldsPOPIII.N_ELEMENTS, element_name);

  if(element_index < 0)
    {
      if(ThisTask == 0)
	printf("POPIII: element not found: %s\n", element_name);

      for(k = 0; k < POPIII_N_MASS_BIN; k++)
	yieldsPOPIII.SPH[SPH_element_index][k] = 0.0;
    }
  else
    {
      yield = (double *) mymalloc(yieldsPOPIII.N_MASS * sizeof(double));

      my_accel_ptr = gsl_interp_accel_alloc();
      my_spline_ptr = gsl_spline_alloc(gsl_interp_linear, yieldsPOPIII.N_MASS);

      for(k = 0; k < yieldsPOPIII.N_MASS; k++)
	yield[k] = yieldsPOPIII.Yield[element_index][k] / pow(10.0, yieldsPOPIII.Mass[k]);

#ifdef TEST_YIELDS
      sprintf(filename, "POPIIIyields.%s.table", element_name);
      outfile = fopen(filename, "w");
      for(k = 0; k < yieldsPOPIII.N_MASS; k++)
	fprintf(outfile, "%g\t%g\n", yieldsPOPIII.Mass[k], yield[k]);
      fclose(outfile);
#endif

      gsl_spline_init(my_spline_ptr, yieldsPOPIII.Mass, yield, yieldsPOPIII.N_MASS);

      for(k = 0; k < POPIII_N_MASS_BIN; k++)
	{
	  if(popiii_yield_mass_bin[k] < yieldsPOPIII.Mass[0])
	    result = yield[0];
	  else if(popiii_yield_mass_bin[k] > yieldsPOPIII.Mass[yieldsPOPIII.N_MASS - 1])
	    result = yield[yieldsPOPIII.N_MASS - 1];
	  else
	    result = gsl_spline_eval(my_spline_ptr, popiii_yield_mass_bin[k], my_accel_ptr);

	  yieldsPOPIII.SPH[SPH_element_index][k] = pow(10.0, popiii_yield_mass_bin[k]) * result;

	  /* above 260 solar masses the star collapses into a black hole with no ejecta */
	  if(popiii_yield_mass_bin[k] > yieldsPOPIII.Mass[yieldsPOPIII.N_MASS - 1])
	    yieldsPOPIII.SPH[SPH_element_index][k] = 0;
	}

#ifdef TEST_YIELDS
      sprintf(filename, "POPIIIyields.%s.spline", element_name);
      outfile = fopen(filename, "w");
      for(k = 0; k < POPIII_N_MASS_BIN; k++)
	fprintf(outfile, "%g\t%g\n", popiii_yield_mass_bin[k], yieldsPOPIII.SPH[SPH_element_index][k]);
      fclose(outfile);
#endif

      myfree(yield);
      gsl_spline_free(my_spline_ptr);
      gsl_interp_accel_free(my_accel_ptr);
    }
#endif /* BG_POPIII */

}

void compute_ejecta()
{
  int i, k;

  double *yield, result;

  gsl_interp_accel *my_accel_ptr;

  gsl_spline *my_spline_ptr;

#ifdef TEST_YIELDS
  FILE *outfile;

  char filename[100];
#endif

  if(ThisTask == 0)
    printf("Computing ejecta\n");

  /* SNII yields */
  yield = (double *) mymalloc(yieldsSNII.N_MASS * sizeof(double));

  my_accel_ptr = gsl_interp_accel_alloc();
  my_spline_ptr = gsl_spline_alloc(gsl_interp_linear, yieldsSNII.N_MASS);

  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      for(k = 0; k < yieldsSNII.N_MASS; k++)
	yield[k] = yieldsSNII.Ejecta[i][k] / pow(10.0, yieldsSNII.Mass[k]);

      gsl_spline_init(my_spline_ptr, yieldsSNII.Mass, yield, yieldsSNII.N_MASS);

      for(k = 0; k < N_MASS_BIN; k++)
	{
	  if(yield_mass_bin[k] < yieldsSNII.Mass[0])
	    result = yield[0];
	  else if(yield_mass_bin[k] > yieldsSNII.Mass[yieldsSNII.N_MASS - 1])
	    result = yield[yieldsSNII.N_MASS - 1];
	  else
	    result = gsl_spline_eval(my_spline_ptr, yield_mass_bin[k], my_accel_ptr);

	  yieldsSNII.Ejecta_SPH[i][k] = pow(10.0, yield_mass_bin[k]) * result;
	}
    }

  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      for(k = 0; k < yieldsSNII.N_MASS; k++)
	yield[k] = yieldsSNII.TotalMetals[i][k] / pow(10.0, yieldsSNII.Mass[k]);

      gsl_spline_init(my_spline_ptr, yieldsSNII.Mass, yield, yieldsSNII.N_MASS);

      for(k = 0; k < N_MASS_BIN; k++)
	{
	  if(yield_mass_bin[k] < yieldsSNII.Mass[0])
	    result = yield[0];
	  else if(yield_mass_bin[k] > yieldsSNII.Mass[yieldsSNII.N_MASS - 1])
	    result = yield[yieldsSNII.N_MASS - 1];
	  else
	    result = gsl_spline_eval(my_spline_ptr, yield_mass_bin[k], my_accel_ptr);

	  yieldsSNII.TotalMetals_SPH[i][k] = pow(10.0, yield_mass_bin[k]) * result;
	}
    }

  myfree(yield);
  gsl_spline_free(my_spline_ptr);
  gsl_interp_accel_free(my_accel_ptr);

#ifdef BG_DUST
  yield = (double *) mymalloc(dustSNII.N_MASS * sizeof(double));

  my_accel_ptr = gsl_interp_accel_alloc();
  my_spline_ptr = gsl_spline_alloc(gsl_interp_linear, dustSNII.N_MASS);

  for(i = 0; i < dustSNII.N_Z; i++)
    {
      for(k = 0; k < dustSNII.N_MASS; k++)
        yield[k] = dustSNII.TotalDust[i][k] / pow(10.0, dustSNII.Mass[k]);

#ifdef TEST_YIELDS
      sprintf(filename, "SNIIyields.dust.table.%d", i);
      outfile = fopen(filename, "w");
      for(k = 0; k < dustSNII.N_MASS; k++)
        fprintf(outfile, "%g\t%g\n", dustSNII.Mass[k], yield[k]);
      fclose(outfile);
#endif

      gsl_spline_init(my_spline_ptr, dustSNII.Mass, yield, dustSNII.N_MASS);

      for(k = 0; k < N_MASS_BIN; k++)
        {
          if(yield_mass_bin[k] < dustSNII.Mass[0])
            result = yield[0];
          else if(yield_mass_bin[k] > dustSNII.Mass[dustSNII.N_MASS - 1])
            result = yield[dustSNII.N_MASS - 1];
          else
            result = gsl_spline_eval(my_spline_ptr, yield_mass_bin[k], my_accel_ptr);

          dustSNII.TotalDust_SPH[i][k] = pow(10.0, yield_mass_bin[k]) * result;
	  }

#ifdef TEST_YIELDS
      sprintf(filename, "SNIIyields.dust.spline.%d", i);
      outfile = fopen(filename, "w");
      for(k = 0; k < N_MASS_BIN; k++)
        fprintf(outfile, "%g\t%g\n", yield_mass_bin[k], dustSNII.TotalDust_SPH[i][k]);
      fclose(outfile);
#endif
    }

  gsl_spline_free(my_spline_ptr);
  gsl_interp_accel_free(my_accel_ptr);
  myfree(yield);
#endif

/* AGB yields */
  yield = (double *) mymalloc(yieldsAGB.N_MASS * sizeof(double));

  my_accel_ptr = gsl_interp_accel_alloc();
  my_spline_ptr = gsl_spline_alloc(gsl_interp_linear, yieldsAGB.N_MASS);

  for(i = 0; i < yieldsAGB.N_Z; i++)
    {
      for(k = 0; k < yieldsAGB.N_MASS; k++)
	yield[k] = yieldsAGB.Ejecta[i][k] / pow(10.0, yieldsAGB.Mass[k]);

      gsl_spline_init(my_spline_ptr, yieldsAGB.Mass, yield, yieldsAGB.N_MASS);

      for(k = 0; k < N_MASS_BIN; k++)
	{
	  if(yield_mass_bin[k] < yieldsAGB.Mass[0])
	    result = yield[0];
	  else if(yield_mass_bin[k] > yieldsAGB.Mass[yieldsAGB.N_MASS - 1])
	    result = yield[yieldsAGB.N_MASS - 1];
	  else
	    result = gsl_spline_eval(my_spline_ptr, yield_mass_bin[k], my_accel_ptr);

	  yieldsAGB.Ejecta_SPH[i][k] = pow(10.0, yield_mass_bin[k]) * result;
	}
    }

  for(i = 0; i < yieldsAGB.N_Z; i++)
    {
      for(k = 0; k < yieldsAGB.N_MASS; k++)
	yield[k] = yieldsAGB.TotalMetals[i][k] / pow(10.0, yieldsAGB.Mass[k]);

      gsl_spline_init(my_spline_ptr, yieldsAGB.Mass, yield, yieldsAGB.N_MASS);

      for(k = 0; k < N_MASS_BIN; k++)
	{
	  if(yield_mass_bin[k] < yieldsAGB.Mass[0])
	    result = yield[0];
	  else if(yield_mass_bin[k] > yieldsAGB.Mass[yieldsAGB.N_MASS - 1])
	    result = yield[yieldsAGB.N_MASS - 1];
	  else
	    result = gsl_spline_eval(my_spline_ptr, yield_mass_bin[k], my_accel_ptr);

	  yieldsAGB.TotalMetals_SPH[i][k] = pow(10.0, yield_mass_bin[k]) * result;
	}
    }

  myfree(yield);
  gsl_spline_free(my_spline_ptr);
  gsl_interp_accel_free(my_accel_ptr);

#ifdef BG_DUST
  yield = (double *) mymalloc(dustAGB.N_MASS * sizeof(double));

  my_accel_ptr = gsl_interp_accel_alloc();
  my_spline_ptr = gsl_spline_alloc(gsl_interp_linear, dustAGB.N_MASS);

  for(i = 0; i < dustAGB.N_Z; i++)
    {
      for(k = 0; k < dustAGB.N_MASS; k++)
        yield[k] = dustAGB.TotalDust[i][k] / pow(10.0, dustAGB.Mass[k]);

#ifdef TEST_YIELDS
      sprintf(filename, "AGByields.dust.table.%d", i);
      outfile = fopen(filename, "w");
      for(k = 0; k < dustAGB.N_MASS; k++)
        fprintf(outfile, "%g\t%g\n", dustAGB.Mass[k], yield[k]);
      fclose(outfile);
#endif

      gsl_spline_init(my_spline_ptr, dustAGB.Mass, yield, dustAGB.N_MASS);

      for(k = 0; k < N_MASS_BIN; k++)
	{
          if(yield_mass_bin[k] < dustAGB.Mass[0])
            result = yield[0];
          else if(yield_mass_bin[k] > dustAGB.Mass[dustAGB.N_MASS - 1])
            result = yield[dustAGB.N_MASS - 1];
          else
            result = gsl_spline_eval(my_spline_ptr, yield_mass_bin[k], my_accel_ptr);

          dustAGB.TotalDust_SPH[i][k] = pow(10.0, yield_mass_bin[k]) * result;
        }

#ifdef TEST_YIELDS
      sprintf(filename, "AGByields.dust.spline.%d", i);
      outfile = fopen(filename, "w");
      for(k = 0; k < N_MASS_BIN; k++)
        fprintf(outfile, "%g\t%g\n", yield_mass_bin[k], dustAGB.TotalDust_SPH[i][k]);
      fclose(outfile);
#endif
    }

  gsl_spline_free(my_spline_ptr);
  gsl_interp_accel_free(my_accel_ptr);
  myfree(yield);
#endif

/* POPIII yields */
#ifdef BG_POPIII
  yield = (double *) mymalloc(yieldsPOPIII.N_MASS * sizeof(double));

  my_accel_ptr = gsl_interp_accel_alloc();
  my_spline_ptr = gsl_spline_alloc(gsl_interp_linear, yieldsPOPIII.N_MASS);

  for(k = 0; k < yieldsPOPIII.N_MASS; k++)
    yield[k] = yieldsPOPIII.Ejecta[k] / pow(10.0, yieldsPOPIII.Mass[k]);

  gsl_spline_init(my_spline_ptr, yieldsPOPIII.Mass, yield, yieldsPOPIII.N_MASS);

  for(k = 0; k < POPIII_N_MASS_BIN; k++)
    {
      if(popiii_yield_mass_bin[k] < yieldsPOPIII.Mass[0])
	result = yield[0];
      else if(popiii_yield_mass_bin[k] > yieldsPOPIII.Mass[yieldsPOPIII.N_MASS - 1])
	result = yield[yieldsPOPIII.N_MASS - 1];
      else
	result = gsl_spline_eval(my_spline_ptr, popiii_yield_mass_bin[k], my_accel_ptr);

	yieldsPOPIII.Ejecta_SPH[k] = pow(10.0, popiii_yield_mass_bin[k]) * result;

	if(popiii_yield_mass_bin[k] > yieldsPOPIII.Mass[yieldsPOPIII.N_MASS - 1])
	  yieldsPOPIII.Ejecta_SPH[k] = 0;
    }

  for(k = 0; k < yieldsPOPIII.N_MASS; k++)
    yield[k] = yieldsPOPIII.TotalMetals[k] / pow(10.0, yieldsPOPIII.Mass[k]);

  gsl_spline_init(my_spline_ptr, yieldsPOPIII.Mass, yield, yieldsPOPIII.N_MASS);

  for(k = 0; k < POPIII_N_MASS_BIN; k++)
    {
      if(popiii_yield_mass_bin[k] < yieldsPOPIII.Mass[0])
	result = yield[0];
      else if(popiii_yield_mass_bin[k] > yieldsPOPIII.Mass[yieldsPOPIII.N_MASS - 1])
	result = yield[yieldsPOPIII.N_MASS - 1];
      else
	result = gsl_spline_eval(my_spline_ptr, popiii_yield_mass_bin[k], my_accel_ptr);

	yieldsPOPIII.TotalMetals_SPH[k] = pow(10.0, popiii_yield_mass_bin[k]) * result;

        if(popiii_yield_mass_bin[k] > yieldsPOPIII.Mass[yieldsPOPIII.N_MASS - 1])
	  yieldsPOPIII.TotalMetals_SPH[k] = 0;
    }

  myfree(yield);
  gsl_spline_free(my_spline_ptr);
  gsl_interp_accel_free(my_accel_ptr);

#ifdef TEST_YIELDS
      sprintf(filename, "POPIIIejecta.spline");
      outfile = fopen(filename, "w");
      for(k = 0; k < POPIII_N_MASS_BIN; k++)
	fprintf(outfile, "%g\t%g\n", popiii_yield_mass_bin[k], yieldsPOPIII.Ejecta_SPH[k]);
      fclose(outfile);

      sprintf(filename, "POPIIItotalmetals.spline");
      outfile = fopen(filename, "w");
      for(k = 0; k < POPIII_N_MASS_BIN; k++)
	fprintf(outfile, "%g\t%g\n", popiii_yield_mass_bin[k], yieldsPOPIII.TotalMetals_SPH[k]);
      fclose(outfile);
#endif

#ifdef BG_DUST
      yield = (double *) mymalloc(dustPOPIII.N_MASS * sizeof(double));

      my_accel_ptr = gsl_interp_accel_alloc();
      my_spline_ptr = gsl_spline_alloc(gsl_interp_linear, dustPOPIII.N_MASS);

      for(k = 0; k < dustPOPIII.N_MASS; k++)
	yield[k] = dustPOPIII.TotalDust[k] / pow(10.0, dustPOPIII.Mass[k]);

#ifdef TEST_YIELDS
      sprintf(filename, "POPIIIyields.dust.table");
      outfile = fopen(filename, "w");
      for(k = 0; k < dustPOPIII.N_MASS; k++)
        fprintf(outfile, "%g\t%g\n", dustPOPIII.Mass[k], yield[k]);
      fclose(outfile);
#endif

      gsl_spline_init(my_spline_ptr, dustPOPIII.Mass, yield, dustPOPIII.N_MASS);

      for(k = 0; k < POPIII_N_MASS_BIN; k++)
	{
	  if(popiii_yield_mass_bin[k] < dustPOPIII.Mass[0])
	    result = yield[0];
	  else if(popiii_yield_mass_bin[k] > dustPOPIII.Mass[dustPOPIII.N_MASS - 1])
	    result = yield[dustPOPIII.N_MASS - 1];
	  else
	    result = gsl_spline_eval(my_spline_ptr, popiii_yield_mass_bin[k], my_accel_ptr);

	  dustPOPIII.TotalDust_SPH[k] = pow(10.0, popiii_yield_mass_bin[k]) * result;

	  if(popiii_yield_mass_bin[k] > dustPOPIII.Mass[dustPOPIII.N_MASS - 1])
	    dustPOPIII.TotalDust_SPH[k] = 0;
	}

#ifdef TEST_YIELDS
      sprintf(filename, "POPIIIyields.dust.spline");
      outfile = fopen(filename, "w");
      for(k = 0; k < POPIII_N_MASS_BIN; k++)
        fprintf(outfile, "%g\t%g\n", popiii_yield_mass_bin[k], dustPOPIII.TotalDust_SPH[k]);
      fclose(outfile);
#endif

      gsl_spline_free(my_spline_ptr);
      gsl_interp_accel_free(my_accel_ptr);
      myfree(yield);
#endif

#endif /* BG_POPIII */

}


void TestYields(void)
{
  int iel;

  double star_Z;		/* metallicty Z = 1-X-Y */

  int nsteps = 1000;

  double time, tmin = 0, tmax = 0.1, dt;	/* time in Gyears */

  double C_Number_of_SNII = 0, C_Number_of_SNIa = 0;

  double Number_of_SNII = 0, Number_of_SNIa = 0;

  MyFloat metals_released[BG_NELEMENTS];

  double initial_metals[BG_NELEMENTS];

  FILE *outfile = 0;

  char filename[100];

  for(iel = 0; iel < BG_NELEMENTS; iel++)
    {
      metals_released[iel] = 0;
      initial_metals[iel] = 0;
    }

  initial_metals[0] = 0.75;
  initial_metals[1] = 0.25;

  star_Z = yieldsSNII.Metallicity[0];
  dt = (tmax - tmin) / ((float) nsteps);
  sprintf(filename, "stellar_yields.test");

  if(ThisTask == 1)
    {
      if(!(outfile = fopen(filename, "w")))
	{
	  printf("can't open file %s\n", filename);
	  endrun(8889);
	}

      fprintf(outfile, "# imf exponent by number = %g\n", All.IMF_Exponent);
      fprintf(outfile, "# time [Gyr] ");

      for(iel = 0; iel < BG_NELEMENTS; iel++)
	fprintf(outfile, " %s\t ", ElementNames[iel]);
    }

  for(time = tmin; time < tmax; time += dt)
    {
      /*      stellar_evolution(time, dt, star_Z, initial_metals, metals_released, &mass_released, &Number_of_SNIa, &Number_of_SNII); */

      C_Number_of_SNII += Number_of_SNII;
      C_Number_of_SNIa += Number_of_SNIa;

      if(ThisTask == 1)
	{
	  printf(" time = %g Number_of_SNII = %g\n", time, C_Number_of_SNII);

	  fprintf(outfile, " %g ", time);

	  for(iel = 0; iel < BG_NELEMENTS; iel++)
	    fprintf(outfile, " %g\t ", metals_released[iel]);

	  fprintf(outfile, "\n");
	}
    }

  if(ThisTask == 1)
    fclose(outfile);
}


#endif
