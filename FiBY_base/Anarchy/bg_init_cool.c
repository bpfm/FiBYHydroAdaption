#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>
#include <mpi.h>

#if defined(BG_COOLING) && !defined(BG_COOLING_OLD_TABLES)

#include "allvars.h"
#include "proto.h"
#include "bg_cooling.h"
#include "bg_proto.h"
#include "bg_vars.h"


/*
 * ----------------------------------------------------------------------
 * Prototype for built-in testing routines
 * ----------------------------------------------------------------------
 */

void DriverTestCool();

void TestCool(float z);

void DriverTestCoolThermalEvolution();


/*
 * ----------------------------------------------------------------------
 * This routine gets the redshifts from a text file
 * ----------------------------------------------------------------------
 */

void GetCoolingRedshifts()
{
  FILE *infile;

  int i = 0;

  char buffer[500], redfilename[500];


  sprintf(redfilename, "%s/redshifts.dat", All.CoolTablePath);
  infile = fopen(redfilename, "r");
  if(infile == NULL)
    puts("GetCoolingRedshifts can't open a file");

  if(fscanf(infile, "%s", buffer) != EOF)
    {
      cooling_N_Redshifts = atoi(buffer);
      cooling_Redshifts = (float *) mymalloc(cooling_N_Redshifts * sizeof(float));

      while(fscanf(infile, "%s", buffer) != EOF)
	{
	  cooling_Redshifts[i] = atof(buffer);
	  i += 1;
	}
    }
  fclose(infile);

#ifdef BG_VERBOSE
  for(i = 0; i < cooling_N_Redshifts; i++)
    printf("Cooling table %2d at redshift %1.3f\n", i + 1, cooling_Redshifts[i]);
  printf("\n");
#endif
}


/*
 * ----------------------------------------------------------------------
 * This routine reads in the header of the cooling table files
 * ----------------------------------------------------------------------
 */

void ReadCoolingHeader(char *fname)
{
  int i;

  hid_t tempfile_id, dataset, datatype;

  herr_t status;


  /* fill the constants */
  tempfile_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(tempfile_id < 0)
    {
      printf("[ReadCoolingHeader()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_density_bins");
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cooling_N_nH);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_temperature_bins");
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cooling_N_Temp);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_helium_fractions");
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cooling_N_He);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_metals");
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cooling_N_Elements);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Number_of_abundances");
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cooling_N_SolarAbundances);
  status = H5Dclose(dataset);

  /* allocate arrays for cooling table header */
  allocate_header_arrays();

  /* fill the arrays */
  dataset = H5Dopen(tempfile_id, "/Solar/Temperature_bins");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling_Temp);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Solar/Hydrogen_density_bins");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling_nH);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Metal_free/Temperature/Energy_density_bins");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling_Therm);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Metal_free/Helium_mass_fraction_bins");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling_HeFrac);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Solar_mass_fractions");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling_SolarAbundances);
  status = H5Dclose(dataset);

  char element_names[cooling_N_Elements][EL_NAME_LENGTH];
  hsize_t string_length = EL_NAME_LENGTH;

  /* names of chemical elements stored in table */
  datatype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(datatype, string_length);
  dataset = H5Dopen(tempfile_id, "/Header/Metal_names");
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, element_names);
  status = H5Dclose(dataset);

  for(i = 0; i < cooling_N_Elements; i++)
    cooling_ElementNames[i] = mystrdup(element_names[i]);

  char solar_abund_names[cooling_N_SolarAbundances][EL_NAME_LENGTH];

  /* assumed solar abundances used in constructing the tables, and corresponding names */
  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Abund_names");
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, solar_abund_names);
  status = H5Dclose(dataset);
  H5Tclose(datatype);

  for(i = 0; i < cooling_N_SolarAbundances; i++)
    cooling_SolarAbundanceNames[i] = mystrdup(solar_abund_names[i]);

  status = H5Fclose(tempfile_id);

  /* Convert to temperature, density and internal energy arrays to log10 */
  for(i = 0; i < cooling_N_Temp; i++)
    {
      cooling_Temp[i] = log10(cooling_Temp[i]);
      cooling_Therm[i] = log10(cooling_Therm[i]);
    }

  for(i = 0; i < cooling_N_nH; i++)
    cooling_nH[i] = log10(cooling_nH[i]);

  printf("Done with cooling table header.\n");
  fflush(stdout);

#ifdef BG_VERBOSE
  printf("\n\nTemperature grid (LOG10) for energy to temperature conversion\n\n");
  for(i = 0; i < cooling_N_Temp; i++)
    printf("Temperature[%3d] = %f\n", i, cooling_Temp[i]);

  printf("\n\nEnergy grid (LOG10) for energy to temperature conversion\n\n");
  for(i = 0; i < cooling_N_Temp; i++)
    printf("Energy[%3d] = %f\n", i, cooling_Therm[i]);

  printf("\n\nDensity grid (LOG10) for energy to temperature conversion\n\n");
  for(i = 0; i < cooling_N_nH; i++)
    printf("Density[%2d] = %f\n", i, cooling_nH[i]);
#endif
}


/*
 * ----------------------------------------------------------------------
 * This routine broadcasts the header infos
 * ----------------------------------------------------------------------
 */

void BroadCastHeader()
{
  int position, bufsize, size, i;

  char *buffer;


  /* get size of the buffer */
  MPI_Pack_size(6, MPI_INT, MPI_COMM_WORLD, &size);
  bufsize = size;

  /* allocate memory for the buffer */
  buffer = (char *) mymalloc(bufsize);

  if(ThisTask == 0)
    {
      position = 0;

      /* pack the array dimensions */
      MPI_Pack(&cooling_N_nH, 1,
	       MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(&cooling_N_Temp, 1,
	       MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(&cooling_N_He, 1,
	       MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(&cooling_N_Elements, 1,
	       MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(&cooling_N_SolarAbundances, 1,
	       MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(&cooling_N_Redshifts, 1,
	       MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
    }

  MPI_Bcast(buffer, bufsize, MPI_PACKED, 0, MPI_COMM_WORLD);

  if(ThisTask != 0)
    {
      position = 0;

      /* unpack the array dimensions */
      MPI_Unpack(buffer, bufsize, &position, &cooling_N_nH, 1,
		 MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, &cooling_N_Temp, 1,
		 MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, &cooling_N_He, 1,
		 MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, &cooling_N_Elements, 1,
		 MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, &cooling_N_SolarAbundances, 1,
		 MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, &cooling_N_Redshifts, 1,
		 MPI_INT, MPI_COMM_WORLD);
    }

  /* free allocated memory */
  myfree(buffer);

  if(ThisTask != 0)
    {
      /* allocate header arrays */
      allocate_header_arrays();

      /* allocate redshift array */
      cooling_Redshifts = (float *) mymalloc(cooling_N_Redshifts * sizeof(float));
    }

  bufsize = 0;
  /* get size of the buffer */
  MPI_Pack_size(cooling_N_Temp,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(cooling_N_nH,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(cooling_N_Temp,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(cooling_N_He,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(cooling_N_SolarAbundances,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(cooling_N_Redshifts, MPI_FLOAT,
		MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(cooling_N_Elements * EL_NAME_LENGTH,
		MPI_CHAR, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(cooling_N_SolarAbundances * EL_NAME_LENGTH,
		MPI_CHAR, MPI_COMM_WORLD, &size);
  bufsize += size;

  /* allocate memory for the buffer */
  buffer = (char *) mymalloc(bufsize);

  if(ThisTask == 0)
    {
      position = 0;

      MPI_Pack(cooling_Temp, cooling_N_Temp,
	       MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(cooling_nH, cooling_N_nH,
	       MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(cooling_Therm, cooling_N_Temp,
	       MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(cooling_HeFrac, cooling_N_He,
	       MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(cooling_SolarAbundances, cooling_N_SolarAbundances,
	       MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(cooling_Redshifts, cooling_N_Redshifts,
	       MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);

      for(i = 0; i < cooling_N_Elements; i++)
	MPI_Pack(cooling_ElementNames[i], EL_NAME_LENGTH,
		 MPI_CHAR, buffer, bufsize, &position, MPI_COMM_WORLD);
      for(i = 0; i < cooling_N_SolarAbundances; i++)
	MPI_Pack(cooling_SolarAbundanceNames[i], EL_NAME_LENGTH,
		 MPI_CHAR, buffer, bufsize, &position, MPI_COMM_WORLD);
    }

  MPI_Bcast(buffer, bufsize, MPI_PACKED, 0, MPI_COMM_WORLD);

  if(ThisTask != 0)
    {
      position = 0;

      MPI_Unpack(buffer, bufsize, &position, cooling_Temp, cooling_N_Temp,
		 MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, cooling_nH, cooling_N_nH,
		 MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, cooling_Therm, cooling_N_Temp,
		 MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, cooling_HeFrac, cooling_N_He,
		 MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, cooling_SolarAbundances, cooling_N_SolarAbundances,
		 MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, cooling_Redshifts, cooling_N_Redshifts,
		 MPI_FLOAT, MPI_COMM_WORLD);

      for(i = 0; i < cooling_N_Elements; i++)
	MPI_Unpack(buffer, bufsize, &position, cooling_ElementNames[i], EL_NAME_LENGTH,
		   MPI_CHAR, MPI_COMM_WORLD);
      for(i = 0; i < cooling_N_SolarAbundances; i++)
	MPI_Unpack(buffer, bufsize, &position, cooling_SolarAbundanceNames[i], EL_NAME_LENGTH,
		   MPI_CHAR, MPI_COMM_WORLD);
    }

  /* free allocated memory */
  myfree(buffer);
}


/*
 * ----------------------------------------------------------------------
 * This routine sets the value of solar metallicity using solar element
 * abundances stored in the cooling tables.
 * ----------------------------------------------------------------------
 */

void set_solar_metallicity(void)
{
  int i;

#ifdef BG_VERBOSE
  if(ThisTask == 0)
    printf("\n\nSolar abundances:\n\n");
#endif

  /* set solar metallicity value */
  for(i = 0, SolarMetallicity = 1; i < cooling_N_SolarAbundances; i++)
    {
#ifdef BG_VERBOSE
      if(ThisTask == 0)
	printf("Solar abundance of %s = %g\n",
	       cooling_SolarAbundanceNames[i],
	       cooling_SolarAbundances[i]);
#endif
      if(strcmp(cooling_SolarAbundanceNames[i], "Hydrogen") == 0)
	SolarMetallicity -= cooling_SolarAbundances[i];

      if(strcmp(cooling_SolarAbundanceNames[i], "Helium") == 0)
	SolarMetallicity -= cooling_SolarAbundances[i];
    }

#ifdef BG_VERBOSE
  if(ThisTask == 0)
    printf("\nSolarMetallicity = %g\n", SolarMetallicity);
#endif
}


/*
 * ----------------------------------------------------------------------
 * Make pointers to element names
 * ----------------------------------------------------------------------
 */

void MakeNamePointers()
{
  int i, j, sili_index = 0;

  ElementNamePointers = (int *) mymalloc(cooling_N_Elements * sizeof(int));
  SolarAbundanceNamePointers = (int *) mymalloc(cooling_N_Elements * sizeof(int));

  for(i = 0; i < BG_NELEMENTS; i++)
    {
      if(strcmp(ElementNames[i], "Silicon") == 0)
	sili_index = i;
    }

  for(i = 0; i < cooling_N_Elements; i++)
    {
      SolarAbundanceNamePointers[i] = -999;
      ElementNamePointers[i] = -999;

      for(j = 0; j < cooling_N_SolarAbundances; j++)
	{
	  if(strcmp(cooling_ElementNames[i], cooling_SolarAbundanceNames[j]) == 0)
	    SolarAbundanceNamePointers[i] = j;
	}

      if(strcmp(cooling_ElementNames[i], "Sulphur") == 0 ||
	 strcmp(cooling_ElementNames[i], "Calcium") == 0)	/* These elements are tracked! */
	ElementNamePointers[i] = -1 * sili_index;
      else
	{
	  for(j = 0; j < BG_NELEMENTS; j++)
	    {
	      if(strcmp(cooling_ElementNames[i], ElementNames[j]) == 0)
		ElementNamePointers[i] = j;
	    }
	}
    }

#ifdef BG_VERBOSE
  if(ThisTask == 0)
    {
      for(i = 0; i < cooling_N_Elements; i++)
	printf("cooling_ElementNames[%d] is cooling_SolarAbundanceNames %s\n", i,
	       cooling_SolarAbundanceNames[SolarAbundanceNamePointers[i]]);
      printf("\n");

      for(i = 0; i < cooling_N_Elements; i++)
	printf("cooling_ElementNames[%d] is ElementNames %s\n", i,
	       ElementNames[ElementNamePointers[i]]);
      printf("\n");
    }
#endif
}


/*
 * ----------------------------------------------------------------------
 * This routine initialize cooling at the beginning of the simulation
 * ----------------------------------------------------------------------
 */

void InitCool()
{
  int i, j, k;

  double Z_Solar_Calcium = 0, Z_Solar_Sulphur = 0, Z_Solar_Silicon = 0;

  char fname[500];


  /* use a general file to get the header */
  if(ThisTask == 0)
    {
      GetCoolingRedshifts();
      sprintf(fname, "%sz_0.000.hdf5", All.CoolTablePath);
      ReadCoolingHeader(fname);
    }

  BroadCastHeader();

  set_solar_metallicity();

  if(ThisTask == 0 && element_present("Silicon") == 0)
    {
      for(i = 0; i < cooling_N_SolarAbundances; i++)
	{
	  if(strcmp(cooling_SolarAbundanceNames[i], "Calcium") == 0)
	    Z_Solar_Calcium = cooling_SolarAbundances[i];

	  if(strcmp(cooling_SolarAbundanceNames[i], "Sulphur") == 0)
	    Z_Solar_Sulphur = cooling_SolarAbundances[i];

	  if(strcmp(cooling_SolarAbundanceNames[i], "Silicon") == 0)
	    Z_Solar_Silicon = cooling_SolarAbundances[i];
	}
    }
  else
    {
      if(ThisTask == 0)
	printf("Silicon is not present: Calcium and Sulphur will not be tracked!");
    }

  cooling_ThermalToTemp = (float ****) mymalloc(2 * sizeof(float ***));

  for(i = 0; i < 2; i++)
    {
      cooling_ThermalToTemp[i] =
	(float ***) mymalloc(cooling_N_He * sizeof(float **));

      for(j = 0; j < cooling_N_He; j++)
	{
	  cooling_ThermalToTemp[i][j] =
	    (float **) mymalloc(cooling_N_nH * sizeof(float *));

	  for(k = 0; k < cooling_N_nH; k++)
	    cooling_ThermalToTemp[i][j][k] =
	      (float *) mymalloc(cooling_N_Temp * sizeof(float));
	}
    }

  cooling_SolarElectronAbundance =
    (float ***) mymalloc(2 * sizeof(float **));

  for(i = 0; i < 2; i++)
    {
      cooling_SolarElectronAbundance[i] = 
	(float **) mymalloc(cooling_N_nH * sizeof(float *));

      for(j = 0; j < cooling_N_nH; j++)
	cooling_SolarElectronAbundance[i][j] =
	  (float *) mymalloc(cooling_N_Temp * sizeof(float));
    }

  /* Allocate for metal cooling */
  cooling_MetalsNetHeating = (float ****) mymalloc(2 * sizeof(float ***));

  for(i = 0; i < 2; i++)
    {
      cooling_MetalsNetHeating[i] =
	(float ***) mymalloc(cooling_N_Elements * sizeof(float **));

      for(j = 0; j < cooling_N_Elements; j++)
	{
	  cooling_MetalsNetHeating[i][j] =
	    (float **) mymalloc(cooling_N_nH * sizeof(float *));

	  for(k = 0; k < cooling_N_nH; k++)
	    cooling_MetalsNetHeating[i][j][k] =
	      (float *) mymalloc(cooling_N_Temp * sizeof(float));
	}
    }

  /* Allocate for H and He cooling */
  cooling_HplusHeNetHeating = (float ****) mymalloc(2 * sizeof(float ***));

  for(i = 0; i < 2; i++)
    {
      cooling_HplusHeNetHeating[i] =
	(float ***) mymalloc(cooling_N_He * sizeof(float **));

      for(j = 0; j < cooling_N_He; j++)
	{
	  cooling_HplusHeNetHeating[i][j] =
	    (float **) mymalloc(cooling_N_nH * sizeof(float *));

	  for(k = 0; k < cooling_N_nH; k++)
	    cooling_HplusHeNetHeating[i][j][k] =
	      (float *) mymalloc(cooling_N_Temp * sizeof(float));
	}
    }

  cooling_HplusHeElectronAbundance = (float ****) mymalloc(2 * sizeof(float ***));

  for(i = 0; i < 2; i++)
    {
      cooling_HplusHeElectronAbundance[i] =
	(float ***) mymalloc(cooling_N_He * sizeof(float **));

      for(j = 0; j < cooling_N_He; j++)
	{
	  cooling_HplusHeElectronAbundance[i][j] =
	    (float **) mymalloc(cooling_N_nH * sizeof(float *));

	  for(k = 0; k < cooling_N_nH; k++)
	    cooling_HplusHeElectronAbundance[i][j][k] =
	      (float *) mymalloc(cooling_N_Temp * sizeof(float));
	}
    }

  if(ThisTask == 0)
    printf("Allocating memory for extra collisional tables\n");

  if(ThisTask == 0)
    printf("%d %d %d %d\n", cooling_N_Elements, cooling_N_He, cooling_N_nH, cooling_N_Temp);

#ifdef BG_COOLING_SHIELDING

  cooling_ThermalToTempCollisional =
    (float ***) mymalloc(cooling_N_He * sizeof(float **));

  for(j = 0; j < cooling_N_He; j++)
    {
      cooling_ThermalToTempCollisional[j] =
	(float **) mymalloc(cooling_N_nH * sizeof(float *));

      for(k = 0; k < cooling_N_nH; k++)
	cooling_ThermalToTempCollisional[j][k] =
	  (float *) mymalloc(cooling_N_Temp * sizeof(float));
    }

  cooling_SolarElectronAbundanceCollisional = 
    (float **) mymalloc(cooling_N_nH * sizeof(float *));

  for(j = 0; j < cooling_N_nH; j++)
    cooling_SolarElectronAbundanceCollisional[j] =
      (float *) mymalloc(cooling_N_Temp * sizeof(float));

  /* Allocate for metal cooling */

  cooling_MetalsNetHeatingCollisional =
    (float ***) mymalloc(cooling_N_Elements * sizeof(float **));

  for(j = 0; j < cooling_N_Elements; j++)
    {
      cooling_MetalsNetHeatingCollisional[j] =
	(float **) mymalloc(cooling_N_nH * sizeof(float *));

      for(k = 0; k < cooling_N_nH; k++)
	cooling_MetalsNetHeatingCollisional[j][k] =
	  (float *) mymalloc(cooling_N_Temp * sizeof(float));
    }

  /* Allocate for H and He cooling */

  cooling_HplusHeNetHeatingCollisional =
    (float ***) mymalloc(cooling_N_He * sizeof(float **));

  for(j = 0; j < cooling_N_He; j++)
    {
      cooling_HplusHeNetHeatingCollisional[j] =
	(float **) mymalloc(cooling_N_nH * sizeof(float *));

      for(k = 0; k < cooling_N_nH; k++)
	cooling_HplusHeNetHeatingCollisional[j][k] =
	  (float *) mymalloc(cooling_N_Temp * sizeof(float));
    }

  cooling_HplusHeElectronAbundanceCollisional =
    (float ***) mymalloc(cooling_N_He * sizeof(float **));
  
  for(j = 0; j < cooling_N_He; j++)
    {
      cooling_HplusHeElectronAbundanceCollisional[j] =
	(float **) mymalloc(cooling_N_nH * sizeof(float *));

      for(k = 0; k < cooling_N_nH; k++)
	cooling_HplusHeElectronAbundanceCollisional[j][k] =
	  (float *) mymalloc(cooling_N_Temp * sizeof(float));
    }
#endif

  MakeNamePointers();

  /* Set H and He reionisation heating to erg/g */
  All.REION_H_Heating_ERGpG = All.REION_H_Heating_EVpH * EV_TO_ERG / PROTONMASS;
  All.REION_He_Heating_ERGpG = All.REION_He_Heating_EVpH * EV_TO_ERG / PROTONMASS;

  if(ThisTask == 0)
    {
      printf("Hydrogen reionisation heating = %e erg/g\n", All.REION_H_Heating_ERGpG);
      printf("  Helium reionisation heating = %e erg/g\n", All.REION_He_Heating_ERGpG);
    }

  cooling_ElementAbundance_SOLAR =
    (double *) mymalloc(cooling_N_Elements * sizeof(double));

  cooling_ElementAbundance_SOLARM1 = 
    (double *) mymalloc(cooling_N_Elements * sizeof(double));

#ifdef BG_MOL_COOLING
  hid_t file_id, dataset;
  herr_t status;

  /* ------------------------------------------------------------------------- */
  /*   read H2 cooling table                                                   */
  /* ------------------------------------------------------------------------- */

  sprintf(fname, "%sH2.hdf5", All.CoolTablePath);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset = H5Dopen(file_id, "Number_of_temperature_bins");
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mol_cooling_H2_N_Temp);
  status = H5Dclose(dataset);

  mol_cooling_H2_Temp = (float *) mymalloc(mol_cooling_H2_N_Temp * sizeof(float));
  mol_cooling_H2_Rate = (float *) mymalloc(mol_cooling_H2_N_Temp * sizeof(float));

  dataset = H5Dopen(file_id, "Temperature");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mol_cooling_H2_Temp);
  status = H5Dclose(dataset);

  dataset = H5Dopen(file_id, "CoolingRate");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mol_cooling_H2_Rate);
  status = H5Dclose(dataset);

  status = H5Fclose(file_id);

  for(i = 0; i < mol_cooling_H2_N_Temp; i++)
    mol_cooling_H2_Temp[i] = log10(mol_cooling_H2_Temp[i]);

  /* ------------------------------------------------------------------------- */
  /*   read HD cooling table                                                   */
  /* ------------------------------------------------------------------------- */

  sprintf(fname, "%sHD.hdf5", All.CoolTablePath);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset = H5Dopen(file_id, "Number_of_temperature_bins");
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mol_cooling_HD_N_Temp);
  status = H5Dclose(dataset);

  mol_cooling_HD_Temp = (float *) mymalloc(mol_cooling_HD_N_Temp * sizeof(float));
  mol_cooling_HD_Rate = (float *) mymalloc(mol_cooling_HD_N_Temp * sizeof(float));

  dataset = H5Dopen(file_id, "Temperature");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mol_cooling_HD_Temp);
  status = H5Dclose(dataset);

  dataset = H5Dopen(file_id, "CoolingRate");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mol_cooling_HD_Rate);
  status = H5Dclose(dataset);

  status = H5Fclose(file_id);

  for(i = 0; i < mol_cooling_HD_N_Temp; i++)
    mol_cooling_HD_Temp[i] = log10(mol_cooling_HD_Temp[i]);
#endif

#if BG_TEST_COOL == 1
  if(ThisTask == 0)
    printf("BG_TEST_COOL option 1: ThisTask = 0 will evaluate cooling rates, then we will stop\n");

  DriverTestCool();

  MPI_Barrier(MPI_COMM_WORLD);
  printf("Test 1 of cooling finished\n");
  endrun(-2);
#endif

#if BG_TEST_COOL == 2
  if(ThisTask == 0)
    printf("BG_TEST_COOL option 2: ThisTask = 0 will evaluate thermal evolution at the mean density\n");

  DriverTestCoolThermalEvolution();

  MPI_Barrier(MPI_COMM_WORLD);
  printf("Test 2 of cooling finished\n");
  endrun(-2);
#endif
}


/*
 * ----------------------------------------------------------------------
 * This routine allocates header arrays for the cooling tables
 * ----------------------------------------------------------------------
 */

void allocate_header_arrays(void)
{
  int i;


  /* bins for interpolation in log10(T [K]) */
  cooling_Temp = (float *) mymalloc(cooling_N_Temp * sizeof(float));

  /* bins for interpolation in hydrogen number density, log10(nh  [/cm^3]) */
  cooling_nH = (float *) mymalloc(cooling_N_nH * sizeof(float));

  /* bins for interpolation in thermal energy per unit mass, log10(u [erg/g]) */
  cooling_Therm = (float *) mymalloc(cooling_N_Temp * sizeof(float));

  /* bins for interpolation in Helium abundance */
  cooling_HeFrac = (float *) mymalloc(cooling_N_He * sizeof(float));

  /* table of solar abundances */
  cooling_SolarAbundances = (float *) mymalloc(cooling_N_SolarAbundances * sizeof(float));

  /* assumed solar abundances used in constructing the tables, and corresponding names */
  cooling_SolarAbundanceNames = (char **) mymalloc(cooling_N_SolarAbundances * sizeof(char *));
  if(ThisTask != 0)
    for(i = 0; i < cooling_N_SolarAbundances; i++)
      cooling_SolarAbundanceNames[i] = (char *) mymalloc(EL_NAME_LENGTH * sizeof(char));

  /* names of chemical elements stored in table */
  cooling_ElementNames = (char **) mymalloc(cooling_N_Elements * sizeof(char *));
  if(ThisTask != 0)
    for(i = 0; i < cooling_N_Elements; i++)
      cooling_ElementNames[i] = (char *) mymalloc(EL_NAME_LENGTH * sizeof(char));
}


/*
 * ======================================================================
 * FUNCTIONS FOR BASIC TESTING OF COOLING ROUTINES - TEST_COOL flag.
 * ======================================================================
 */

/*
 * ----------------------------------------------------------------------
 * This routine tests the interpolation of the cooling tables
 * ----------------------------------------------------------------------
 */

void DriverTestCool()
{
  int iz;

  int z_index;

  float z1, z2, dz;


  /* we use redshifts at exact values from the tables, */
  /* as well as values in between */
  for(iz = cooling_N_Redshifts - 2; iz > 0; iz--)
    {
      z1 = cooling_Redshifts[iz + 1];	/* collisional table */
      z2 = cooling_Redshifts[iz];	/* collisional table */

      get_redshift_index(z1, &z_index, &dz);
      LoadCoolingTables(z_index);

      /*
	TestCool(z1 - 0.1*(z1-z2));
      */

      if(ThisTask == 0)
	TestCool(z1);
    }
}


void TestCool(float z)
{
  FILE *outfile;

  int i, j, k, index1, index2, H_index, He_index, z_index;

  char filename[500];

  double Hefrac, ener, dens, temp;

  static double ParticleMetalFraction[BG_NELEMENTS];

  float dz;


  /* given Helium abundance */
  Hefrac = cooling_HeFrac[0];

  get_redshift_index(z, &z_index, &dz);

  H_index = element_index("Hydrogen");
  He_index = element_index("Helium");

  /* loop over elements */
  for(i = 0; i < BG_NELEMENTS; i++)
    {
      /* put abundance of element i to solar value, zero all others  */
      for(j = 0; j < BG_NELEMENTS; j++)
	ParticleMetalFraction[j] = 0;

      index1 = get_element_index(cooling_ElementNames,
				 cooling_N_Elements,
				 ElementNames[i]);

      if(index1 >= 0)  /* index < 0 for Helium which is not a metal */
	{
	  index2 = get_element_index(cooling_SolarAbundanceNames,
				     cooling_N_SolarAbundances,
				     cooling_ElementNames[index1]);

	  if(index2 < 0 || index2 >= cooling_N_SolarAbundances)
	    {
	      printf(" sorry: element %s not found \n", cooling_ElementNames[index1]);
	      endrun(-12344);
	    }

	  sprintf(filename, "testcool.z=%1.3f.Helium=%1.3f.Element=%s",
		  z, Hefrac, ElementNames[i]);

	  outfile = fopen(filename, "w");

	  ParticleMetalFraction[H_index] = 1.0 - Hefrac - cooling_SolarAbundances[index2];
	  ParticleMetalFraction[He_index] = Hefrac;

	  ParticleMetalFraction[i] = cooling_SolarAbundances[index2];

	  /* write temperature density grid */
	  fprintf(outfile, "# element = %s\n", ElementNames[i]);
	  fprintf(outfile, "# redshift = %e\n", z);
	  fprintf(outfile, "# Helium abundance = %e\n", Hefrac);
	  fprintf(outfile, "# cooling_N_Temp = %d\n", cooling_N_Temp);
	  fprintf(outfile, "# cooling_N_nH = %d\n", cooling_N_nH);

	  for(j = 0; j < cooling_N_Temp - 1; j++)
	    {
	      ener = cooling_Therm[j];

	      for(k = 0; k < cooling_N_nH - 1; k++)
		{
		  dens = cooling_nH[k];

		  temp = convert_u_to_temp(dz, pow(10.0, ener), pow(10.0, dens), Hefrac);

		  fprintf(outfile, " %e  %e  %e\n", temp, pow(10.0, dens),
			  CoolingRate(pow(10.0, ener), pow(10.0, dens), z, dz, ParticleMetalFraction, 0, 0));
		}
	    }

	  fclose(outfile);
	}
    }
}


/*
 * ----------------------------------------------------------------------
 * This routine tests the thermal evolution of gas at the critical density
 * ----------------------------------------------------------------------
 */

void DriverTestCoolThermalEvolution()
{
  FILE *outfile = 0;

  char filename[500];

  int i, HeIndex, HIndex;

  int z_index;

  float z, dz;

  double aini = 1.0 / 50, a, dloga = 0.01, uinit = 5.0e10, delta = 0.0, rho, u;

  double RhoCrit, hubble_z, dtime, deltaz;

  double element_metallicity[BG_NELEMENTS];

  double Hefrac = 0.248, Hfrac = 0.752, temp;

  double uc = 0;


  for(i = 0; i < BG_NELEMENTS; i++)
    element_metallicity[i] = 0;

  HeIndex = element_index("Helium");
  HIndex = element_index("Hydrogen");

  element_metallicity[HeIndex] = Hefrac;
  element_metallicity[HIndex] = Hfrac;

  /* crit. density in code units */
  RhoCrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  /* crit. density in physical g/cm^3 */
  RhoCrit *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

  u = uinit;

  if(ThisTask == 0)
    {
      sprintf(filename, "thermevol.data");

      if(!(outfile = fopen(filename, "w")))
	{
	  printf("can't open file %s\n", filename);
	  endrun(888);
	}

      fprintf(outfile, "# delta = %g\n", delta);
      fprintf(outfile, "# uinit = %g\n", uinit);

      fprintf(outfile, "# metallicity = ");

      for(i = 0; i < BG_NELEMENTS; i++)
	fprintf(outfile, " %g\t", element_metallicity[i]);

      fprintf(outfile, "\n");
      fprintf(outfile, "# z T[K] u[erg/g] rho[g/cm^3]\n");
    }

#ifdef BG_COOLING_SHIELDING
  LoadCollisionalCoolingTables();
#endif

  double time = 0;

  for(a = aini; log(a) + dloga < 0; a += a * dloga)
    {
      z = 1 / a - 1;

      //deltaz = 1 / exp(log(a) + dloga) - 1 / a;
      deltaz = 1 / a - 1 / exp(log(a) - dloga);

      if(ThisTask == 0)
	{
	  printf("  REDSHIFT = %f\n", z);
	  fflush(stdout);
	}

      u *= pow(((1 + z + deltaz) / (1 + z)), 2);

      get_redshift_index(z, &z_index, &dz);
      LoadCoolingTables(z_index);

      /* physical baryon density at this z, for this over density, in g/cm^3 */
      rho = (1 + delta) * RhoCrit * All.OmegaBaryon * pow(1 + z, 3);

      /* rho *= 1e6; */

      /* da = -deltaz / pow(1 + z, 2); */

      hubble_z = (All.HubbleParam * HUBBLE) *
	sqrt(All.Omega0 * pow(1 + z, 3) +
	     (1 - All.OmegaLambda - All.Omega0) * pow(1 + z, 2) + All.OmegaLambda);

      dtime = dloga / hubble_z;
      time += dtime;

      printf("time = %g, dtime = %g, a = %g, dz = %g\n", time / SEC_PER_YEAR, dtime / SEC_PER_YEAR, a, deltaz);

#ifdef BG_MOL_COOLING
      u = DoCooling(u, rho, dtime, deltaz, z, dz, element_metallicity, 0, 0);
#else
      u = DoCooling(u, rho, dtime, deltaz, z, dz, element_metallicity);
#endif
      uc = uc + Reionization(z, deltaz, 0);

      //temp = convert_u_to_temp(dz, u, rho, Hefrac);
      temp = convert_u_to_temp(dz, u, rho * Hfrac / PROTONMASS, Hefrac);

      if(ThisTask == 0)
	fprintf(outfile, " %g \t %g \t %g \t %g \t %g\n", z, temp, u, rho, uc);
    }

  if(ThisTask == 0)
    fclose(outfile);
}

#endif
