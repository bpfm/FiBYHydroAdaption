#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <hdf5.h>
#include <mpi.h>

#if defined(BG_COOLING) && defined(BG_COOLING_OLD_TABLES)

#include "allvars.h"
#include "proto.h"
#include "bg_cooling_old.h"
#include "bg_proto.h"


/*
 * ----------------------------------------------------------------------
 * Get cooling table redshift index
 * ----------------------------------------------------------------------
 */

void get_redshift_index(float z, int *z_index, float *dz)
{
  int i, iz;

  static int first_call = 0;

  static int previous = -1;


  if(first_call == 0)
    {
      first_call = 1;
      previous = cooling_N_Redshifts - 2;

      /* this routine assumes cooling_redshifts table is in increasing order. Test this. */
      for(i = 0; i < cooling_N_Redshifts - 2; i++)
	if(cooling_Redshifts[i + 1] < cooling_Redshifts[i])
	  {
	    printf("[get_redshift_index]: table should be in increasing order\n");
	    endrun(-4);
	  }
    }

  /* before the earliest redshift or before hydrogen reionization, flag for collisional cooling */
  if(z > All.REION_H_ZCenter)
    {
      *z_index = cooling_N_Redshifts;
      *dz = 0.0;
    }
  /* from reionization use the cooling tables */
  else if(z > cooling_Redshifts[cooling_N_Redshifts - 1] && z <= All.REION_H_ZCenter)
    {
      *z_index = cooling_N_Redshifts + 1;
      *dz = 0.0;
    }
  /* at the end, just use the last value */
  else if(z <= cooling_Redshifts[0])
    {
      *z_index = 0;
      *dz = 0.0;
    }
  else
    {
      /* start at the previous index and search */
      for(iz = previous; iz >= 0; iz--)
	{
	  if(z > cooling_Redshifts[iz])
	    {
	      *dz = (z - cooling_Redshifts[iz]) /
		(cooling_Redshifts[iz + 1] - cooling_Redshifts[iz]);

	      previous = *z_index = iz;

	      break;
	    }
	}
    }
}

 
/*
 * ----------------------------------------------------------------------
 * Get the cooling table for photoionized cooling (before redshift ~9)
 * ----------------------------------------------------------------------
 */

void GetNoComptTable()
{
  hid_t file_id, dataset;

  herr_t status;

  char fname[100], set_name[100];

  int Hes, specs, j, k;

  float heating_rate[cooling_N_nH][cooling_N_Temp];
  float cooling_rate[cooling_N_nH][cooling_N_Temp];
  float temperature[cooling_N_nH][cooling_N_Temp];
  float electron_abundance[cooling_N_nH][cooling_N_Temp];


  sprintf(fname, "%sz_8.989nocompton.hdf5", All.CoolTablePath);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  printf("Redshift 1 %ld %s\n", (long int) file_id, fname);
  fflush(stdout);

  /* For normal elements */
  for(specs = 0; specs < cooling_N_Elements; specs++)
    {
      sprintf(set_name, "/%s/Heating", cooling_ElementNames[specs]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, heating_rate);
      status = H5Dclose(dataset);

      sprintf(set_name, "/%s/Cooling", cooling_ElementNames[specs]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling_rate);
      status = H5Dclose(dataset);

      for(j = 0; j < cooling_N_Temp; j++)
	for(k = 0; k < cooling_N_nH; k++)
	  cooling_MetalsNetHeating[0][specs][k][j] = heating_rate[k][j] - cooling_rate[k][j];
    }

  /* Helium */
  for(Hes = 0; Hes < cooling_N_He; Hes++)
    {
      sprintf(set_name, "/Metal_free/%s/Heating", cooling_HeNames[Hes]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, heating_rate);
      status = H5Dclose(dataset);

      sprintf(set_name, "/Metal_free/%s/Cooling", cooling_HeNames[Hes]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling_rate);
      status = H5Dclose(dataset);

      sprintf(set_name, "/Metal_free/%s/Temperature", cooling_HeNames[Hes]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temperature);
      status = H5Dclose(dataset);

      sprintf(set_name, "/Metal_free/%s/Electron_density_over_n_h", cooling_HeNames[Hes]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, electron_abundance);
      status = H5Dclose(dataset);

      for(j = 0; j < cooling_N_Temp; j++)
	{
	  for(k = 0; k < cooling_N_nH; k++)
	    {
	      cooling_HplusHeNetHeating[0][Hes][k][j] = heating_rate[k][j] - cooling_rate[k][j];
	      cooling_ThermalToTemp[0][Hes][k][j] = log10(temperature[k][j]);
	      cooling_CollisionalElectronAbundance[Hes][k][j] = electron_abundance[k][j];
	    }
	}
    }

  status = H5Fclose(file_id);
}


/*
 * ----------------------------------------------------------------------
 * Get the cooling table for collisional cooling (before reionisation)
 * ----------------------------------------------------------------------
 */

void GetCollisTable()
{
  hid_t file_id, dataset;

  herr_t status;

  char fname[100], set_name[100];

  int Hes, specs, j, k;

  float heating_rate[cooling_N_nH][cooling_N_Temp];
  float cooling_rate[cooling_N_nH][cooling_N_Temp];
  float temperature[cooling_N_nH][cooling_N_Temp];
  float electron_abundance[cooling_N_nH][cooling_N_Temp];


  sprintf(fname, "%sz_photodis.hdf5", All.CoolTablePath);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* For normal elements */
  for(specs = 0; specs < cooling_N_Elements; specs++)
    {
      sprintf(set_name, "/%s/Heating", cooling_ElementNames[specs]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, heating_rate);
      status = H5Dclose(dataset);

      sprintf(set_name, "/%s/Cooling", cooling_ElementNames[specs]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling_rate);
      status = H5Dclose(dataset);

      for(j = 0; j < cooling_N_Temp; j++)
	for(k = 0; k < cooling_N_nH; k++)
	  cooling_MetalsNetHeating[0][specs][k][j] = heating_rate[k][j] - cooling_rate[k][j];
    }

  /* Helium */
  for(Hes = 0; Hes < cooling_N_He; Hes++)
    {
      sprintf(set_name, "/Metal_free/%s/Heating", cooling_HeNames[Hes]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, heating_rate);
      status = H5Dclose(dataset);

      sprintf(set_name, "/Metal_free/%s/Cooling", cooling_HeNames[Hes]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling_rate);
      status = H5Dclose(dataset);

      sprintf(set_name, "/Metal_free/%s/Temperature", cooling_HeNames[Hes]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temperature);
      status = H5Dclose(dataset);

      sprintf(set_name, "/Metal_free/%s/Electron_density_over_n_h", cooling_HeNames[Hes]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, electron_abundance);
      status = H5Dclose(dataset);

      for(j = 0; j < cooling_N_Temp; j++)
	{
	  for(k = 0; k < cooling_N_nH; k++)
	    {
	      cooling_HplusHeNetHeating[0][Hes][k][j] = heating_rate[k][j] - cooling_rate[k][j];
	      cooling_ThermalToTemp[0][Hes][k][j] = log10(temperature[k][j]);
	      cooling_CollisionalElectronAbundance[Hes][k][j] = electron_abundance[k][j];
	    }
	}
    }

  status = H5Fclose(file_id);
}


/*
 * ----------------------------------------------------------------------
 * Get the cooling tables that bound the given redshift
 * ----------------------------------------------------------------------
 */

void GetCoolingTables(int highz_table_index, int lowz_table_index)
{
  hid_t file_id, dataset;

  herr_t status;

  char fname[100], set_name[100];

  int Hes, specs, i, j, skipme = 0;

  float heating_rate[cooling_N_nH][cooling_N_Temp];
  float cooling_rate[cooling_N_nH][cooling_N_Temp];
  float temperature[cooling_N_nH][cooling_N_Temp];

  static int loaded_highz_table_index = -1;
  static int loaded_lowz_table_index = -1;


  /* check to see if we can use the last one  */
  if(loaded_lowz_table_index == highz_table_index)
    {
      skipme = 1;
      for(specs = 0; specs < cooling_N_Elements; specs++)
	{
	  for(i = 0; i < cooling_N_nH; i++)
	    {
	      for(j = 0; j < cooling_N_Temp; j++)
		cooling_MetalsNetHeating[1][specs][i][j] =
		  cooling_MetalsNetHeating[0][specs][i][j];
	    }
	}

      for(specs = 0; specs < cooling_N_He; specs++)
	{
	  for(i = 0; i < cooling_N_nH; i++)
	    {
	      for(j = 0; j < cooling_N_Temp; j++)
		{
		  cooling_HplusHeNetHeating[1][specs][i][j] =
		    cooling_HplusHeNetHeating[0][specs][i][j];
		  cooling_ThermalToTemp[1][specs][i][j] =
		    cooling_ThermalToTemp[0][specs][i][j];
		}
	    }
	}
    }

  sprintf(fname, "%sz_%1.3f.hdf5", All.CoolTablePath, cooling_Redshifts[lowz_table_index]);

#ifdef BG_VERBOSE
  printf("Loading cooling table %s (actual redshift = %f)\n", fname, 1 / All.Time - 1);
  printf("First redshift, index = %d\n", lowz_table_index);
  fflush(stdout);
#endif

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  printf("Redshift 1 %ld %s\n", (long int) file_id, fname);
  fflush(stdout);

  if(file_id < 0)
    {
      printf("[GetCoolingTables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  /* Redshift number one */
  /* For normal elements */
  for(specs = 0; specs < cooling_N_Elements; specs++)
    {
      sprintf(set_name, "/%s/Heating", cooling_ElementNames[specs]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, heating_rate);
      status = H5Dclose(dataset);

      sprintf(set_name, "/%s/Cooling", cooling_ElementNames[specs]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling_rate);
      status = H5Dclose(dataset);

      for(i = 0; i < cooling_N_nH; i++)
	{
	  for(j = 0; j < cooling_N_Temp; j++)
	    cooling_MetalsNetHeating[0][specs][i][j] = heating_rate[i][j] - cooling_rate[i][j];
	}
    }

  /* Helium */
  for(Hes = 0; Hes < cooling_N_He; Hes++)
    {
      sprintf(set_name, "/Metal_free/%s/Heating", cooling_HeNames[Hes]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, heating_rate);
      status = H5Dclose(dataset);

      sprintf(set_name, "/Metal_free/%s/Cooling", cooling_HeNames[Hes]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling_rate);
      status = H5Dclose(dataset);

      sprintf(set_name, "/Metal_free/%s/Temperature", cooling_HeNames[Hes]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temperature);
      status = H5Dclose(dataset);

      for(i = 0; i < cooling_N_nH; i++)
	{
	  for(j = 0; j < cooling_N_Temp; j++)
	    {
	      cooling_HplusHeNetHeating[0][Hes][i][j] = heating_rate[i][j] - cooling_rate[i][j];
	      cooling_ThermalToTemp[0][Hes][i][j] = log10(temperature[i][j]);
	    }
	}
    }

  status = H5Fclose(file_id);

  /*redshift two */
  if(skipme == 0)
    {
      sprintf(fname, "%sz_%1.3f.hdf5", All.CoolTablePath, cooling_Redshifts[highz_table_index]);

#ifdef BG_VERBOSE
      printf("Loading cooling table %s (actual redshift = %f)\n", fname, 1 / All.Time - 1);
      printf("Second redshift, index = %d\n", highz_table_index);
      fflush(stdout);
#endif

      file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

      printf("Redshift 2 %ld %s\n", (long int) file_id, fname);
      fflush(stdout);

      if(file_id < 0)
	{
	  printf("[GetCoolingTables()]: unable to open file %s\n", fname);
	  endrun(101);
	}

      /* For normal elements */
      for(specs = 0; specs < cooling_N_Elements; specs++)
	{
	  sprintf(set_name, "/%s/Heating", cooling_ElementNames[specs]);
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, heating_rate);
	  status = H5Dclose(dataset);

	  sprintf(set_name, "/%s/Cooling", cooling_ElementNames[specs]);
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling_rate);
	  status = H5Dclose(dataset);

	  for(i = 0; i < cooling_N_nH; i++)
	    {
	      for(j = 0; j < cooling_N_Temp; j++)
		cooling_MetalsNetHeating[1][specs][i][j] = heating_rate[i][j] - cooling_rate[i][j];
	    }
	}

      /* Helium */
      for(Hes = 0; Hes < cooling_N_He; Hes++)
	{
	  sprintf(set_name, "/Metal_free/%s/Heating", cooling_HeNames[Hes]);
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, heating_rate);
	  status = H5Dclose(dataset);

	  sprintf(set_name, "/Metal_free/%s/Cooling", cooling_HeNames[Hes]);
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling_rate);
	  status = H5Dclose(dataset);

	  sprintf(set_name, "/Metal_free/%s/Temperature", cooling_HeNames[Hes]);
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temperature);
	  status = H5Dclose(dataset);

	  for(i = 0; i < cooling_N_nH; i++)
	    {
	      for(j = 0; j < cooling_N_Temp; j++)
		{
		  cooling_HplusHeNetHeating[1][Hes][i][j] = heating_rate[i][j] - cooling_rate[i][j];
		  cooling_ThermalToTemp[1][Hes][i][j] = log10(temperature[i][j]);
		}
	    }
	}

      status = H5Fclose(file_id);
    }

  loaded_highz_table_index = highz_table_index;
  loaded_lowz_table_index = lowz_table_index;
}


/*
 * ----------------------------------------------------------------------
 * This routine broadcasts the cooling tables and the electron abundance
 * for (mode = 0); only the cooling tables for (mode /= 0)
 * ----------------------------------------------------------------------
 */

void BroadcastCoolingTables(int mode)
{
  int position, bufsize, size, i, j, k;

  char *buffer;


  /* get size of the buffer  */
  MPI_Pack_size(2 * cooling_N_Temp * cooling_N_nH * cooling_N_Elements,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize = size;
  MPI_Pack_size(2 * cooling_N_Temp * cooling_N_nH * cooling_N_He,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(2 * cooling_N_Temp * cooling_N_nH * cooling_N_He,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;

  if(mode == 0)
    {
      MPI_Pack_size(cooling_N_Temp * cooling_N_nH * cooling_N_He,
		    MPI_FLOAT, MPI_COMM_WORLD, &size);
      bufsize += size;
    }

  /* allocate memory for the buffer  */
  buffer = (char *) mymalloc(bufsize);

  if(ThisTask == 0)
    {
      position = 0;
      for(i = 0; i < 2; i++)
	for(j = 0; j < cooling_N_Elements; j++)
	  for(k = 0; k < cooling_N_nH; k++)
	    MPI_Pack(cooling_MetalsNetHeating[i][j][k], cooling_N_Temp,
		     MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);

      for(i = 0; i < 2; i++)
	for(j = 0; j < cooling_N_He; j++)
	  for(k = 0; k < cooling_N_nH; k++)
	    MPI_Pack(cooling_HplusHeNetHeating[i][j][k], cooling_N_Temp,
		     MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);

      for(i = 0; i < 2; i++)
	for(j = 0; j < cooling_N_He; j++)
	  for(k = 0; k < cooling_N_nH; k++)
	    MPI_Pack(cooling_ThermalToTemp[i][j][k], cooling_N_Temp,
		     MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);

      if(mode == 0)
	for(j = 0; j < cooling_N_He; j++)
	  for(k = 0; k < cooling_N_nH; k++)
	    MPI_Pack(cooling_CollisionalElectronAbundance[j][k], cooling_N_Temp,
		     MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
    }

  MPI_Bcast(buffer, bufsize, MPI_PACKED, 0, MPI_COMM_WORLD);

  if(ThisTask != 0)
    {
      position = 0;
      for(i = 0; i < 2; i++)
	for(j = 0; j < cooling_N_Elements; j++)
	  for(k = 0; k < cooling_N_nH; k++)
	    MPI_Unpack(buffer, bufsize, &position, cooling_MetalsNetHeating[i][j][k],
		       cooling_N_Temp, MPI_FLOAT, MPI_COMM_WORLD);

      for(i = 0; i < 2; i++)
	for(j = 0; j < cooling_N_He; j++)
	  for(k = 0; k < cooling_N_nH; k++)
	    MPI_Unpack(buffer, bufsize, &position, cooling_HplusHeNetHeating[i][j][k],
		       cooling_N_Temp, MPI_FLOAT, MPI_COMM_WORLD);

      for(i = 0; i < 2; i++)
	for(j = 0; j < cooling_N_He; j++)
	  for(k = 0; k < cooling_N_nH; k++)
	    MPI_Unpack(buffer, bufsize, &position, cooling_ThermalToTemp[i][j][k],
		       cooling_N_Temp, MPI_FLOAT, MPI_COMM_WORLD);

      if(mode == 0)
	for(j = 0; j < cooling_N_He; j++)
	  for(k = 0; k < cooling_N_nH; k++)
	    MPI_Unpack(buffer, bufsize, &position, cooling_CollisionalElectronAbundance[j][k],
		       cooling_N_Temp, MPI_FLOAT, MPI_COMM_WORLD);
    }

  myfree(buffer);
}


/*
 * ----------------------------------------------------------------------
 * Check whether we need to load new tables
 * ----------------------------------------------------------------------
 */

void LoadCoolingTables(int z_index)
{
  static int highz_table_index = -1;
  static int lowz_table_index = -1;

  if(lowz_table_index == z_index)
    return;
  if(z_index >= cooling_N_Redshifts)		/* load collision only (high z) or high z reion table */
    {
      if(z_index == cooling_N_Redshifts + 1)
	{
	  if(ThisTask == 0)
	    GetNoComptTable();
	}
      else
	{
	  if(ThisTask == 0)
	    GetCollisTable();
	}
      BroadcastCoolingTables(0);
      highz_table_index = z_index;
      lowz_table_index = z_index;
    }
  else
    {
      lowz_table_index = z_index;
      highz_table_index = z_index + 1;

      if(ThisTask == 0)
	GetCoolingTables(highz_table_index, lowz_table_index);

      BroadcastCoolingTables(1);
    }
}

#endif
