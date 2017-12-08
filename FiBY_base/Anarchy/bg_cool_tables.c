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

  char fname[500], set_name[500];

  int specs, i, j, k;

  float net_cooling_rate[cooling_N_Temp][cooling_N_nH];
  float electron_abundance[cooling_N_Temp][cooling_N_nH];

  float temperature[cooling_N_He][cooling_N_Temp][cooling_N_nH];
  float he_net_cooling_rate[cooling_N_He][cooling_N_Temp][cooling_N_nH];
  float he_electron_abundance[cooling_N_He][cooling_N_Temp][cooling_N_nH];


  sprintf(fname, "%sz_8.989nocompton.hdf5", All.CoolTablePath);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  printf("Redshift 1 %ld %s\n", (long int) file_id, fname);
  fflush(stdout);

  /* For normal elements */
  for(specs = 0; specs < cooling_N_Elements; specs++)
    {
      sprintf(set_name, "/%s/Net_Cooling", cooling_ElementNames[specs]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, net_cooling_rate);
      status = H5Dclose(dataset);

      for(j = 0; j < cooling_N_Temp; j++)
	for(k = 0; k < cooling_N_nH; k++)
	  cooling_MetalsNetHeating[0][specs][k][j] = - net_cooling_rate[j][k];
    }

  /* Helium */
  sprintf(set_name, "/Metal_free/Net_Cooling");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, he_net_cooling_rate);
  status = H5Dclose(dataset);

  sprintf(set_name, "/Metal_free/Temperature/Temperature");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temperature);
  status = H5Dclose(dataset);

  sprintf(set_name, "/Metal_free/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, he_electron_abundance);
  status = H5Dclose(dataset);

  for(i = 0; i < cooling_N_He; i++)
    for(j = 0; j < cooling_N_Temp; j++)
      for(k = 0; k < cooling_N_nH; k++)
	{
	  cooling_HplusHeNetHeating[0][i][k][j] = - he_net_cooling_rate[i][j][k];
	  cooling_HplusHeElectronAbundance[0][i][k][j] = he_electron_abundance[i][j][k];
	  cooling_ThermalToTemp[0][i][k][j] = log10(temperature[i][j][k]);
	}

  sprintf(set_name, "/Solar/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, electron_abundance);
  status = H5Dclose(dataset);

  for(i = 0; i < cooling_N_Temp; i++)
    for(j = 0; j < cooling_N_nH; j++)
      cooling_SolarElectronAbundance[0][j][i] = electron_abundance[i][j];

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

  char fname[500], set_name[500];

  int specs, i, j, k;

  float net_cooling_rate[cooling_N_Temp][cooling_N_nH];
  float electron_abundance[cooling_N_Temp][cooling_N_nH];

  float temperature[cooling_N_He][cooling_N_Temp][cooling_N_nH];
  float he_net_cooling_rate[cooling_N_He][cooling_N_Temp][cooling_N_nH];
  float he_electron_abundance[cooling_N_He][cooling_N_Temp][cooling_N_nH];


  sprintf(fname, "%sz_photodis.hdf5", All.CoolTablePath);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* For normal elements */
  for(specs = 0; specs < cooling_N_Elements; specs++)
    {
      sprintf(set_name, "/%s/Net_Cooling", cooling_ElementNames[specs]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, net_cooling_rate);
      status = H5Dclose(dataset);

      for(j = 0; j < cooling_N_Temp; j++)
	for(k = 0; k < cooling_N_nH; k++)
	  cooling_MetalsNetHeating[0][specs][k][j] = - net_cooling_rate[j][k];
    }

  /* Helium */
  sprintf(set_name, "/Metal_free/Net_Cooling");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, he_net_cooling_rate);
  status = H5Dclose(dataset);

  sprintf(set_name, "/Metal_free/Temperature/Temperature");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temperature);
  status = H5Dclose(dataset);

  sprintf(set_name, "/Metal_free/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, he_electron_abundance);
  status = H5Dclose(dataset);

  for(i = 0; i < cooling_N_He; i++)
    for(j = 0; j < cooling_N_Temp; j++)
      for(k = 0; k < cooling_N_nH; k++)
	{
	  cooling_HplusHeNetHeating[0][i][k][j] = - he_net_cooling_rate[i][j][k];
	  cooling_HplusHeElectronAbundance[0][i][k][j] = he_electron_abundance[i][j][k];
	  cooling_ThermalToTemp[0][i][k][j] = log10(temperature[i][j][k]);
	}

  sprintf(set_name, "/Solar/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, electron_abundance);
  status = H5Dclose(dataset);

  for(i = 0; i < cooling_N_Temp; i++)
    for(j = 0; j < cooling_N_nH; j++)
      cooling_SolarElectronAbundance[0][j][i] = electron_abundance[i][j];

  status = H5Fclose(file_id);
}

#ifdef BG_COOLING_SHIELDING
/*
 * ----------------------------------------------------------------------
 * Get the cooling table for collisional cooling for self-shielding case
 * ----------------------------------------------------------------------
 */

void GetCollisTableShielding()
{
  hid_t file_id, dataset;

  herr_t status;

  char fname[500], set_name[500];

  int specs, i, j, k;

  float net_cooling_rate[cooling_N_Temp][cooling_N_nH];
  float electron_abundance[cooling_N_Temp][cooling_N_nH];

  float temperature[cooling_N_He][cooling_N_Temp][cooling_N_nH];
  float he_net_cooling_rate[cooling_N_He][cooling_N_Temp][cooling_N_nH];
  float he_electron_abundance[cooling_N_He][cooling_N_Temp][cooling_N_nH];

  if(ThisTask == 0)
    printf("Loading collisional tables for self-shielding approximation...");

  sprintf(fname, "%sz_photodis.hdf5", All.CoolTablePath);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* For normal elements */
  for(specs = 0; specs < cooling_N_Elements; specs++)
    {
      sprintf(set_name, "/%s/Net_Cooling", cooling_ElementNames[specs]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, net_cooling_rate);
      status = H5Dclose(dataset);

      for(j = 0; j < cooling_N_Temp; j++)
	for(k = 0; k < cooling_N_nH; k++)
	  cooling_MetalsNetHeatingCollisional[specs][k][j] = - net_cooling_rate[j][k];
    }

  /* Helium */
  sprintf(set_name, "/Metal_free/Net_Cooling");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, he_net_cooling_rate);
  status = H5Dclose(dataset);

  sprintf(set_name, "/Metal_free/Temperature/Temperature");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temperature);
  status = H5Dclose(dataset);

  sprintf(set_name, "/Metal_free/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, he_electron_abundance);
  status = H5Dclose(dataset);

  for(i = 0; i < cooling_N_He; i++)
    for(j = 0; j < cooling_N_Temp; j++)
      for(k = 0; k < cooling_N_nH; k++)
	{
	  cooling_HplusHeNetHeatingCollisional[i][k][j] = - he_net_cooling_rate[i][j][k];
	  cooling_HplusHeElectronAbundanceCollisional[i][k][j] = he_electron_abundance[i][j][k];
	  cooling_ThermalToTempCollisional[i][k][j] = log10(temperature[i][j][k]);
	}

  sprintf(set_name, "/Solar/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, electron_abundance);
  status = H5Dclose(dataset);

  for(i = 0; i < cooling_N_Temp; i++)
    for(j = 0; j < cooling_N_nH; j++)
      cooling_SolarElectronAbundanceCollisional[j][i] = electron_abundance[i][j];

  status = H5Fclose(file_id);

  if(ThisTask == 0)
    printf(" done.\n");
}
#endif

/*
 * ----------------------------------------------------------------------
 * Get the cooling tables that bound the given redshift
 * ----------------------------------------------------------------------
 */

void GetCoolingTables(int highz_table_index, int lowz_table_index)
{
  hid_t file_id, dataset;

  herr_t status;

  char fname[500], set_name[500];

  int specs, i, j, k, skipme = 0;

  float net_cooling_rate[cooling_N_Temp][cooling_N_nH];
  float electron_abundance[cooling_N_Temp][cooling_N_nH];

  float temperature[cooling_N_He][cooling_N_Temp][cooling_N_nH];
  float he_net_cooling_rate[cooling_N_He][cooling_N_Temp][cooling_N_nH];
  float he_electron_abundance[cooling_N_He][cooling_N_Temp][cooling_N_nH];

  static int loaded_highz_table_index = -1;
  static int loaded_lowz_table_index = -1;


  /* check to see if we can use the last one  */
  if(loaded_lowz_table_index == highz_table_index)
    {
      skipme = 1;

      for(specs = 0; specs < cooling_N_Elements; specs++)
	for(i = 0; i < cooling_N_nH; i++)
	  for(j = 0; j < cooling_N_Temp; j++)
	    cooling_MetalsNetHeating[1][specs][i][j] =
	      cooling_MetalsNetHeating[0][specs][i][j];

      for(specs = 0; specs < cooling_N_He; specs++)
	for(i = 0; i < cooling_N_nH; i++)
	  for(j = 0; j < cooling_N_Temp; j++)
	    {
	      cooling_HplusHeNetHeating[1][specs][i][j] =
		cooling_HplusHeNetHeating[0][specs][i][j];

	      cooling_HplusHeElectronAbundance[1][specs][i][j] =
		cooling_HplusHeElectronAbundance[0][specs][i][j];

	      cooling_ThermalToTemp[1][specs][i][j] =
		cooling_ThermalToTemp[0][specs][i][j];
	    }

      for(i = 0; i < cooling_N_nH; i++)
	for(j = 0; j < cooling_N_Temp; j++)
	  cooling_SolarElectronAbundance[1][i][j] =
	    cooling_SolarElectronAbundance[0][i][j];
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


  /*
   * Redshift number one
   */

  /* For normal elements */
  for(specs = 0; specs < cooling_N_Elements; specs++)
    {
      sprintf(set_name, "/%s/Net_Cooling", cooling_ElementNames[specs]);
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, net_cooling_rate);
      status = H5Dclose(dataset);

      for(i = 0; i < cooling_N_nH; i++)
	for(j = 0; j < cooling_N_Temp; j++)
	  cooling_MetalsNetHeating[0][specs][i][j] = - net_cooling_rate[j][i];
    }

  /* Helium */
  sprintf(set_name, "/Metal_free/Net_Cooling");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, he_net_cooling_rate);
  status = H5Dclose(dataset);

  sprintf(set_name, "/Metal_free/Temperature/Temperature");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temperature);
  status = H5Dclose(dataset);

  sprintf(set_name, "/Metal_free/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, he_electron_abundance);
  status = H5Dclose(dataset);

  for(i = 0; i < cooling_N_He; i++)
    for(j = 0; j < cooling_N_Temp; j++)
      for(k = 0; k < cooling_N_nH; k++)
	{
	  cooling_HplusHeNetHeating[0][i][k][j] = - he_net_cooling_rate[i][j][k];
	  cooling_HplusHeElectronAbundance[0][i][k][j] = he_electron_abundance[i][j][k];
	  cooling_ThermalToTemp[0][i][k][j] = log10(temperature[i][j][k]);
	}

  sprintf(set_name, "/Solar/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, electron_abundance);
  status = H5Dclose(dataset);

  for(i = 0; i < cooling_N_Temp; i++)
    for(j = 0; j < cooling_N_nH; j++)
      cooling_SolarElectronAbundance[0][j][i] = electron_abundance[i][j];

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
	  sprintf(set_name, "/%s/Net_Cooling", cooling_ElementNames[specs]);
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, net_cooling_rate);
	  status = H5Dclose(dataset);

	  for(i = 0; i < cooling_N_nH; i++)
	    for(j = 0; j < cooling_N_Temp; j++)
	      cooling_MetalsNetHeating[1][specs][i][j] = - net_cooling_rate[j][i];
	}

      /* Helium */
      sprintf(set_name, "/Metal_free/Net_Cooling");
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, he_net_cooling_rate);
      status = H5Dclose(dataset);

      sprintf(set_name, "/Metal_free/Temperature/Temperature");
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temperature);
      status = H5Dclose(dataset);

      sprintf(set_name, "/Metal_free/Electron_density_over_n_h");
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, he_electron_abundance);
      status = H5Dclose(dataset);

      for(i = 0; i < cooling_N_He; i++)
	for(j = 0; j < cooling_N_Temp; j++)
	  for(k = 0; k < cooling_N_nH; k++)
	    {
	      cooling_HplusHeNetHeating[1][i][k][j] = - he_net_cooling_rate[i][j][k];
	      cooling_HplusHeElectronAbundance[1][i][k][j] = he_electron_abundance[i][j][k];
	      cooling_ThermalToTemp[1][i][k][j] = log10(temperature[i][j][k]);
	    }

      sprintf(set_name, "/Solar/Electron_density_over_n_h");
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, electron_abundance);
      status = H5Dclose(dataset);

      for(i = 0; i < cooling_N_Temp; i++)
	for(j = 0; j < cooling_N_nH; j++)
	  cooling_SolarElectronAbundance[1][j][i] = electron_abundance[i][j];

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

/*
void BroadcastCoolingTables(int mode)
*/
void BroadcastCoolingTables(void)
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
  MPI_Pack_size(2 * cooling_N_Temp * cooling_N_nH * cooling_N_He,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(2 * cooling_N_Temp * cooling_N_nH,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;

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
	    MPI_Pack(cooling_HplusHeElectronAbundance[i][j][k], cooling_N_Temp,
		     MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);

      for(i = 0; i < 2; i++)
	for(j = 0; j < cooling_N_He; j++)
	  for(k = 0; k < cooling_N_nH; k++)
	    MPI_Pack(cooling_ThermalToTemp[i][j][k], cooling_N_Temp,
		     MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);

      for(i = 0; i < 2; i++)
	for(j = 0; j < cooling_N_nH; j++)
	  MPI_Pack(cooling_SolarElectronAbundance[i][j], cooling_N_Temp,
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
	    MPI_Unpack(buffer, bufsize, &position, cooling_HplusHeElectronAbundance[i][j][k],
		       cooling_N_Temp, MPI_FLOAT, MPI_COMM_WORLD);

      for(i = 0; i < 2; i++)
	for(j = 0; j < cooling_N_He; j++)
	  for(k = 0; k < cooling_N_nH; k++)
	    MPI_Unpack(buffer, bufsize, &position, cooling_ThermalToTemp[i][j][k],
		       cooling_N_Temp, MPI_FLOAT, MPI_COMM_WORLD);

      for(i = 0; i < 2; i++)
	for(j = 0; j < cooling_N_nH; j++)
	    MPI_Unpack(buffer, bufsize, &position, cooling_SolarElectronAbundance[i][j],
		       cooling_N_Temp, MPI_FLOAT, MPI_COMM_WORLD);
    }

  myfree(buffer);
}

#ifdef BG_COOLING_SHIELDING
void BroadcastCoolingTablesShielding(void)
{
  int position, bufsize, size, j, k;

  char *buffer;

  if(ThisTask == 0)
    printf("Broadcast collisional tables for self-shielding approximation...");

  /* get size of the buffer  */
  MPI_Pack_size(cooling_N_Temp * cooling_N_nH * cooling_N_Elements,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize = size;
  MPI_Pack_size(cooling_N_Temp * cooling_N_nH * cooling_N_He,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(cooling_N_Temp * cooling_N_nH * cooling_N_He,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(cooling_N_Temp * cooling_N_nH * cooling_N_He,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(cooling_N_Temp * cooling_N_nH,
		MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;

  /* allocate memory for the buffer  */
  buffer = (char *) mymalloc(bufsize);

  if(ThisTask == 0)
    {
      position = 0;

      for(j = 0; j < cooling_N_Elements; j++)
	for(k = 0; k < cooling_N_nH; k++)
	  MPI_Pack(cooling_MetalsNetHeatingCollisional[j][k], cooling_N_Temp,
		   MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);

      for(j = 0; j < cooling_N_He; j++)
	for(k = 0; k < cooling_N_nH; k++)
	  MPI_Pack(cooling_HplusHeNetHeatingCollisional[j][k], cooling_N_Temp,
		   MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);

      for(j = 0; j < cooling_N_He; j++)
	for(k = 0; k < cooling_N_nH; k++)
	  MPI_Pack(cooling_HplusHeElectronAbundanceCollisional[j][k], cooling_N_Temp,
		   MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);

      for(j = 0; j < cooling_N_He; j++)
	for(k = 0; k < cooling_N_nH; k++)
	  MPI_Pack(cooling_ThermalToTempCollisional[j][k], cooling_N_Temp,
		   MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);

      for(j = 0; j < cooling_N_nH; j++)
	MPI_Pack(cooling_SolarElectronAbundanceCollisional[j], cooling_N_Temp,
		 MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
    }

  MPI_Bcast(buffer, bufsize, MPI_PACKED, 0, MPI_COMM_WORLD);

  if(ThisTask != 0)
    {
      position = 0;

      for(j = 0; j < cooling_N_Elements; j++)
	for(k = 0; k < cooling_N_nH; k++)
	  MPI_Unpack(buffer, bufsize, &position, cooling_MetalsNetHeatingCollisional[j][k],
		     cooling_N_Temp, MPI_FLOAT, MPI_COMM_WORLD);

      for(j = 0; j < cooling_N_He; j++)
	for(k = 0; k < cooling_N_nH; k++)
	  MPI_Unpack(buffer, bufsize, &position, cooling_HplusHeNetHeatingCollisional[j][k],
		     cooling_N_Temp, MPI_FLOAT, MPI_COMM_WORLD);

      for(j = 0; j < cooling_N_He; j++)
	for(k = 0; k < cooling_N_nH; k++)
	  MPI_Unpack(buffer, bufsize, &position, cooling_HplusHeElectronAbundanceCollisional[j][k],
		     cooling_N_Temp, MPI_FLOAT, MPI_COMM_WORLD);

      for(j = 0; j < cooling_N_He; j++)
	for(k = 0; k < cooling_N_nH; k++)
	  MPI_Unpack(buffer, bufsize, &position, cooling_ThermalToTempCollisional[j][k],
		     cooling_N_Temp, MPI_FLOAT, MPI_COMM_WORLD);

      for(j = 0; j < cooling_N_nH; j++)
	MPI_Unpack(buffer, bufsize, &position, cooling_SolarElectronAbundanceCollisional[j],
		   cooling_N_Temp, MPI_FLOAT, MPI_COMM_WORLD);
    }

  myfree(buffer);

  if(ThisTask == 0)
    printf(" done\n");
}
#endif

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
      BroadcastCoolingTables();
      highz_table_index = z_index;
      lowz_table_index = z_index;
    }
  else
    {
      lowz_table_index = z_index;
      highz_table_index = z_index + 1;

      if(ThisTask == 0)
	GetCoolingTables(highz_table_index, lowz_table_index);

      BroadcastCoolingTables();
    }
}

#ifdef BG_COOLING_SHIELDING
void LoadCollisionalCoolingTables()
{
  GetCollisTableShielding();
  BroadcastCoolingTablesShielding();
}
#endif

#endif
