#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#include "bg_proto.h"
#include "bg_yields.h"

#ifdef BG_STELLAR_EVOLUTION

double *yield_mass_bin, *stellar_yield, *integrand;
double *imf_by_number, *imf_mass_bin, *imf_mass_bin_log10;

#ifdef BG_POPIII
double *popiii_yield_mass_bin, *popiii_stellar_yield, *popiii_integrand;
double *popiii_imf_by_number, *popiii_imf_mass_bin, *popiii_imf_mass_bin_log10;
#endif

#ifdef BG_DOUBLE_IMF
double *imf_by_number1;
#endif

struct YieldTypeIa yieldsSNIa;

struct YieldTypeII_and_AGB yieldsSNII, yieldsAGB;

#ifdef BG_DUST
struct DustTypeII_and_AGB dustSNII, dustAGB;
#endif

struct Lifetime_Table Lifetimes;

#ifdef BG_POPIII
struct YieldPOPIII yieldsPOPIII;

struct POPIIILifetime_Table POPIIILifetimes;

#ifdef BG_DUST
struct DustPOPIII dustPOPIII;

#define DUST_YIELD_Z_NAME_LENGTH 20
#endif

#define POPIII_YIELD_EL_NAME_LENGTH 20
#endif /* BG_POPIII */

/*
 * ----------------------------------------------------------------------
 * This routine reads in the yield tables (HDF5 format). It is called by
 * processor 0. Processor 0 will then broadcast the tables.
 * ----------------------------------------------------------------------
 */

int read_yield_tables()
{
  int i, j, k;

  char fname[256], setname[100];

  hid_t file_id, dataset, datatype;

  /* --------------------------------- */
  /* read in table dimensions for SNIa */
  /* --------------------------------- */

  sprintf(fname, "%s/%s", All.YieldTablePath, SNIa_yieldname);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(file_id, "Number_of_species");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsSNIa.N_ELEMENTS);
  H5Dclose(dataset);

  H5Fclose(file_id);

  /* --------------------------------- */
  /* read in table dimensions for SNII */
  /* --------------------------------- */

  sprintf(fname, "%s/%s", All.YieldTablePath, SNII_yieldname);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(file_id, "Number_of_species");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsSNII.N_ELEMENTS);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "Number_of_masses");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsSNII.N_MASS);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "Number_of_metallicities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsSNII.N_Z);
  H5Dclose(dataset);

  H5Fclose(file_id);

#ifdef BG_DUST
  sprintf(fname, "%s/%s", All.YieldTablePath, SNII_dustname);

  printf("%s\n", fname);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(file_id, "Number_of_masses");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dustSNII.N_MASS);
  H5Dclose(dataset);

  printf("SNII n_Mass = %d\n", dustSNII.N_MASS);

  dataset = H5Dopen(file_id, "Number_of_metallicities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dustSNII.N_Z);
  H5Dclose(dataset);

  printf("SNII n_Z = %d\n", dustSNII.N_Z);

  H5Fclose(file_id);
#endif

  /* -------------------------------- */
  /* read in table dimensions for AGB */
  /* -------------------------------- */

  sprintf(fname, "%s/%s", All.YieldTablePath, AGB_yieldname);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(file_id, "Number_of_species");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsAGB.N_ELEMENTS);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "Number_of_masses");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsAGB.N_MASS);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "Number_of_metallicities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsAGB.N_Z);
  H5Dclose(dataset);

  H5Fclose(file_id);

#ifdef BG_DUST
  sprintf(fname, "%s/%s", All.YieldTablePath, AGB_dustname);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(file_id, "Number_of_masses");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dustAGB.N_MASS);
  H5Dclose(dataset);

  printf("AGB n_Mass = %d\n", dustAGB.N_MASS);

  dataset = H5Dopen(file_id, "Number_of_metallicities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dustAGB.N_Z);
  H5Dclose(dataset);

  printf("AGB n_Z = %d\n", dustAGB.N_Z);

  H5Fclose(file_id);
#endif

  /* ----------------------------------- */
  /* read in table dimensions for POPIII */
  /* ----------------------------------- */
#ifdef BG_POPIII
  sprintf(fname, "%s/%s", All.YieldTablePath, POPIII_yieldname);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(file_id, "Number_of_species");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsPOPIII.N_ELEMENTS);
  H5Dclose(dataset);

  /* printf("[POPIII] Number of elements = %d\n", yieldsPOPIII.N_ELEMENTS); */

  dataset = H5Dopen(file_id, "Number_of_masses");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsPOPIII.N_MASS);
  H5Dclose(dataset);

  /* printf("[POPIII] Number of masses = %d\n", yieldsPOPIII.N_MASS); */

  H5Fclose(file_id);

#ifdef BG_DUST
  sprintf(fname, "%s/%s", All.YieldTablePath, POPIII_dustname);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(file_id, "Number_of_masses");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dustPOPIII.N_MASS);
  H5Dclose(dataset);

  printf("POPIII n_Mass = %d\n", dustPOPIII.N_MASS);

  H5Fclose(file_id);
#endif
#endif

  /* ------------------------------------------- */
  /* read in table dimensions for Lifetime table */
  /* ------------------------------------------- */

  sprintf(fname, "%s/%s", All.YieldTablePath, Lifetime_name);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(file_id, "Number_of_masses");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lifetimes.N_MASS);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "Number_of_metallicities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lifetimes.N_Z);
  H5Dclose(dataset);

  H5Fclose(file_id);

#ifdef BG_POPIII
  sprintf(fname, "%s/%s", All.YieldTablePath, POPIIILifetime_name);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(file_id, "Number_of_masses");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &POPIIILifetimes.N_MASS);
  H5Dclose(dataset);

  H5Fclose(file_id);
#endif

  /* -------------------------------- */
  /* allocate memory for yield tables */
  /* -------------------------------- */

  allocate_yield_tables();

  /* ------------------------ */
  /* read in SNIa yield table */
  /* ------------------------ */

  sprintf(fname, "%s/%s", All.YieldTablePath, SNIa_yieldname);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  /* allocate element name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names");
  H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, yieldsSNIa.ElementName);
  H5Dclose(dataset);
  H5Tclose(datatype);

  for(i = 0; i < yieldsSNIa.N_ELEMENTS; i++)
    yieldsSNIa.ElementName[i] = mystrdup(yieldsSNIa.ElementName[i]);

  dataset = H5Dopen(file_id, "Yield");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, yieldsSNIa.Yield);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "Total_Metals");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsSNIa.TotalMetals_SPH);
  H5Dclose(dataset);

  H5Fclose(file_id);

  /* -------------------------------------------------------- */
  /* read in SNII mass and metallicity arrays and yield table */
  /* -------------------------------------------------------- */

  sprintf(fname, "%s/%s", All.YieldTablePath, SNII_yieldname);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  /* allocate element name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names");
  H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, yieldsSNII.ElementName);
  H5Dclose(dataset);
  H5Tclose(datatype);

  for(i = 0; i < yieldsSNII.N_ELEMENTS; i++)
    yieldsSNII.ElementName[i] = mystrdup(yieldsSNII.ElementName[i]);

  dataset = H5Dopen(file_id, "Masses");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, yieldsSNII.Mass);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "Metallicities");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, yieldsSNII.Metallicity);
  H5Dclose(dataset);

  double tempyield1[yieldsSNII.N_ELEMENTS][yieldsSNII.N_MASS];

  double tempej1[yieldsSNII.N_MASS], tempmet1[yieldsSNII.N_MASS];

  char *tempname1[yieldsSNII.N_Z];

  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Yield_names");
  H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempname1);
  H5Dclose(dataset);
  H5Tclose(datatype);

  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      sprintf(setname, "/Yields/%s/Yield", tempname1[i]);
      dataset = H5Dopen(file_id, setname);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempyield1);
      H5Dclose(dataset);
      sprintf(setname, "/Yields/%s/Ejected_mass", tempname1[i]);
      dataset = H5Dopen(file_id, setname);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempej1);
      H5Dclose(dataset);
      sprintf(setname, "/Yields/%s/Total_Metals", tempname1[i]);
      dataset = H5Dopen(file_id, setname);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempmet1);
      H5Dclose(dataset);

      for(k = 0; k < yieldsSNII.N_MASS; k++)
	{
	  yieldsSNII.Ejecta[i][k] = tempej1[k];
	  yieldsSNII.TotalMetals[i][k] = tempmet1[k];

	  for(j = 0; j < yieldsSNII.N_ELEMENTS; j++)
	    yieldsSNII.Yield[i][j][k] = tempyield1[j][k];
	}
    }

  H5Fclose(file_id);

#ifdef BG_DUST
  sprintf(fname, "%s/%s", All.YieldTablePath, SNII_dustname);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(file_id, "Masses");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dustSNII.Mass);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "Metallicities");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dustSNII.Metallicity);
  H5Dclose(dataset);

  printf("Before yields\n");

  /*                                                                                                                                                                                                                             
    Dust yield table was made with IDL which does not fully support HDF5                                                                                                                                                         
    definition of strings of variable length. This trick was also applied                                                                                                                                                        
    to the reading of the new cooling tables.                                                                                                                                                                                    
  */

  char z_names[dustSNII.N_Z][DUST_YIELD_Z_NAME_LENGTH];
  hsize_t z_string_length = DUST_YIELD_Z_NAME_LENGTH;

  double tempdust[dustSNII.N_MASS];
  char *tempzname[dustSNII.N_Z];

  /* read yield name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, z_string_length);
  dataset = H5Dopen(file_id, "Yield_names");
  H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, z_names);
  H5Dclose(dataset);
  H5Tclose(datatype);

  for(i = 0; i < dustSNII.N_Z; i++)
    tempzname[i] = mystrdup(z_names[i]);

  for(i = 0; i < dustSNII.N_Z; i++)
    {
      printf("%s\n", tempzname[i]);

      sprintf(setname, "/Yields/%s/Total_Dust", tempzname[i]);
      dataset = H5Dopen(file_id, setname);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempdust);
      H5Dclose(dataset);

      for(k = 0; k < dustSNII.N_MASS; k++)
        dustSNII.TotalDust[i][k] = tempdust[k];
    }


  printf("After yields\n");

  H5Fclose(file_id);
#endif

  /* ------------------------------------------------------- */
  /* read in AGB mass and metallicity arrays and yield table */
  /* ------------------------------------------------------- */

  sprintf(fname, "%s/%s", All.YieldTablePath, AGB_yieldname);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  /* allocate element name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names");
  H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, yieldsAGB.ElementName);
  H5Dclose(dataset);
  H5Tclose(datatype);

  for(i = 0; i < yieldsAGB.N_ELEMENTS; i++)
    yieldsAGB.ElementName[i] = mystrdup(yieldsAGB.ElementName[i]);

  dataset = H5Dopen(file_id, "Masses");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, yieldsAGB.Mass);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "Metallicities");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, yieldsAGB.Metallicity);
  H5Dclose(dataset);

  double tempyield2[yieldsAGB.N_ELEMENTS][yieldsAGB.N_MASS];

  double tempej2[yieldsAGB.N_MASS], tempmet2[yieldsAGB.N_MASS];

  char *tempname2[yieldsAGB.N_Z];

  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Yield_names");
  H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempname2);
  H5Dclose(dataset);
  H5Tclose(datatype);

  for(i = 0; i < yieldsAGB.N_Z; i++)
    {
      sprintf(setname, "/Yields/%s/Yield", tempname2[i]);
      dataset = H5Dopen(file_id, setname);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempyield2);
      H5Dclose(dataset);
      sprintf(setname, "/Yields/%s/Ejected_mass", tempname2[i]);
      dataset = H5Dopen(file_id, setname);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempej2);
      H5Dclose(dataset);
      sprintf(setname, "/Yields/%s/Total_Metals", tempname2[i]);
      dataset = H5Dopen(file_id, setname);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempmet2);
      H5Dclose(dataset);

      for(k = 0; k < yieldsAGB.N_MASS; k++)
	{
	  yieldsAGB.Ejecta[i][k] = tempej2[k];
	  yieldsAGB.TotalMetals[i][k] = tempmet2[k];

	  for(j = 0; j < yieldsAGB.N_ELEMENTS; j++)
	    yieldsAGB.Yield[i][j][k] = tempyield2[j][k];
	}
    }

  H5Fclose(file_id);

#ifdef BG_DUST
  sprintf(fname, "%s/%s", All.YieldTablePath, AGB_dustname);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(file_id, "Masses");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dustAGB.Mass);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "Metallicities");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dustAGB.Metallicity);
  H5Dclose(dataset);

  /*                                                                                                                                                                                                                             
    Dust yield table was made with IDL which does not fully support HDF5                                                                                                                                                         
    definition of strings of variable length. This trick was also applied                                                                                                                                                        
    to the reading of the new cooling tables.                                                                                                                                                                                    
  */

  char z_names1[dustAGB.N_Z][DUST_YIELD_Z_NAME_LENGTH];
  hsize_t z_string_length1 = DUST_YIELD_Z_NAME_LENGTH;

  double tempdust1[dustAGB.N_MASS];
  char *tempzname1[dustAGB.N_Z];

  /* read yield name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, z_string_length1);
  dataset = H5Dopen(file_id, "Yield_names");
  H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, z_names1);
  H5Dclose(dataset);
  H5Tclose(datatype);

  for(i = 0; i < dustAGB.N_Z; i++)
    tempzname1[i] = mystrdup(z_names1[i]);

  for(i = 0; i < dustAGB.N_Z; i++)
    {
      printf("%s\n", tempzname1[i]);

      sprintf(setname, "/Yields/%s/Total_Dust", tempzname1[i]);

      dataset = H5Dopen(file_id, setname);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempdust1);
      H5Dclose(dataset);

      for(k = 0; k < dustAGB.N_MASS; k++)
        dustAGB.TotalDust[i][k] = tempdust1[k];
    }

  H5Fclose(file_id);
#endif

  /* ---------------------------------------------------------- */
  /* read in POPIII mass and metallicity arrays and yield table */
  /* ---------------------------------------------------------- */
#ifdef BG_POPIII
  sprintf(fname, "%s/%s", All.YieldTablePath, POPIII_yieldname);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  /*
  if(ThisTask == 0) printf("[POPIII] el. names\n");
  fflush(stdout);
  */

  /*
    POPIII yield table was made with IDL which does not fully support HDF5
    definition of strings of variable length. This trick was also applied
    to the reading of the new cooling tables.
  */

  char element_names[yieldsPOPIII.N_ELEMENTS][POPIII_YIELD_EL_NAME_LENGTH];
  hsize_t string_length = POPIII_YIELD_EL_NAME_LENGTH;

  /* allocate element name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, string_length);
  dataset = H5Dopen(file_id, "Species_names");
  H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, element_names);
  H5Dclose(dataset);
  H5Tclose(datatype);

  for(i = 0; i < yieldsPOPIII.N_ELEMENTS; i++)
    yieldsPOPIII.ElementName[i] = mystrdup(element_names[i]);

  /*
  for(i = 0; i < yieldsPOPIII.N_ELEMENTS; i++)
    printf("[POPIII] element_name[%d] = %s\n", i, yieldsPOPIII.ElementName[i]);
  fflush(stdout);

  if(ThisTask == 0) printf("[POPIII] masses\n");
  fflush(stdout);
  */

  dataset = H5Dopen(file_id, "Masses");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, yieldsPOPIII.Mass);
  H5Dclose(dataset);

  /*
  for(i = 0; i < yieldsPOPIII.N_MASS; i++)
    printf("[POPIII] mass[%d] = %f\n", i, yieldsPOPIII.Mass[i]);
  fflush(stdout);

  if(ThisTask == 0) printf("[POPIII] ejecta\n");
  fflush(stdout);
  */

  dataset = H5Dopen(file_id, "Ejected_mass");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, yieldsPOPIII.Ejecta);
  H5Dclose(dataset);

  /*
  for(i = 0; i < yieldsPOPIII.N_MASS; i++)
    printf("[POPIII] ejected_mass[%d] = %f\n", i, yieldsPOPIII.Ejecta[i]);
  fflush(stdout);

  if(ThisTask == 0) printf("[POPIII] total metals\n");
  fflush(stdout);
  */

  dataset = H5Dopen(file_id, "Total_Metals");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, yieldsPOPIII.TotalMetals);
  H5Dclose(dataset);

  /*
  for(i = 0; i < yieldsPOPIII.N_MASS; i++)
    printf("[POPIII] total_metals[%d] = %f\n", i, yieldsPOPIII.TotalMetals[i]);
  fflush(stdout);

  if(ThisTask == 0) printf("[POPIII] yields\n");
  fflush(stdout);
  */

  /*
    Again IDL does stupid things. It saves the 2-dimensional yield array
    swapping the dimensions. C reads it in right though.
  */

  double tempyield3[yieldsPOPIII.N_MASS][yieldsPOPIII.N_ELEMENTS];

  dataset = H5Dopen(file_id, "Yield");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempyield3);
  H5Dclose(dataset);

  /*
  for(j = 0; j < yieldsPOPIII.N_ELEMENTS; j++)
    printf("[POPIII] tmp_yield[0][%d] = %e\n", j, tempyield3[0][j]);
  fflush(stdout);

  for(j = 0; j < yieldsPOPIII.N_ELEMENTS; j++)
    printf("[POPIII] tmp_yield[%d][%d] = %e\n", yieldsPOPIII.N_MASS-1, j, tempyield3[yieldsPOPIII.N_MASS-1][j]);
  fflush(stdout);
  */

  for(k = 0; k < yieldsPOPIII.N_ELEMENTS; k++)
    for(j = 0; j < yieldsPOPIII.N_MASS; j++)
      yieldsPOPIII.Yield[k][j] = tempyield3[j][k];

  /*
  for(j = 0; j < yieldsPOPIII.N_ELEMENTS; j++)
    printf("[POPIII] yield[%d][0] = %e\n", j, yieldsPOPIII.Yield[j][0]);
  fflush(stdout);

  for(j = 0; j < yieldsPOPIII.N_ELEMENTS; j++)
    printf("[POPIII] yield[%d][%d] = %e\n", j, yieldsPOPIII.N_MASS-1, yieldsPOPIII.Yield[j][yieldsPOPIII.N_MASS-1]);
  fflush(stdout);
  */

  H5Fclose(file_id);

#ifdef BG_DUST
  sprintf(fname, "%s/%s", All.YieldTablePath, POPIII_dustname);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(file_id, "Masses");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dustPOPIII.Mass);
H5Dclose(dataset);

for(i = 0; i < dustPOPIII.N_MASS; i++)
    printf("[POPIII] mass[%d] = %f\n", i, dustPOPIII.Mass[i]);
  fflush(stdout);

dataset = H5Dopen(file_id, "Total_Dust");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dustPOPIII.TotalDust);
  H5Dclose(dataset);

  for(i = 0; i < dustPOPIII.N_MASS; i++)
    printf("[POPIII] total_dust[%d] = %f\n", i, dustPOPIII.TotalDust[i]);
  fflush(stdout);

  H5Fclose(file_id);
#endif
#endif

  /* ---------------------- */
  /* read in Lifetime table */
  /* ---------------------- */

  sprintf(fname, "%s/%s", All.YieldTablePath, Lifetime_name);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(file_id, "Masses");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Lifetimes.Mass);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "Metallicities");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Lifetimes.Metallicity);
  H5Dclose(dataset);

  double temptime[Lifetimes.N_Z][Lifetimes.N_MASS];

  dataset = H5Dopen(file_id, "Lifetimes");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temptime);
  H5Dclose(dataset);

  for(i = 0; i < Lifetimes.N_Z; i++)
    for(j = 0; j < Lifetimes.N_MASS; j++)
      Lifetimes.Dyingtime[i][j] = log10(temptime[i][j]);

  H5Fclose(file_id);

#ifdef BG_POPIII
  sprintf(fname, "%s/%s", All.YieldTablePath, POPIIILifetime_name);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("[read_yield_tables()]: unable to open file %s\n", fname);
      endrun(101);
    }

  dataset = H5Dopen(file_id, "Masses");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, POPIIILifetimes.Mass);
  H5Dclose(dataset);

  /*
  for(i = 0; i < POPIIILifetimes.N_MASS; i++)
    printf("%f\n", POPIIILifetimes.Mass[i]);
  */

  dataset = H5Dopen(file_id, "Lifetimes");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, POPIIILifetimes.Dyingtime);
  H5Dclose(dataset);

  /*
  for(i = 0; i < POPIIILifetimes.N_MASS; i++)
    printf("%f\n", POPIIILifetimes.Dyingtime[i]);
  */

  for(i = 0; i < POPIIILifetimes.N_MASS; i++)
   POPIIILifetimes.Dyingtime[i] = log10(POPIIILifetimes.Dyingtime[i]);

  /*
  for(i = 0; i < POPIIILifetimes.N_MASS; i++)
    printf("%f\n", POPIIILifetimes.Dyingtime[i]);
  */

  H5Fclose(file_id);
#endif

  return 0;
}


/*
 * ----------------------------------------------------------------------
 * This routine broadcast the yield table dimensions.
 * ----------------------------------------------------------------------
 */

int bcast_yield_table_dim()
{
  int position, bufsize, packsize;

  char *buffer;

  /* --------------------------------------------------- */
  /* broadcast SNIa, SNII and AGB yield table dimensions */
  /* --------------------------------------------------- */

  /* get size of the buffer */
  /*
#ifdef BG_POPIII
  MPI_Pack_size(12, MPI_INT, MPI_COMM_WORLD, &bufsize);
#else
  MPI_Pack_size(9, MPI_INT, MPI_COMM_WORLD, &bufsize);
#endif
  */

  /* get size of the buffer */
  packsize = 9;
#ifdef BG_DUST
  packsize += 4; /* for dust from SNII and AGB, N_MASS and N_Z */
#endif
#ifdef BG_POPIII
  packsize += 3;
#ifdef BG_DUST
  packsize += 1; /* for dust from POPIII stars, N_MASS */
#endif
#endif

  MPI_Pack_size(packsize, MPI_INT, MPI_COMM_WORLD, &bufsize);

  /* allocate memory for the buffer */
  buffer = (char *) mymalloc(bufsize);

  if(ThisTask == 0)
    {
      position = 0;

      /* pack SNIa table dimension */
      MPI_Pack(&yieldsSNIa.N_ELEMENTS, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);

      /* pack SNII table dimensions */
      MPI_Pack(&yieldsSNII.N_ELEMENTS, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(&yieldsSNII.N_MASS, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(&yieldsSNII.N_Z, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
#ifdef BG_DUST
      MPI_Pack(&dustSNII.N_MASS, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(&dustSNII.N_Z, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
#endif

      /* pack AGB stars table dimensions */
      MPI_Pack(&yieldsAGB.N_ELEMENTS, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(&yieldsAGB.N_MASS, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(&yieldsAGB.N_Z, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
#ifdef BG_DUST
      MPI_Pack(&dustAGB.N_MASS, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(&dustAGB.N_Z, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
#endif

      /* pack POPIII stars table dimensions */
#ifdef BG_POPIII
      MPI_Pack(&yieldsPOPIII.N_ELEMENTS, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(&yieldsPOPIII.N_MASS, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
#ifdef BG_DUST
      MPI_Pack(&dustPOPIII.N_MASS, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
#endif
#endif
      /* pack Lifetime table dimensions */
      MPI_Pack(&Lifetimes.N_MASS, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(&Lifetimes.N_Z, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
#ifdef BG_POPIII
      MPI_Pack(&POPIIILifetimes.N_MASS, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
#endif
    }

  MPI_Bcast(buffer, bufsize, MPI_PACKED, 0, MPI_COMM_WORLD);

  if(ThisTask != 0)
    {
      position = 0;

      /* unpack SNIa table dimension */
      MPI_Unpack(buffer, bufsize, &position, &yieldsSNIa.N_ELEMENTS, 1, MPI_INT, MPI_COMM_WORLD);

      /* unpack SNII table dimensions */
      MPI_Unpack(buffer, bufsize, &position, &yieldsSNII.N_ELEMENTS, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, &yieldsSNII.N_MASS, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, &yieldsSNII.N_Z, 1, MPI_INT, MPI_COMM_WORLD);
#ifdef BG_DUST
      MPI_Unpack(buffer, bufsize, &position, &dustSNII.N_MASS, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, &dustSNII.N_Z, 1, MPI_INT, MPI_COMM_WORLD);
#endif

      /* unpack AGB stars table dimensions */
      MPI_Unpack(buffer, bufsize, &position, &yieldsAGB.N_ELEMENTS, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, &yieldsAGB.N_MASS, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, &yieldsAGB.N_Z, 1, MPI_INT, MPI_COMM_WORLD);
#ifdef BG_DUST
      MPI_Unpack(buffer, bufsize, &position, &dustAGB.N_MASS, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, &dustAGB.N_Z, 1, MPI_INT, MPI_COMM_WORLD);
#endif

      /* unpack POPIII stars table dimensions */
#ifdef BG_POPIII
      MPI_Unpack(buffer, bufsize, &position, &yieldsPOPIII.N_ELEMENTS, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, &yieldsPOPIII.N_MASS, 1, MPI_INT, MPI_COMM_WORLD);
#ifdef BG_DUST
      MPI_Unpack(buffer, bufsize, &position, &dustPOPIII.N_MASS, 1, MPI_INT, MPI_COMM_WORLD);
#endif
#endif
      /* unpack Lifetime table dimensions */
      MPI_Unpack(buffer, bufsize, &position, &Lifetimes.N_MASS, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, &Lifetimes.N_Z, 1, MPI_INT, MPI_COMM_WORLD);
#ifdef BG_POPIII
      MPI_Unpack(buffer, bufsize, &position, &POPIIILifetimes.N_MASS, 1, MPI_INT, MPI_COMM_WORLD);
#endif
    }

  myfree(buffer);

  return 0;
}


/*
 * ----------------------------------------------------------------------
 * This routine broadcast the yield table values
 * ----------------------------------------------------------------------
 */

int bcast_yield_tables()
{
  int i, j, position, bufsize, size;

  char *buffer;

  /* --------------------------- */
  /* broadcast SNIa yield tables */
  /* --------------------------- */

  bufsize = 0;

  /* buffer size for SNIa */
  MPI_Pack_size(yieldsSNIa.N_ELEMENTS * EL_NAME_LENGTH, MPI_CHAR, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsSNIa.N_ELEMENTS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;

  /* buffer size for SNII */
  MPI_Pack_size(yieldsSNII.N_ELEMENTS * EL_NAME_LENGTH, MPI_CHAR, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsSNII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsSNII.N_Z, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsSNII.N_Z * yieldsSNII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsSNII.N_Z * yieldsSNII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsSNII.N_ELEMENTS * yieldsSNII.N_Z * yieldsSNII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD,
		&size);
  bufsize += size;
#ifdef BG_DUST
  MPI_Pack_size(dustSNII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(dustSNII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(dustSNII.N_Z * dustSNII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
#endif

  /* buffer size for AGB stars */
  MPI_Pack_size(yieldsAGB.N_ELEMENTS * EL_NAME_LENGTH, MPI_CHAR, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsAGB.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsAGB.N_Z, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsAGB.N_Z * yieldsAGB.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsAGB.N_Z * yieldsAGB.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsAGB.N_ELEMENTS * yieldsAGB.N_Z * yieldsAGB.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
#ifdef BG_DUST
  MPI_Pack_size(dustAGB.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(dustAGB.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(dustAGB.N_Z * dustAGB.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
#endif

  /* buffer size POPIII */
#ifdef BG_POPIII
  MPI_Pack_size(yieldsPOPIII.N_ELEMENTS * POPIII_YIELD_EL_NAME_LENGTH, MPI_CHAR, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsPOPIII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsPOPIII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsPOPIII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(yieldsPOPIII.N_ELEMENTS * yieldsPOPIII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
#ifdef BG_DUST
  MPI_Pack_size(dustPOPIII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(dustPOPIII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
#endif
#endif
  /* buffer size for Lifetimes */
  MPI_Pack_size(Lifetimes.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(Lifetimes.N_Z, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(Lifetimes.N_Z * Lifetimes.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
#ifdef BG_POPIII
  MPI_Pack_size(POPIIILifetimes.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(POPIIILifetimes.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  bufsize += size;
#endif

  buffer = (char *) mymalloc(bufsize);

  if(ThisTask == 0)
    {
      position = 0;

      /* pack SNIa table */
      for(i = 0; i < yieldsSNIa.N_ELEMENTS; i++)
	MPI_Pack(yieldsSNIa.ElementName[i], EL_NAME_LENGTH, MPI_CHAR, buffer, bufsize, &position,
		 MPI_COMM_WORLD);
      MPI_Pack(yieldsSNIa.Yield, yieldsSNIa.N_ELEMENTS, MPI_DOUBLE, buffer, bufsize, &position,
	       MPI_COMM_WORLD);
      MPI_Pack(&yieldsSNIa.TotalMetals_SPH, 1, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);

      /* pack SNII table */
      for(i = 0; i < yieldsSNII.N_ELEMENTS; i++)
	MPI_Pack(yieldsSNII.ElementName[i], EL_NAME_LENGTH, MPI_CHAR, buffer, bufsize, &position,
		 MPI_COMM_WORLD);
      MPI_Pack(yieldsSNII.Mass, yieldsSNII.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(yieldsSNII.Metallicity, yieldsSNII.N_Z, MPI_DOUBLE, buffer, bufsize, &position,
	       MPI_COMM_WORLD);
      for(i = 0; i < yieldsSNII.N_Z; i++)
	{
	  MPI_Pack(yieldsSNII.Ejecta[i], yieldsSNII.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
	  MPI_Pack(yieldsSNII.TotalMetals[i], yieldsSNII.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
	  for(j = 0; j < yieldsSNII.N_ELEMENTS; j++)
	    MPI_Pack(yieldsSNII.Yield[i][j], yieldsSNII.N_MASS, MPI_DOUBLE, buffer, bufsize, &position,
		     MPI_COMM_WORLD);
	}
#ifdef BG_DUST
      MPI_Pack(dustSNII.Mass, dustSNII.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(dustSNII.Metallicity, dustSNII.N_Z, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
      for(i = 0; i < dustSNII.N_Z; i++)
        MPI_Pack(dustSNII.TotalDust[i], dustSNII.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
#endif

      /* pack AGB stars table */
      for(i = 0; i < yieldsAGB.N_ELEMENTS; i++)
	MPI_Pack(yieldsAGB.ElementName[i], EL_NAME_LENGTH, MPI_CHAR, buffer, bufsize, &position,
		 MPI_COMM_WORLD);
      MPI_Pack(yieldsAGB.Mass, yieldsAGB.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(yieldsAGB.Metallicity, yieldsAGB.N_Z, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
      for(i = 0; i < yieldsAGB.N_Z; i++)
	{
	  MPI_Pack(yieldsAGB.Ejecta[i], yieldsAGB.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
	  MPI_Pack(yieldsAGB.TotalMetals[i], yieldsAGB.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
	  for(j = 0; j < yieldsAGB.N_ELEMENTS; j++)
	    MPI_Pack(yieldsAGB.Yield[i][j], yieldsAGB.N_MASS, MPI_DOUBLE, buffer, bufsize, &position,
		     MPI_COMM_WORLD);
	}
#ifdef BG_DUST
      MPI_Pack(dustAGB.Mass, dustAGB.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(dustAGB.Metallicity, dustAGB.N_Z, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
      for(i = 0; i < dustAGB.N_Z; i++)
        MPI_Pack(dustAGB.TotalDust[i], dustAGB.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
#endif

      /* pack POPIII table */
#ifdef BG_POPIII
      for(i = 0; i < yieldsPOPIII.N_ELEMENTS; i++)
	MPI_Pack(yieldsPOPIII.ElementName[i], POPIII_YIELD_EL_NAME_LENGTH, MPI_CHAR, buffer, bufsize, &position,
		 MPI_COMM_WORLD);
      MPI_Pack(yieldsPOPIII.Mass, yieldsPOPIII.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(yieldsPOPIII.Ejecta, yieldsPOPIII.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(yieldsPOPIII.TotalMetals, yieldsPOPIII.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
      for(i = 0; i < yieldsPOPIII.N_ELEMENTS; i++)
	MPI_Pack(yieldsPOPIII.Yield[i], yieldsPOPIII.N_MASS, MPI_DOUBLE, buffer, bufsize, &position,
		 MPI_COMM_WORLD);
#endif
#ifdef BG_DUST
      MPI_Pack(dustPOPIII.Mass, dustPOPIII.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(dustPOPIII.TotalDust, dustPOPIII.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
#endif
      /* pack Lifetimes table */
      MPI_Pack(Lifetimes.Mass, Lifetimes.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(Lifetimes.Metallicity, Lifetimes.N_Z, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
      for(i = 0; i < Lifetimes.N_Z; i++)
	MPI_Pack(Lifetimes.Dyingtime[i], Lifetimes.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
#ifdef BG_POPIII
      MPI_Pack(POPIIILifetimes.Mass, POPIIILifetimes.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack(POPIIILifetimes.Dyingtime, POPIIILifetimes.N_MASS, MPI_DOUBLE, buffer, bufsize, &position, MPI_COMM_WORLD);
#endif
    }

  MPI_Bcast(buffer, bufsize, MPI_PACKED, 0, MPI_COMM_WORLD);

  if(ThisTask != 0)
    {
      position = 0;

      /* unpack SNIa table */
      for(i = 0; i < yieldsSNIa.N_ELEMENTS; i++)
	MPI_Unpack(buffer, bufsize, &position, yieldsSNIa.ElementName[i], EL_NAME_LENGTH, MPI_CHAR,
		   MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, yieldsSNIa.Yield, yieldsSNIa.N_ELEMENTS, MPI_DOUBLE,
		 MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, &yieldsSNIa.TotalMetals_SPH, 1, MPI_DOUBLE, MPI_COMM_WORLD);

      /* unpack SNII table */
      for(i = 0; i < yieldsSNII.N_ELEMENTS; i++)
	MPI_Unpack(buffer, bufsize, &position, yieldsSNII.ElementName[i], EL_NAME_LENGTH, MPI_CHAR,
		   MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, yieldsSNII.Mass, yieldsSNII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, yieldsSNII.Metallicity, yieldsSNII.N_Z, MPI_DOUBLE,
		 MPI_COMM_WORLD);
      for(i = 0; i < yieldsSNII.N_Z; i++)
	{
	  MPI_Unpack(buffer, bufsize, &position, yieldsSNII.Ejecta[i], yieldsSNII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
	  MPI_Unpack(buffer, bufsize, &position, yieldsSNII.TotalMetals[i], yieldsSNII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
	  for(j = 0; j < yieldsSNII.N_ELEMENTS; j++)
	    MPI_Unpack(buffer, bufsize, &position, yieldsSNII.Yield[i][j], yieldsSNII.N_MASS, MPI_DOUBLE,
		       MPI_COMM_WORLD);
	}
#ifdef BG_DUST
      MPI_Unpack(buffer, bufsize, &position, dustSNII.Mass, dustSNII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, dustSNII.Metallicity, dustSNII.N_Z, MPI_DOUBLE, MPI_COMM_WORLD);
      for(i = 0; i < dustSNII.N_Z; i++)
        MPI_Unpack(buffer, bufsize, &position, dustSNII.TotalDust[i], dustSNII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
#endif

      /* unpack AGB stars table */
      for(i = 0; i < yieldsAGB.N_ELEMENTS; i++)
	MPI_Unpack(buffer, bufsize, &position, yieldsAGB.ElementName[i], EL_NAME_LENGTH, MPI_CHAR,
		   MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, yieldsAGB.Mass, yieldsAGB.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, yieldsAGB.Metallicity, yieldsAGB.N_Z, MPI_DOUBLE,
		 MPI_COMM_WORLD);
      for(i = 0; i < yieldsAGB.N_Z; i++)
	{
	  MPI_Unpack(buffer, bufsize, &position, yieldsAGB.Ejecta[i], yieldsAGB.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
	  MPI_Unpack(buffer, bufsize, &position, yieldsAGB.TotalMetals[i], yieldsAGB.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
	  for(j = 0; j < yieldsAGB.N_ELEMENTS; j++)
	    MPI_Unpack(buffer, bufsize, &position, yieldsAGB.Yield[i][j], yieldsAGB.N_MASS, MPI_DOUBLE,
		       MPI_COMM_WORLD);
	}
#ifdef BG_DUST
      MPI_Unpack(buffer, bufsize, &position, dustAGB.Mass, dustAGB.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, dustAGB.Metallicity, dustAGB.N_Z, MPI_DOUBLE, MPI_COMM_WORLD);
      for(i = 0; i < dustAGB.N_Z; i++)
        MPI_Unpack(buffer, bufsize, &position, dustAGB.TotalDust[i], dustAGB.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
#endif

      /* unpack POPIII table */
#ifdef BG_POPIII
      for(i = 0; i < yieldsPOPIII.N_ELEMENTS; i++)
        MPI_Unpack(buffer, bufsize, &position, yieldsPOPIII.ElementName[i], POPIII_YIELD_EL_NAME_LENGTH, MPI_CHAR,
                   MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, yieldsPOPIII.Mass, yieldsPOPIII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, yieldsPOPIII.Ejecta, yieldsPOPIII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, yieldsPOPIII.TotalMetals, yieldsPOPIII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
      for(i = 0; i < yieldsPOPIII.N_ELEMENTS; i++)
	MPI_Unpack(buffer, bufsize, &position, yieldsPOPIII.Yield[i], yieldsPOPIII.N_MASS, MPI_DOUBLE,
		   MPI_COMM_WORLD);
#ifdef BG_DUST
      MPI_Unpack(buffer, bufsize, &position, dustPOPIII.Mass, dustPOPIII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, dustPOPIII.TotalDust, dustPOPIII.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
#endif
#endif
      /* unpack Lifetimes table */
      MPI_Unpack(buffer, bufsize, &position, Lifetimes.Mass, Lifetimes.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, Lifetimes.Metallicity, Lifetimes.N_Z, MPI_DOUBLE,
		 MPI_COMM_WORLD);
      for(i = 0; i < Lifetimes.N_Z; i++)
	MPI_Unpack(buffer, bufsize, &position, Lifetimes.Dyingtime[i], Lifetimes.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
#ifdef BG_POPIII
      MPI_Unpack(buffer, bufsize, &position, POPIIILifetimes.Mass, POPIIILifetimes.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufsize, &position, POPIIILifetimes.Dyingtime, POPIIILifetimes.N_MASS, MPI_DOUBLE, MPI_COMM_WORLD);
#endif
    }

  myfree(buffer);

  /*
  if(ThisTask == 1)
    for(j = 0; j < yieldsPOPIII.N_ELEMENTS; j++)
      printf("[POPIII] yield[%d][0] = %e\n", j, yieldsPOPIII.Yield[j][0]);
  fflush(stdout);

  if(ThisTask == 1)
    for(j = 0; j < yieldsPOPIII.N_ELEMENTS; j++)
      printf("[POPIII] yield[%d][%d] = %e\n", j, yieldsPOPIII.N_MASS-1, yieldsPOPIII.Yield[j][yieldsPOPIII.N_MASS-1]);
  fflush(stdout);
  */

  /* debug */
#ifdef BG_DUST
  if(ThisTask == 1)
    for(j = 0; j < dustPOPIII.N_MASS; j++)
      printf("[POPIII] dust[%d] = %e\n", j, dustPOPIII.TotalDust[j]);
  fflush(stdout);
#endif

  return 0;
}


/*
 * ----------------------------------------------------------------------
 * This routine allocates the yield table arrays.
 * ----------------------------------------------------------------------
 */

int allocate_yield_tables()
{
  int i, j;

  /* ------------------- */
  /* allocate SNIa table */
  /* ------------------- */

  /* allocate element name array */
  yieldsSNIa.ElementName = (char **) mymalloc(yieldsSNIa.N_ELEMENTS * sizeof(char *));

  if(ThisTask != 0)
    for(i = 0; i < yieldsSNIa.N_ELEMENTS; i++)
      yieldsSNIa.ElementName[i] = (char *) mymalloc(EL_NAME_LENGTH * sizeof(char));

  /* allocate yield array */
  yieldsSNIa.Yield = (double *) mymalloc(yieldsSNIa.N_ELEMENTS * sizeof(double));

  /* ------------------- */
  /* allocate SNII table */
  /* ------------------- */

  /* allocate element name array */
  yieldsSNII.ElementName = (char **) mymalloc(yieldsSNII.N_ELEMENTS * sizeof(char *));

  if(ThisTask != 0)
    for(i = 0; i < yieldsSNII.N_ELEMENTS; i++)
      yieldsSNII.ElementName[i] = (char *) mymalloc(EL_NAME_LENGTH * sizeof(char));

  /* allocate mass array */
  yieldsSNII.Mass = (double *) mymalloc(yieldsSNII.N_MASS * sizeof(double));

  /* allocate metallicity array */
  yieldsSNII.Metallicity = (double *) mymalloc(yieldsSNII.N_Z * sizeof(double));

  /* allocate yield array */
  yieldsSNII.Yield = (double ***) mymalloc(yieldsSNII.N_Z * sizeof(double **));
  yieldsSNII.Ejecta = (double **) mymalloc(yieldsSNII.N_Z * sizeof(double *));
  yieldsSNII.TotalMetals = (double **) mymalloc(yieldsSNII.N_Z * sizeof(double *));

  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      yieldsSNII.Yield[i] = (double **) mymalloc(yieldsSNII.N_ELEMENTS * sizeof(double *));
      yieldsSNII.Ejecta[i] = (double *) mymalloc(yieldsSNII.N_MASS * sizeof(double));
      yieldsSNII.TotalMetals[i] = (double *) mymalloc(yieldsSNII.N_MASS * sizeof(double));

      for(j = 0; j < yieldsSNII.N_ELEMENTS; j++)
	yieldsSNII.Yield[i][j] = (double *) mymalloc(yieldsSNII.N_MASS * sizeof(double));
    }

#ifdef BG_DUST
  /* allocate mass array */
  dustSNII.Mass = (double *) mymalloc(dustSNII.N_MASS * sizeof(double));

  /* allocate metallicity array */
  dustSNII.Metallicity = (double *) mymalloc(dustSNII.N_Z * sizeof(double));

  /* allocate yield array */
  dustSNII.TotalDust = (double **) mymalloc(dustSNII.N_Z * sizeof(double *));

  for(i = 0; i < dustSNII.N_Z; i++)
    dustSNII.TotalDust[i] = (double *) mymalloc(dustSNII.N_MASS * sizeof(double));
#endif

  /* ------------------ */
  /* allocate AGB table */
  /* ------------------ */

  /* allocate element name array */
  yieldsAGB.ElementName = (char **) mymalloc(yieldsAGB.N_ELEMENTS * sizeof(char *));

  if(ThisTask != 0)
    for(i = 0; i < yieldsAGB.N_ELEMENTS; i++)
      yieldsAGB.ElementName[i] = (char *) mymalloc(EL_NAME_LENGTH * sizeof(char));

  /* allocate mass array */
  yieldsAGB.Mass = (double *) mymalloc(yieldsAGB.N_MASS * sizeof(double));

  /* allocate metallicity array */
  yieldsAGB.Metallicity = (double *) mymalloc(yieldsAGB.N_Z * sizeof(double));

  /* allocate yield array */
  yieldsAGB.Yield = (double ***) mymalloc(yieldsAGB.N_Z * sizeof(double **));
  yieldsAGB.Ejecta = (double **) mymalloc(yieldsAGB.N_Z * sizeof(double *));
  yieldsAGB.TotalMetals = (double **) mymalloc(yieldsAGB.N_Z * sizeof(double *));

  for(i = 0; i < yieldsAGB.N_Z; i++)
    {
      yieldsAGB.Yield[i] = (double **) mymalloc(yieldsAGB.N_ELEMENTS * sizeof(double *));
      yieldsAGB.Ejecta[i] = (double *) mymalloc(yieldsAGB.N_MASS * sizeof(double));
      yieldsAGB.TotalMetals[i] = (double *) mymalloc(yieldsAGB.N_MASS * sizeof(double));

      for(j = 0; j < yieldsAGB.N_ELEMENTS; j++)
	yieldsAGB.Yield[i][j] = (double *) mymalloc(yieldsAGB.N_MASS * sizeof(double));
    }

#ifdef BG_DUST
  /* allocate mass array */
  dustAGB.Mass = (double *) mymalloc(dustAGB.N_MASS * sizeof(double));

  /* allocate metallicity array */
  dustAGB.Metallicity = (double *) mymalloc(dustAGB.N_Z * sizeof(double));

  /* allocate yield array */
  dustAGB.TotalDust = (double **) mymalloc(dustAGB.N_Z * sizeof(double *));

  for(i = 0; i < dustAGB.N_Z; i++)
    dustAGB.TotalDust[i] = (double *) mymalloc(dustAGB.N_MASS * sizeof(double));
#endif

  /* --------------------- */
  /* allocate POPIII table */
  /* --------------------- */
#ifdef BG_POPIII
  /* allocate element name array */
  yieldsPOPIII.ElementName = (char **) mymalloc(yieldsPOPIII.N_ELEMENTS * sizeof(char *));

  if(ThisTask != 0)
    for(i = 0; i < yieldsPOPIII.N_ELEMENTS; i++)
      yieldsPOPIII.ElementName[i] = (char *) mymalloc(POPIII_YIELD_EL_NAME_LENGTH * sizeof(char));

  /* allocate mass array */
  yieldsPOPIII.Mass = (double *) mymalloc(yieldsPOPIII.N_MASS * sizeof(double));

  /* allocate yield array */
  yieldsPOPIII.Yield = (double **) mymalloc(yieldsPOPIII.N_ELEMENTS * sizeof(double *));
  yieldsPOPIII.Ejecta = (double *) mymalloc(yieldsPOPIII.N_MASS * sizeof(double));
  yieldsPOPIII.TotalMetals = (double *) mymalloc(yieldsPOPIII.N_MASS * sizeof(double));

  for(i = 0; i < yieldsPOPIII.N_ELEMENTS; i++)
    yieldsPOPIII.Yield[i] = (double *) mymalloc(yieldsPOPIII.N_MASS * sizeof(double));
#ifdef BG_DUST
  /* allocate mass array */
  dustPOPIII.Mass = (double *) mymalloc(dustPOPIII.N_MASS * sizeof(double));

  /* allocate yield array */
  dustPOPIII.TotalDust = (double *) mymalloc(dustPOPIII.N_MASS * sizeof(double));
#endif
#endif

  /* ----------------------- */
  /* allocate Lifetime table */
  /* ----------------------- */

  /* allocate mass array */
  Lifetimes.Mass = (double *) mymalloc(Lifetimes.N_MASS * sizeof(double));

  /* allocate metallicity array */
  Lifetimes.Metallicity = (double *) mymalloc(Lifetimes.N_Z * sizeof(double));

  /* allocate lifetime array */
  Lifetimes.Dyingtime = (double **) mymalloc(Lifetimes.N_Z * sizeof(double *));

  for(i = 0; i < Lifetimes.N_Z; i++)
    Lifetimes.Dyingtime[i] = (double *) mymalloc(Lifetimes.N_MASS * sizeof(double));

#ifdef BG_POPIII
  /* allocate mass array */
  POPIIILifetimes.Mass = (double *) mymalloc(POPIIILifetimes.N_MASS * sizeof(double));

  /* allocate lifetime array */
  POPIIILifetimes.Dyingtime = (double *) mymalloc(POPIIILifetimes.N_MASS * sizeof(double));
#endif

  return 0;
}

#endif
