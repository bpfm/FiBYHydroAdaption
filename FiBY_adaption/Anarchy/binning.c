#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <hdf5.h>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "allvars.h"
#include "proto.h"
#include "bg_proto.h"
#include "bg_vars.h"

#define N_Zbins 13
#define N_Tbins 14
#define N_Dbins 20

#if defined(BG_SFR) && defined(BG_STELLAR_EVOLUTION) && defined(BG_OUTPUT_GRID)

void binning(void)
{
  hid_t hdf5_file = 0, dataset = 0;	/* file and dataset handles  */

  hid_t datatype = 0, dataspace = 0, hdf5_headergrp = 0;	/* dataspace handle  */

  hid_t hdf5_paramgrp = 0, hdf5_unitsgrp = 0, hdf5_constgrp = 0;

  hid_t hdf5_gridgrp = 0, hdf5_elementgrp = 0;

  /* dataset dimensions  */
  hsize_t dim1d[1] = { 1 };
  hsize_t dim1d_temp[1] = { N_Tbins };
  hsize_t dim1d_dens[1] = { N_Dbins };
  hsize_t dim1d_metallicity[1] = { N_Zbins };
  hsize_t dim2d_temp[2] = { N_Zbins, N_Dbins };
  hsize_t dim2d_dens[2] = { N_Zbins, N_Tbins };
  hsize_t dim2d_metallicity[2] = { N_Tbins, N_Dbins };
  hsize_t dim3d[3] = { N_Zbins, N_Tbins, N_Dbins };

  herr_t status;

  char hdfname[100], buf[500];

  char comment[] =
    "This grid output represents the total fraction of the baryon mass distributed as per its metalicity, temperature, and density (indexed in that order).  The bins should be taken as seperators, thus the extreme bins are open ended. Enjoy!";

  int i, j, k, l, n, index_Z, index_T, index_D;

  double time;

  double dTm1, dDm1, dZm1, RhoCrit, rho;

  double baryon_mass = 0, global_baryon_mass;

  double metal_mass = 0, global_metal_mass;

  double llZ = -5.0, ulZ = 0.5, llT = 2.0, ulT = 8.0, llD = -2.0, ulD = 7.0;

  double this_element[N_Zbins][N_Tbins][N_Dbins];

  double mass_bins_sf[N_Zbins][N_Tbins][N_Dbins],
    mass_bins_nsf[N_Zbins][N_Tbins][N_Dbins],
    metal_mass_bins_sf[N_Tbins][N_Dbins],
    metal_mass_bins_nsf[N_Tbins][N_Dbins],
    el_mass_bins_sf[N_Zbins][N_Tbins][N_Dbins][BG_NELEMENTS],
    el_mass_bins_nsf[N_Zbins][N_Tbins][N_Dbins][BG_NELEMENTS],
#ifdef BG_SNIA_IRON
    snia_iron_bins_sf[N_Zbins][N_Tbins][N_Dbins], snia_iron_bins_nsf[N_Zbins][N_Tbins][N_Dbins],
#endif
#ifdef BG_Z_WEIGHTED_REDSHIFT
    z_weighted_redshift_sf[N_Tbins][N_Dbins], z_weighted_redshift_nsf[N_Tbins][N_Dbins],
#endif
   
    mean_temp_sf[N_Zbins][N_Dbins],
    mean_dens_sf[N_Zbins][N_Tbins],
    mean_metallicity_sf[N_Tbins][N_Dbins],
    mean_temp_nsf[N_Zbins][N_Dbins],
    mean_dens_nsf[N_Zbins][N_Tbins],
    mean_metallicity_nsf[N_Tbins][N_Dbins], Zdivid[N_Zbins - 1], Tdivid[N_Tbins - 1], Ddivid[N_Dbins - 1];

  double mean_temp_dens_sf[N_Dbins],
    mean_dens_temp_sf[N_Tbins],
    mean_temp_met_sf[N_Zbins],
    mean_dens_met_sf[N_Zbins],
    global_mean_temp_dens_sf[N_Dbins],
    global_mean_dens_temp_sf[N_Tbins], global_mean_temp_met_sf[N_Zbins], global_mean_dens_met_sf[N_Zbins];

  double mean_temp_dens_nsf[N_Dbins],
    mean_dens_temp_nsf[N_Tbins],
    mean_temp_met_nsf[N_Zbins],
    mean_dens_met_nsf[N_Zbins],
    global_mean_temp_dens_nsf[N_Dbins],
    global_mean_dens_temp_nsf[N_Tbins], global_mean_temp_met_nsf[N_Zbins], global_mean_dens_met_nsf[N_Zbins];

  double global_mass_bins_sf[N_Zbins][N_Tbins][N_Dbins],
    global_mass_bins_nsf[N_Zbins][N_Tbins][N_Dbins],
    global_metal_mass_bins_sf[N_Tbins][N_Dbins],
    global_metal_mass_bins_nsf[N_Tbins][N_Dbins],
    global_el_mass_bins_sf[N_Zbins][N_Tbins][N_Dbins][BG_NELEMENTS],
    global_el_mass_bins_nsf[N_Zbins][N_Tbins][N_Dbins][BG_NELEMENTS];
#ifdef BG_SNIA_IRON
  double global_snia_iron_bins_sf[N_Zbins][N_Tbins][N_Dbins],
    global_snia_iron_bins_nsf[N_Zbins][N_Tbins][N_Dbins];
#endif
  double global_mean_dens_sf[N_Zbins][N_Tbins], global_mean_dens_nsf[N_Zbins][N_Tbins];

  double global_mean_temp_sf[N_Zbins][N_Dbins], global_mean_temp_nsf[N_Zbins][N_Dbins];

  double global_mean_metallicity_sf[N_Tbins][N_Dbins], global_mean_metallicity_nsf[N_Tbins][N_Dbins];

  double global_mean_dens_mass_sf[N_Zbins][N_Tbins], global_mean_dens_mass_nsf[N_Zbins][N_Tbins];

  double global_mean_temp_mass_sf[N_Zbins][N_Dbins], global_mean_temp_mass_nsf[N_Zbins][N_Dbins];

  double global_mean_metallicity_mass_sf[N_Tbins][N_Dbins],
    global_mean_metallicity_mass_nsf[N_Tbins][N_Dbins];
#ifdef BG_Z_WEIGHTED_REDSHIFT
  double global_z_weighted_redshift_sf[N_Tbins][N_Dbins], global_z_weighted_redshift_nsf[N_Tbins][N_Dbins];
#endif

  double log_Z, log_T, log_D, Z;

  double a3inv = 1, redshift = 0;

  double temp;

  if(ThisTask == 0)
    {
      sprintf(buf, "%s/grid", All.OutputDir);
      mkdir(buf, 02755);
    }

  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      redshift = 1 / All.Time - 1;
    }
  else
    {
      a3inv = 1;
      redshift = 0;
    }

  /* update header informations */
/*   header.time = All.Time; */
/*   header.redshift = redshift; */



  /* determine particle numbers of each type in file */

/*   if(ThisTask == writeTask) */
/*     { */
/*       for(n = 0; n < 6; n++) */
/* 	ntot_type[n] = n_type[n]; */

/*       for(task = writeTask + 1; task <= lastTask; task++) */
/* 	{ */
/* 	  MPI_Recv(&nn[0], 6, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status); */
/* 	  for(n = 0; n < 6; n++) */
/* 	    ntot_type[n] += nn[n]; */
/* 	} */

/*       for(task = writeTask + 1; task <= lastTask; task++) */
/* 	MPI_Send(&ntot_type[0], 6, MPI_INT, task, TAG_N, MPI_COMM_WORLD); */
/*     } */
/*   else */
/*     { */
/*       MPI_Send(&n_type[0], 6, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD); */
/*       MPI_Recv(&ntot_type[0], 6, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status); */
/*     } */



  /* fill file header */

/*   for(n = 0; n < 6; n++) */
/*     { */
/*       header.npart[n] = ntot_type[n]; */
/*       header.npartTotal[n] = (unsigned int) ntot_type_all[n]; */
/*       header.npartTotalHighWord[n] = (unsigned int) (ntot_type_all[n] >> 32); */
/*     } */

  for(n = 0; n < 6; n++)
    header.mass[n] = All.MassTable[n];

  header.time = All.Time;

  if(All.ComovingIntegrationOn)
    header.redshift = redshift;
  else
    header.redshift = 0;

  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;
  header.num_files = All.NumFilesPerSnapshot;
  header.BoxSize = All.BoxSize;
  header.Omega0 = All.Omega0;
  header.OmegaLambda = All.OmegaLambda;
  header.HubbleParam = All.HubbleParam;


  RhoCrit = All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);


  /* Initialise the grid  */
  for(i = 0; i < N_Zbins; i++)
    for(j = 0; j < N_Tbins; j++)
      for(k = 0; k < N_Dbins; k++)
	{
	  mass_bins_sf[i][j][k] = 0;
	  mass_bins_nsf[i][j][k] = 0;
#ifdef BG_SNIA_IRON
	  snia_iron_bins_sf[i][j][k] = 0;
	  snia_iron_bins_nsf[i][j][k] = 0;
#endif
	  for(l = 0; l < BG_NELEMENTS; l++)
	    {
	      el_mass_bins_sf[i][j][k][l] = 0;
	      el_mass_bins_nsf[i][j][k][l] = 0;
	    }
	}

  for(i = 0; i < N_Zbins; i++)
    for(j = 0; j < N_Tbins; j++)
      {
	mean_dens_sf[i][j] = 0;
	mean_dens_nsf[i][j] = 0;
      }

  for(i = 0; i < N_Zbins; i++)
    for(k = 0; k < N_Dbins; k++)
      {
	mean_temp_sf[i][k] = 0;
	mean_temp_nsf[i][k] = 0;
      }

  for(j = 0; j < N_Tbins; j++)
    for(k = 0; k < N_Dbins; k++)
      {
	mean_metallicity_sf[j][k] = 0;
	mean_metallicity_nsf[j][k] = 0;
	metal_mass_bins_sf[j][k] = 0;
	metal_mass_bins_nsf[j][k] = 0;
#ifdef BG_Z_WEIGHTED_REDSHIFT
	z_weighted_redshift_sf[j][k] = 0;
	z_weighted_redshift_nsf[j][k] = 0;
#endif
      }

  for(k = 0; k < N_Dbins; k++)
    {
      mean_temp_dens_sf[k] = 0;
      mean_temp_dens_nsf[k] = 0;
    }

  for(j = 0; j < N_Tbins; j++)
    {
      mean_dens_temp_sf[j] = 0;
      mean_dens_temp_nsf[j] = 0;
    }

  for(i = 0; i < N_Zbins; i++)
    {
      mean_temp_met_sf[i] = 0;
      mean_dens_met_sf[i] = 0;
      mean_temp_met_nsf[i] = 0;
      mean_dens_met_nsf[i] = 0;
    }


  /* Initialise the binning intervals  */
  dTm1 = 1 / (ulT - llT);
  dDm1 = 1 / (ulD - llD);
  dZm1 = 1 / (ulZ - llZ);

  /* compute global baryon (gas + stars) and metal (gas) mass */
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0 || P[i].Type == 4)
	baryon_mass += P[i].Mass;

      if(P[i].Type == 0)
	metal_mass += P[i].Mass * SphP[i].Metallicity;
    }

  MPI_Allreduce(&baryon_mass, &global_baryon_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&metal_mass, &global_metal_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* Loop over every particle, converting the entropy to temperature, and binning them  */
  for(i = 0; i < N_gas; i++)
    {
      drift_particle(i, All.Ti_Current);
#ifdef BG_METALSMOOTHING
      Z = SphP[i].MetallicitySmoothed / ZSOLAR;
#else
      Z = SphP[i].Metallicity / ZSOLAR;
#endif
      if(Z > 0)
	log_Z = log10(Z);
      else
	log_Z = llZ;

      rho = SphP[i].d.Density;
      log_D = log10(rho / RhoCrit);
      log_T = log10(bg_get_temperature(i));

      /* Metallicity  */
      if(log_Z < llZ)
	index_Z = 0;
      else if(log_Z > ulZ)
	index_Z = N_Zbins - 1;
      else
	index_Z = (int) ceil((N_Zbins - 2) * (log_Z - llZ) * dZm1);

      /* Temperature  */
      if(log_T < llT)
	index_T = 0;
      else if(log_T > ulT)
	index_T = N_Tbins - 1;
      else
	index_T = (int) ceil((N_Tbins - 2) * (log_T - llT) * dTm1);

      /* Density  */
      if(log_D < llD)
	index_D = 0;
      else if(log_D > ulD)
	index_D = N_Dbins - 1;
      else
	index_D = (int) ceil((N_Dbins - 2) * (log_D - llD) * dDm1);

      /* Put particle in the grid  */
      if(SphP[i].Sfr == 0)
	{
	  mass_bins_nsf[index_Z][index_T][index_D] += P[i].Mass;
#ifdef BG_SNIA_IRON
	  snia_iron_bins_nsf[index_Z][index_T][index_D] += SphP[i].IronFromSNIa;
#endif
	  /* mass weighted temperature */
	  mean_temp_nsf[index_Z][index_D] += P[i].Mass * bg_get_temperature(i);

	  /* mass weighted density */
	  mean_dens_nsf[index_Z][index_T] += P[i].Mass * SphP[i].d.Density;

	  /* mass weighted metallicity */
	  mean_metallicity_nsf[index_T][index_D] += P[i].Mass * SphP[i].Metallicity;

	  /* metal mass */
	  metal_mass_bins_nsf[index_T][index_D] += P[i].Mass * SphP[i].Metallicity;

	  /* elements mass */
	  for(l = 0; l < BG_NELEMENTS; l++)
	    el_mass_bins_nsf[index_Z][index_T][index_D][l] += SphP[i].Metals[l];

#ifdef BG_Z_WEIGHTED_REDSHIFT
	  z_weighted_redshift_nsf[index_T][index_D] +=
	    P[i].Mass * SphP[i].Metallicity * SphP[i].MetallicityWeightedRedshift;
#endif

	  mean_temp_dens_nsf[index_D] += P[i].Mass * bg_get_temperature(i);
	  mean_dens_temp_nsf[index_T] += P[i].Mass * SphP[i].d.Density;
	  mean_temp_met_nsf[index_Z] += P[i].Mass * bg_get_temperature(i);
	  mean_dens_met_nsf[index_Z] += P[i].Mass * SphP[i].d.Density;
	}
      else
	{
	  mass_bins_sf[index_Z][index_T][index_D] += P[i].Mass;
#ifdef BG_SNIA_IRON
	  snia_iron_bins_sf[index_Z][index_T][index_D] += SphP[i].IronFromSNIa;
#endif
	  /* mass weighted temperature */
	  mean_temp_sf[index_Z][index_D] += P[i].Mass * bg_get_temperature(i);

	  /* mass weighted density */
	  mean_dens_sf[index_Z][index_T] += P[i].Mass * SphP[i].d.Density;

	  /* mass weighted metallicity */
	  mean_metallicity_sf[index_T][index_D] += P[i].Mass * SphP[i].Metallicity;

	  /* metal mass */
	  metal_mass_bins_sf[index_T][index_D] += P[i].Mass * SphP[i].Metallicity;

	  /* elements mass */
	  for(l = 0; l < BG_NELEMENTS; l++)
	    el_mass_bins_sf[index_Z][index_T][index_D][l] += SphP[i].Metals[l];

#ifdef BG_Z_WEIGHTED_REDSHIFT
	  z_weighted_redshift_sf[index_T][index_D] +=
	    P[i].Mass * SphP[i].Metallicity * SphP[i].MetallicityWeightedRedshift;
#endif

	  mean_temp_dens_sf[index_D] += P[i].Mass * bg_get_temperature(i);
	  mean_dens_temp_sf[index_T] += P[i].Mass * SphP[i].d.Density;
	  mean_temp_met_sf[index_Z] += P[i].Mass * bg_get_temperature(i);
	  mean_dens_met_sf[index_Z] += P[i].Mass * SphP[i].d.Density;
	}
    }

  /* gas mass */
  MPI_Reduce(mass_bins_nsf, global_mass_bins_nsf, N_Zbins * N_Tbins * N_Dbins,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(mass_bins_sf, global_mass_bins_sf, N_Zbins * N_Tbins * N_Dbins,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* metal mass */
  MPI_Reduce(metal_mass_bins_nsf, global_metal_mass_bins_nsf, N_Tbins * N_Dbins,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(metal_mass_bins_sf, global_metal_mass_bins_sf, N_Tbins * N_Dbins,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* elements mass */
  MPI_Reduce(el_mass_bins_nsf, global_el_mass_bins_nsf, N_Zbins * N_Tbins * N_Dbins * BG_NELEMENTS,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(el_mass_bins_sf, global_el_mass_bins_sf, N_Zbins * N_Tbins * N_Dbins * BG_NELEMENTS,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

#ifdef BG_Z_WEIGHTED_REDSHIFT
  MPI_Reduce(z_weighted_redshift_sf, global_z_weighted_redshift_sf, N_Tbins * N_Dbins,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(z_weighted_redshift_nsf, global_z_weighted_redshift_nsf, N_Tbins * N_Dbins,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

  /* iron from SNIa */
#ifdef BG_SNIA_IRON
  MPI_Reduce(snia_iron_bins_nsf, global_snia_iron_bins_nsf, N_Zbins * N_Tbins * N_Dbins,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(snia_iron_bins_sf, global_snia_iron_bins_sf, N_Zbins * N_Tbins * N_Dbins,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

  /* mean temperatures */
  MPI_Reduce(mean_temp_nsf, global_mean_temp_nsf, N_Zbins * N_Dbins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(mean_temp_sf, global_mean_temp_sf, N_Zbins * N_Dbins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(mean_temp_dens_nsf, global_mean_temp_dens_nsf, N_Dbins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(mean_temp_dens_sf, global_mean_temp_dens_sf, N_Dbins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(mean_temp_met_nsf, global_mean_temp_met_nsf, N_Zbins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(mean_temp_met_sf, global_mean_temp_met_sf, N_Zbins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* mean densities */
  MPI_Reduce(mean_dens_nsf, global_mean_dens_nsf, N_Zbins * N_Tbins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(mean_dens_sf, global_mean_dens_sf, N_Zbins * N_Tbins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(mean_dens_temp_nsf, global_mean_dens_temp_nsf, N_Tbins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(mean_dens_temp_sf, global_mean_dens_temp_sf, N_Tbins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(mean_dens_met_nsf, global_mean_dens_met_nsf, N_Zbins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(mean_dens_met_sf, global_mean_dens_met_sf, N_Zbins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* mean metallicity */
  MPI_Reduce(mean_metallicity_nsf, global_mean_metallicity_nsf, N_Tbins * N_Dbins,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(mean_metallicity_sf, global_mean_metallicity_sf, N_Tbins * N_Dbins,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


  /* Output */
  if(ThisTask == 0)
    {
      for(i = 0; i < N_Zbins; i++)
	for(k = 0; k < N_Dbins; k++)
	  {
	    global_mean_temp_mass_nsf[i][k] = 0;
	    global_mean_temp_mass_sf[i][k] = 0;
	  }

      for(i = 0; i < N_Zbins; i++)
	for(j = 0; j < N_Tbins; j++)
	  {
	    global_mean_dens_mass_nsf[i][j] = 0;
	    global_mean_dens_mass_sf[i][j] = 0;
	  }

      for(j = 0; j < N_Tbins; j++)
	for(k = 0; k < N_Dbins; k++)
	  {
	    global_mean_metallicity_mass_nsf[j][k] = 0;
	    global_mean_metallicity_mass_sf[j][k] = 0;
	  }

      for(i = 0; i < N_Zbins; i++)
	for(j = 0; j < N_Tbins; j++)
	  for(k = 0; k < N_Dbins; k++)
	    {
	      global_mean_temp_mass_nsf[i][k] += global_mass_bins_nsf[i][j][k];
	      global_mean_dens_mass_nsf[i][j] += global_mass_bins_nsf[i][j][k];
	      global_mean_metallicity_mass_nsf[j][k] += global_mass_bins_nsf[i][j][k];

	      global_mean_temp_mass_sf[i][k] += global_mass_bins_sf[i][j][k];
	      global_mean_dens_mass_sf[i][j] += global_mass_bins_sf[i][j][k];
	      global_mean_metallicity_mass_sf[j][k] += global_mass_bins_sf[i][j][k];
	    }

      double mass_Z_nsf[N_Zbins], mass_T_nsf[N_Tbins], mass_D_nsf[N_Dbins];

      double mass_Z_sf[N_Zbins], mass_T_sf[N_Tbins], mass_D_sf[N_Dbins];

      for(i = 0; i < N_Zbins; i++)
	{
	  mass_Z_nsf[i] = 0;
	  mass_Z_sf[i] = 0;
	}

      for(j = 0; j < N_Tbins; j++)
	{
	  mass_T_nsf[j] = 0;
	  mass_T_sf[j] = 0;
	}

      for(k = 0; k < N_Dbins; k++)
	{
	  mass_D_nsf[k] = 0;
	  mass_D_sf[k] = 0;
	}

      for(i = 0; i < N_Zbins; i++)
	for(j = 0; j < N_Tbins; j++)
	  for(k = 0; k < N_Dbins; k++)
	    {
	      mass_Z_nsf[i] += global_mass_bins_nsf[i][j][k];
	      mass_Z_sf[i] += global_mass_bins_sf[i][j][k];
	    }

      for(i = 0; i < N_Zbins; i++)
	for(j = 0; j < N_Tbins; j++)
	  for(k = 0; k < N_Dbins; k++)
	    {
	      mass_T_nsf[j] += global_mass_bins_nsf[i][j][k];
	      mass_T_sf[j] += global_mass_bins_sf[i][j][k];
	    }

      for(i = 0; i < N_Zbins; i++)
	for(j = 0; j < N_Tbins; j++)
	  for(k = 0; k < N_Dbins; k++)
	    {
	      mass_D_nsf[k] += global_mass_bins_nsf[i][j][k];
	      mass_D_sf[k] += global_mass_bins_sf[i][j][k];
	    }

      /* convert element masses to mass fractions */
      for(l = 0; l < BG_NELEMENTS; l++)
	for(i = 0; i < N_Zbins; i++)
	  for(j = 0; j < N_Tbins; j++)
	    for(k = 0; k < N_Dbins; k++)
	      {
		if(global_mass_bins_nsf[i][j][k] > 0)
		  global_el_mass_bins_nsf[i][j][k][l] /= global_mass_bins_nsf[i][j][k];
		if(global_mass_bins_sf[i][j][k] > 0)
		  global_el_mass_bins_sf[i][j][k][l] /= global_mass_bins_sf[i][j][k];
	      }

      /* convert iron from SNIa mass to mass fraction */
      for(i = 0; i < N_Zbins; i++)
	for(j = 0; j < N_Tbins; j++)
	  for(k = 0; k < N_Dbins; k++)
	    {
#ifdef BG_SNIA_IRON
	      if(global_mass_bins_nsf[i][j][k] > 0)
		global_snia_iron_bins_nsf[i][j][k] /= global_mass_bins_nsf[i][j][k];
	      if(global_mass_bins_sf[i][j][k] > 0)
		global_snia_iron_bins_sf[i][j][k] /= global_mass_bins_sf[i][j][k];
#endif
	    }


      /* convert to mass fractions */
      /* this must be done after normalizing element masses */
      for(i = 0; i < N_Zbins; i++)
	for(j = 0; j < N_Tbins; j++)
	  for(k = 0; k < N_Dbins; k++)
	    {
	      global_mass_bins_nsf[i][j][k] /= global_baryon_mass;
	      global_mass_bins_sf[i][j][k] /= global_baryon_mass;
#ifdef BG_SNIA_IRON
	      global_snia_iron_bins_nsf[i][j][k] /= global_baryon_mass;
	      global_snia_iron_bins_sf[i][j][k] /= global_baryon_mass;
#endif
	    }


      /* normalise 1D arrays */
      for(i = 0; i < N_Zbins; i++)
	{
	  if(mass_Z_sf[i] > 0)
	    {
	      global_mean_temp_met_sf[i] /= mass_Z_sf[i];
	      global_mean_dens_met_sf[i] /= mass_Z_sf[i];
	    }
	  if(mass_Z_nsf[i] > 0)
	    {
	      global_mean_temp_met_nsf[i] /= mass_Z_nsf[i];
	      global_mean_dens_met_nsf[i] /= mass_Z_nsf[i];
	    }
	}

      for(j = 0; j < N_Tbins; j++)
	{
	  if(mass_T_sf[j] > 0)
	    global_mean_dens_temp_sf[j] /= mass_T_sf[j];
	  if(mass_T_nsf[j] > 0)
	    global_mean_dens_temp_nsf[j] /= mass_T_nsf[j];
	}

      for(k = 0; k < N_Dbins; k++)
	{
	  if(mass_D_sf[k] > 0)
	    global_mean_temp_dens_sf[k] /= mass_D_sf[k];
	  if(mass_D_nsf[k] > 0)
	    global_mean_temp_dens_nsf[k] /= mass_D_nsf[k];
	}


      /* normalise mean density */
      for(i = 0; i < N_Zbins; i++)
	for(j = 0; j < N_Tbins; j++)
	  {
	    if(global_mean_dens_mass_sf[i][j] > 0)
	      global_mean_dens_sf[i][j] /= global_mean_dens_mass_sf[i][j];
	    if(global_mean_dens_mass_nsf[i][j] > 0)
	      global_mean_dens_nsf[i][j] /= global_mean_dens_mass_nsf[i][j];
	  }

      /* normalise mean temperature */
      for(i = 0; i < N_Zbins; i++)
	for(k = 0; k < N_Dbins; k++)
	  {
	    if(global_mean_temp_mass_sf[i][k] > 0)
	      global_mean_temp_sf[i][k] /= global_mean_temp_mass_sf[i][k];
	    if(global_mean_temp_mass_nsf[i][k] > 0)
	      global_mean_temp_nsf[i][k] /= global_mean_temp_mass_nsf[i][k];
	  }

      /* normalise mean metallicity */
      for(j = 0; j < N_Tbins; j++)
	for(k = 0; k < N_Dbins; k++)
	  {
	    if(global_mean_metallicity_mass_sf[j][k] > 0)
	      global_mean_metallicity_sf[j][k] /= global_mean_metallicity_mass_sf[j][k];
	    if(global_mean_metallicity_mass_nsf[j][k] > 0)
	      global_mean_metallicity_nsf[j][k] /= global_mean_metallicity_mass_nsf[j][k];
#ifdef BG_Z_WEIGHTED_REDSHIFT
	    if(global_metal_mass_bins_sf[j][k] > 0)
	      global_z_weighted_redshift_sf[j][k] /= global_metal_mass_bins_sf[j][k];
	    if(global_metal_mass_bins_nsf[j][k] > 0)
	      global_z_weighted_redshift_nsf[j][k] /= global_metal_mass_bins_nsf[j][k];
#endif
	  }


      /* convert metal mass to metal mass fraction */
      /* this must be done after normalizing the metal mass weighted quantities */
      for(j = 0; j < N_Tbins; j++)
	for(k = 0; k < N_Dbins; k++)
	  {
	    if(global_mean_metallicity_mass_sf[j][k] > 0)
	      global_metal_mass_bins_sf[j][k] /= global_mean_metallicity_mass_sf[j][k];
	    if(global_mean_metallicity_mass_nsf[j][k] > 0)
	      global_metal_mass_bins_nsf[j][k] /= global_mean_metallicity_mass_nsf[j][k];
	  }


      if(All.ComovingIntegrationOn)
	{
	  if(redshift >= 10.0)
	    sprintf(hdfname, "%s/grid/massdist_z%1.3f.hdf5", All.OutputDir, redshift);
	  else
	    sprintf(hdfname, "%s/grid/massdist_z0%1.3f.hdf5", All.OutputDir, redshift);
	}
      else
	sprintf(hdfname, "%s/grid/massdist_t%1.3f.hdf5", All.OutputDir, All.Time);

      hdf5_file = H5Fcreate(hdfname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      time = bg_get_elapsed_time(0, All.Time, 1);	/* note: this is valid for a flat universe! */


      /* write usual header info */
      hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);
      hdf5_paramgrp = H5Gcreate(hdf5_file, "/Parameters", 0);
      hdf5_unitsgrp = H5Gcreate(hdf5_file, "/Units", 0);
      hdf5_constgrp = H5Gcreate(hdf5_file, "/Constants", 0);

      write_header_attributes_in_hdf5(hdf5_headergrp);
      write_parameters_attributes_in_hdf5(hdf5_paramgrp);
      write_units_attributes_in_hdf5(hdf5_unitsgrp);
      write_constants_attributes_in_hdf5(hdf5_constgrp);

      H5Gclose(hdf5_headergrp);
      H5Gclose(hdf5_paramgrp);
      H5Gclose(hdf5_unitsgrp);
      H5Gclose(hdf5_constgrp);


      /* write grid data */
      hdf5_gridgrp = H5Gcreate(hdf5_file, "/Grid", 0);
      /* */
      temp = ZSOLAR;
      dataspace = H5Screate_simple(0, dim1d, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "ReferenceMetallicity_SOLAR",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &temp);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      /* Comment */
      dataspace = H5Screate(H5S_SCALAR);
      datatype = H5Tcopy(H5T_C_S1);
      status = H5Tset_size(datatype, sizeof(comment));
      dataset = H5Dcreate(hdf5_gridgrp, "Comment", datatype, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, comment);
      H5Sclose(dataspace);
      H5Dclose(dataset);
      H5Tclose(datatype);

      /* Bins */

      /* Metal */
      for(i = 0; i < N_Zbins - 1; i++)
	Zdivid[i] = llZ + i * (ulZ - llZ) / (N_Zbins - 2);

      dim1d[0] = N_Zbins - 1;
      dataspace = H5Screate_simple(1, dim1d, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "LogMetallicityBins_SOLAR",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Zdivid);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      /* Temperature */
      for(i = 0; i < N_Tbins - 1; i++)
	Tdivid[i] = llT + i * (ulT - llT) / (N_Tbins - 2);

      dim1d[0] = N_Tbins - 1;
      dataspace = H5Screate_simple(1, dim1d, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "LogTemperatureBins_K", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Tdivid);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      /* Density */
      for(i = 0; i < N_Dbins - 1; i++)
	Ddivid[i] = llD + i * (ulD - llD) / (N_Dbins - 2);

      dim1d[0] = N_Dbins - 1;
      dataspace = H5Screate_simple(1, dim1d, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "LogDensityBins_MEANBARYONDENSITY",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Ddivid);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      /* Actual data - NSF */
      dataspace = H5Screate_simple(3, dim3d, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassDistribution_NonStarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mass_bins_nsf);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      dataspace = H5Screate_simple(2, dim2d_metallicity, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MetalMassDistribution_NonStarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status =
	H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_metal_mass_bins_nsf);
      H5Dclose(dataset);
      H5Sclose(dataspace);
#ifdef BG_Z_WEIGHTED_REDSHIFT
      dataspace = H5Screate_simple(2, dim2d_metallicity, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MetallicityWeightedRedshift_NonStarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status =
	H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_z_weighted_redshift_nsf);
      H5Dclose(dataset);
      H5Sclose(dataspace);
#endif
      dataspace = H5Screate_simple(2, dim2d_dens, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassWeightedMeanDensity_NonStarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mean_dens_nsf);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      dataspace = H5Screate_simple(1, dim1d_temp, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassWeightedMeanDensity_TBins_NonStarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mean_dens_temp_nsf);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      dataspace = H5Screate_simple(1, dim1d_metallicity, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassWeightedMeanDensity_ZBins_NonStarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mean_dens_met_nsf);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      dataspace = H5Screate_simple(2, dim2d_temp, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassWeightedMeanTemperature_NonStarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mean_temp_nsf);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      dataspace = H5Screate_simple(1, dim1d_dens, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassWeightedMeanTemperature_DBins_NonStarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mean_temp_dens_nsf);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      dataspace = H5Screate_simple(1, dim1d_metallicity, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassWeightedMeanTemperature_ZBins_NonStarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mean_temp_met_nsf);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      dataspace = H5Screate_simple(2, dim2d_metallicity, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassWeightedMeanMetallicity_NonStarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status =
	H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mean_metallicity_nsf);
      H5Dclose(dataset);
      H5Sclose(dataspace);
#ifdef BG_SNIA_IRON
      dataspace = H5Screate_simple(3, dim3d, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "IronFromSNIa_NonStarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_snia_iron_bins_nsf);
      H5Dclose(dataset);
      H5Sclose(dataspace);
#endif

      /* ======================================================================== */

      /* ElementMassDistribution_NonStarFormingParticles */
      hdf5_elementgrp = H5Gcreate(hdf5_file, "/Grid/ElementMassDistribution_NonStarFormingParticles", 0);

      for(l = 0; l < BG_NELEMENTS; l++)
	{
	  for(i = 0; i < N_Zbins; i++)
	    for(j = 0; j < N_Tbins; j++)
	      for(k = 0; k < N_Dbins; k++)
		this_element[i][j][k] = global_el_mass_bins_nsf[i][j][k][l];

	  dataspace = H5Screate(H5S_SIMPLE);
	  H5Sset_extent_simple(dataspace, 3, dim3d, NULL);
	  dataset =
	    H5Dcreate(hdf5_elementgrp, ElementNames[l], H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
	  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this_element);
	  H5Dclose(dataset);
	  H5Sclose(dataspace);
	}

      H5Gclose(hdf5_elementgrp);

      /* ======================================================================== */

      /* Actual data - SF */
      dataspace = H5Screate_simple(3, dim3d, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassDistribution_StarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mass_bins_sf);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      dataspace = H5Screate_simple(2, dim2d_metallicity, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MetalMassDistribution_StarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_metal_mass_bins_sf);
      H5Dclose(dataset);
      H5Sclose(dataspace);
#ifdef BG_Z_WEIGHTED_REDSHIFT
      dataspace = H5Screate_simple(2, dim2d_metallicity, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MetallicityWeightedRedshift_StarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status =
	H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_z_weighted_redshift_sf);
      H5Dclose(dataset);
      H5Sclose(dataspace);
#endif
      dataspace = H5Screate_simple(2, dim2d_dens, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassWeightedMeanDensity_StarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mean_dens_sf);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      dataspace = H5Screate_simple(1, dim1d_temp, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassWeightedMeanDensity_TBins_StarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mean_dens_temp_sf);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      dataspace = H5Screate_simple(1, dim1d_metallicity, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassWeightedMeanDensity_ZBins_StarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mean_dens_met_sf);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      dataspace = H5Screate_simple(2, dim2d_temp, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassWeightedMeanTemperature_StarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mean_temp_sf);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      dataspace = H5Screate_simple(1, dim1d_dens, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassWeightedMeanTemperature_DBins_StarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mean_temp_dens_sf);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      dataspace = H5Screate_simple(1, dim1d_metallicity, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassWeightedMeanTemperature_ZBins_StarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mean_temp_met_sf);
      H5Dclose(dataset);
      H5Sclose(dataspace);

      dataspace = H5Screate_simple(2, dim2d_metallicity, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "MassWeightedMeanMetallicity_StarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status =
	H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_mean_metallicity_sf);
      H5Dclose(dataset);
      H5Sclose(dataspace);
#ifdef BG_SNIA_IRON
      dataspace = H5Screate_simple(3, dim3d, NULL);
      dataset = H5Dcreate(hdf5_gridgrp, "IronFromSNIa_StarFormingParticles",
			  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_snia_iron_bins_sf);
      H5Dclose(dataset);
      H5Sclose(dataspace);
#endif

      /* ======================================================================== */

      /* ElementMassDistribution_StarFormingParticles */
      hdf5_elementgrp = H5Gcreate(hdf5_file, "/Grid/ElementMassDistribution_StarFormingParticles", 0);

      for(l = 0; l < BG_NELEMENTS; l++)
	{
	  for(i = 0; i < N_Zbins; i++)
	    for(j = 0; j < N_Tbins; j++)
	      for(k = 0; k < N_Dbins; k++)
		this_element[i][j][k] = global_el_mass_bins_sf[i][j][k][l];

	  dataspace = H5Screate(H5S_SIMPLE);
	  H5Sset_extent_simple(dataspace, 3, dim3d, NULL);
	  dataset =
	    H5Dcreate(hdf5_elementgrp, ElementNames[l], H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
	  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this_element);
	  H5Dclose(dataset);
	  H5Sclose(dataspace);
	}

      H5Gclose(hdf5_elementgrp);

      /* ======================================================================== */

      H5Gclose(hdf5_gridgrp);

      H5Fclose(hdf5_file);
    }

#ifdef USE_HDF5_FIX
  H5close();
  hdf5_memory_cleanup();
#endif


}

#endif
