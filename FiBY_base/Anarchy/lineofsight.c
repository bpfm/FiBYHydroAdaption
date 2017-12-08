#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <sys/stat.h>
#include <sys/types.h>

#if defined(BG_SFR) || defined(BG_COOLING)

#include "allvars.h"
#include "proto.h"
#include "bg_vars.h"
#include "bg_proto.h"
#include "bg_vars.h"

#ifdef HAVE_HDF5

#include <hdf5.h>
#ifndef DOUBLEPRECISION		/* default is single-precision */
#define H5T_FLOAT H5T_NATIVE_FLOAT
#endif

#if (DOUBLEPRECISION == 1)	/* everything double-precision */
#define H5T_FLOAT H5T_NATIVE_DOUBLE
#endif

#if (DOUBLEPRECISION == 2)	/* mixed precision */
#define H5T_FLOAT H5T_NATIVE_FLOAT
#endif

#endif /* HAVE_HDF5 */

#ifdef OUTPUTLINEOFSIGHT

#ifdef OUTPUTLINEOFSIGHT_SPECTRUM
#define  PIXELS 512		/* number of bins along line of sight */
#else
#define  PIXELS 1
#endif

#define  N_LOS  100		/* number of lines of sight selected  */
#define  number_of_records 128

#ifdef MyMaxLosBuffer
#define MaxLosBuffer MyMaxLosBuffer	/* maximum size of los buffer */
#else
#define MaxLosBuffer -1
#endif
#define InLine       -99	/* label sph particles in los by putting itype = InLine */

static double H_a, Wmax;


struct line_of_sight
{
  int xaxis, yaxis, zaxis;
  double Xpos, Ypos;		/* relative position of line-of-sight on face of box */
  double BoxSize, Wmax, Time;

  /* total gas density */
  double Rho[PIXELS];
  double Vpec[PIXELS];
  double Temp[PIXELS];
  double Metallicity[PIXELS];

  /* neutral hydrogen */
  double RhoHI[PIXELS];
  double NHI[PIXELS];
  double VpecHI[PIXELS];
  double TempHI[PIXELS];
  double TauHI[PIXELS];

  /* HeII quantities */
  double RhoHeII[PIXELS];
  double NHeII[PIXELS];
  double VpecHeII[PIXELS];
  double TempHeII[PIXELS];
  double TauHeII[PIXELS];
}
 *Los, *LosGlobal;


/* if you change this structure you must change the definition
   of data space in memory (L. 683)!!!!!!!!!!!!!!!!!!!!!!!!!!! */
#ifdef MyMaxLosBuffer
struct line_of_sight_particles
{
  MyFloat Pos[3];
  MyFloat Vel[3];
  MyFloat Hsml;
  MyFloat Utherm;
  MyFloat Temperature;
  MyFloat Mass;
  MyFloat Density;
  MyFloat Sfr;
#ifdef BG_STELLAR_EVOLUTION
  MyFloat Metallicity;
#ifdef BG_SNIA_IRON
  MyFloat IronFromSNIa;
#endif
  MyFloat Metals[BG_NELEMENTS];
#else
  MyFloat Metallicity;
#endif
}
particles[MaxLosBuffer];
#else
struct line_of_sight_particles
{
  MyFloat Pos[3];
  MyFloat Vel[3];
  MyFloat Hsml;
  MyFloat Utherm;
  MyFloat Temperature;
  MyFloat Mass;
  MyFloat Density;
  MyFloat Sfr;
#ifdef BG_STELLAR_EVOLUTION
  MyFloat Metallicity;
#ifdef BG_SNIA_IRON
  MyFloat IronFromSNIa;
#endif
  MyFloat Metals[BG_NELEMENTS];
#else
  MyFloat Metallicity;
#endif
}
 *particles;
#endif

#ifdef HAVE_HDF5
void write_hdf5_lineofsight_particles(struct line_of_sight_particles *particles, int pc, int pcsum,
				      int count_total, hid_t hdf5_dset[number_of_records]);
void write_los_properties_in_hdf5(hid_t hdf5_particle_grp, int count_total);

hid_t element_dset[BG_NELEMENTS];

hid_t selement_dset[BG_NELEMENTS];
#endif


void extract_lineofsight_particles();

void lineofsight_output(void)
{
  char buf[500];

  int n, s;

#ifdef MEMDEBUG
  size_t nbefore;

  nbefore = check_for_largest_memory_block();
  printf("Task=%d can alloacte maximal %g MB before los-output (ceiling is %g MB)\n", ThisTask,
	 nbefore / (1024.0 * 1024.0), (nbefore + AllocatedBytes) / (1024.0 * 1024.0));
#endif


  H_a =
    All.Omega0 / (All.Time * All.Time * All.Time) + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time)
#ifdef DARKENERGY
    + DarkEnergy_a(a);
#else
    + All.OmegaLambda;
#endif
  H_a = All.Hubble * sqrt(H_a);
  Wmax = All.Time * H_a * All.BoxSize;


  if(ThisTask == 0)
    {
      sprintf(buf, "%s/los", All.OutputDir);
      mkdir(buf, 02755);
    }

  Los = (struct line_of_sight *) mymalloc(sizeof(struct line_of_sight));

  LosGlobal = (struct line_of_sight *) mymalloc(sizeof(struct line_of_sight));

  for(n = 0, s = 0; n < N_LOS; n++)
    {
      if(s + 3 >= RNDTABLE)
	{
	  set_random_numbers();
	  s = 0;
	}

      Los->zaxis = (int) (3.0 * get_random_number(s++));
      switch (Los->zaxis)
	{
	case 2:
	  Los->xaxis = 0;
	  Los->yaxis = 1;
	  break;
	case 0:
	  Los->xaxis = 1;
	  Los->yaxis = 2;
	  break;
	case 1:
	  Los->xaxis = 2;
	  Los->yaxis = 0;
	  break;
	}

      Los->Xpos = All.BoxSize * get_random_number(s++);
      Los->Ypos = All.BoxSize * get_random_number(s++);

#ifdef OUTPUTLINEOFSIGHT_SPECTRUM
      add_along_lines_of_sight();
      sum_over_processors_and_normalize();
      absorb_along_lines_of_sight();
      output_lines_of_sight(n);
#endif

#ifdef OUTPUTLINEOFSIGHT_PARTICLES
      find_particles_and_save_them(n);
#endif
    }

  myfree(LosGlobal);
  myfree(Los);

#ifdef USE_HDF5_FIX
  H5close();
  hdf5_memory_cleanup();
#endif

#ifdef MEMDEBUG
  size_t nafter;

  nafter = check_for_largest_memory_block();
  if(nafter < nbefore)
    printf("Task=%d can allocate %g MB after los-output, LESS than before (%g)\n", ThisTask,
	   nafter / (1024.0 * 1024.0), nbefore / (1024.0 * 1024.0));
#endif
}

void extract_lineofsight_particles()
	  /* extract particle in sightline into 'particles' array */
	  /* exit loop if MaxLosBuffer extracted, if MaxLosBuffer > 0 */
{
  int n, nsent = 0;

  double a3inv, ascale;

  double rho, redshift, He_frac, eint, n_H, XH;

  int j, k, z_index;

  float d_z;

  static int Hydrogen_index, Helium_index;

  static int first_call = 0;

  if(first_call == 0)
    {
      first_call = 1;
      if(ThisTask == 0)
	{
	  printf(" sightline buffer has maximum size  = %d\n", MaxLosBuffer);
	}
#ifdef BG_COOLING
      Hydrogen_index = element_index("Hydrogen");
      Helium_index = element_index("Helium");
#endif
    }

  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      ascale = All.Time;
      redshift = 1 / ascale - 1;
    }
  else
    {
      a3inv = ascale = 1;
      redshift = 0;
    }

  /* get d_z for u to temp conversion */
  get_redshift_index(redshift, &z_index, &d_z);

  for(n = 0; n < N_gas; n++)
    {
      if(P[n].Type == InLine)
	{
	  P[n].Type = 0;	/* reset to type gas */

	  particles[nsent].Mass = P[n].Mass;
	  particles[nsent].Density = SphP[n].d.Density;
	  particles[nsent].Sfr = SphP[n].Sfr;

	  for(k = 0; k < 3; k++)
	    particles[nsent].Pos[k] = P[n].Pos[k];

	  for(k = 0; k < 3; k++)
	    particles[nsent].Vel[k] = P[n].Vel[k];

	  particles[nsent].Hsml = PPP[n].Hsml;
	  particles[nsent].Utherm = SphP[n].Entropy / GAMMA_MINUS1 *
	    pow(SphP[n].d.Density * a3inv, GAMMA_MINUS1);

#ifdef BG_STELLAR_EVOLUTION
	  for(j = 0; j < BG_NELEMENTS; j++)
	    particles[nsent].Metals[j] = SphP[n].Metals[j];

	  particles[nsent].Metallicity = SphP[n].Metallicity;
#ifdef BG_SNIA_IRON
	  particles[nsent].IronFromSNIa = SphP[n].IronFromSNIa;
#endif
#else
	  particles[nsent].Metallicity = P[n].Metallicity;
#endif

	  /* Convert to metal fraction */
#ifdef BG_STELLAR_EVOLUTION
	  for(j = 0; j < BG_NELEMENTS; j++)
	    particles[nsent].Metals[j] /= particles[nsent].Mass;

#ifdef BG_SNIA_IRON
	  particles[nsent].IronFromSNIa /= particles[nsent].Mass;
#endif
#endif

#ifdef BG_STELLAR_EVOLUTION
	  XH = particles[nsent].Metals[Hydrogen_index];
	  He_frac = particles[nsent].Metals[Helium_index];
#else
	  XH = 0.25;
	  He_frac = 1. - XH;
#endif

	  eint = SphP[n].Entropy / GAMMA_MINUS1 * pow(SphP[n].d.Density * a3inv, GAMMA_MINUS1);
	  eint *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

	  rho = SphP[n].d.Density * a3inv;
	  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

	  /* convert Hydrogen mass fraction in Hydrogen number density */
	  n_H = rho * XH / PROTONMASS;

	  particles[nsent].Temperature = convert_u_to_temp(d_z, eint, n_H, He_frac);

	  nsent++;
	  if(MaxLosBuffer > 0)
	    {
	      if(nsent >= MaxLosBuffer)
		break;
	    }
	}
    }
  if(nsent > 0)
    nsent--;
}

void find_particles_and_save_them(int num)
{
  int n, count_local, this_count, maxcount, *countlist, extracted, counttot, rep, idata, ndata = 0;

  MyFloat dx, dy, r2;

  char fname[1000];

  MPI_Status status;

  FILE *fd = 0;


#ifdef HAVE_HDF5
  static hid_t hdf5_file = 0;

  hid_t hdf5_headergrp = 0, hdf5_paramgrp = 0;

  hid_t hdf5_constgrp = 0, hdf5_unitsgrp = 0;

  hid_t hdf5_particle_grp = 0, hdf5_dataspace = 0, hdf5_attribute = 0;

  hid_t hdf5_dset[number_of_records], element_grp = 0;

  hid_t hdf5_dataspace_in_file[number_of_records];

  hsize_t dims[2];

  int rank = 0, pcsum = 0, pc = 0;

  char buf[500], part_grp_name[500];

  int i_los, count_total = 0;

  enum iofields blocknr;
#endif

  double a3inv, ascale;

  double redshift;

  int i, j, z_index;

  float d_z;

  int nthis_send, nsent;

  static int first_call = 0;

  static int Hydrogen_index, Helium_index;

  if(first_call == 0)
    {
      first_call = 1;
#ifdef BG_COOLING
      Hydrogen_index = element_index("Hydrogen");
      Helium_index = element_index("Helium");
#endif
    }

#ifdef HAVE_HDF5
  for(i = 0; i < number_of_records; i++)
    hdf5_dset[i] = 0;
#endif

  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      ascale = All.Time;
      redshift = 1 / ascale - 1;
    }
  else
    {
      a3inv = ascale = 1;
      redshift = 0;
    }

  /* update header informations */
  //  header.time = All.Time;
  //  header.redshift = redshift;





  /* determine particle numbers of each type in file */

/*   if(ThisTask == 0) */
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






  /* get d_z for u to temp conversion */
  get_redshift_index(redshift, &z_index, &d_z);

  countlist = (int *) mymalloc(sizeof(int) * NTask);


  for(n = 0, count_local = 0; n < N_gas; n++)
    {
      if(P[n].Type == 0)
	{
	  dx = los_periodic(P[n].Pos[Los->xaxis] - Los->Xpos);

	  if(dx < PPP[n].Hsml)
	    {
	      dy = los_periodic(P[n].Pos[Los->yaxis] - Los->Ypos);

	      if(dy < PPP[n].Hsml)
		{
		  r2 = dx * dx + dy * dy;

		  if(r2 < PPP[n].Hsml * PPP[n].Hsml)
		    {
		      P[n].Type = InLine;

		      count_local++;
		    }
		}
	    }
	}
    }


  MPI_Gather(&count_local, 1, MPI_INT, countlist, 1, MPI_INT, 0, MPI_COMM_WORLD);

#ifndef MyMaxLosBuffer
  if(MaxLosBuffer > 0)
    particles = mymalloc(sizeof(struct line_of_sight_particles) * MaxLosBuffer);
  else
    {
      maxcount = count_local;
      if(ThisTask == 0)
	for(rep = 0; rep < NTask; rep++)
	  if(maxcount < countlist[rep])
	    maxcount = countlist[rep];
      if(ThisTask == 0)
	printf(" allocation size = %d\n", maxcount);
      particles = mymalloc(sizeof(struct line_of_sight_particles) * maxcount);
    }
#endif

#ifdef HAVE_HDF5
  /* total number of particles in this sightline */
  if(ThisTask == 0)
    {
      for(n = 0, count_total = 0; n < NTask; n++)
	count_total += countlist[n];
    }
#endif

  if(ThisTask == 0)
    {
      sprintf(fname, "%s/los/part_los_z%05.3f_%03d.dat", All.OutputDir, 1 / All.Time - 1, num);

      if(All.SnapFormat == 3)
	{
	  sprintf(fname, "%s/los/part_los_z%05.3f", All.OutputDir, 1 / All.Time - 1);
	  sprintf(buf, "%s.hdf5", fname);

	  if(num == 0)
	    {
	      hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	      hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);
	      hdf5_paramgrp = H5Gcreate(hdf5_file, "/Parameters", 0);
	      hdf5_unitsgrp = H5Gcreate(hdf5_file, "/Units", 0);
	      hdf5_constgrp = H5Gcreate(hdf5_file, "/Constants", 0);

	      /* usual header information */
	      write_header_attributes_in_hdf5(hdf5_headergrp);
	      write_parameters_attributes_in_hdf5(hdf5_paramgrp);
	      write_units_attributes_in_hdf5(hdf5_unitsgrp);
	      write_constants_attributes_in_hdf5(hdf5_constgrp);

	      /* number of sight lines in this file */
	      hdf5_dataspace = H5Screate(H5S_SCALAR);
	      hdf5_attribute = H5Acreate(hdf5_headergrp, "Number_of_sight_lines",
					 H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);

	      i_los = N_LOS;
	      H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &i_los);
	      H5Aclose(hdf5_attribute);
	      H5Sclose(hdf5_dataspace);
	      H5Gclose(hdf5_constgrp);
	      H5Gclose(hdf5_unitsgrp);
	      H5Gclose(hdf5_headergrp);
	      H5Gclose(hdf5_paramgrp);
	    }
	}
      else
	{
	  if(!(fd = fopen(fname, "w")))
	    {
	      printf("can't open file `%s`\n", fname);
	      endrun(112);
	    }

	  my_fwrite(&count_local, sizeof(int), 1, fd);	/* will be overwritten later */
	  my_fwrite(&LosGlobal->xaxis, sizeof(int), 1, fd);
	  my_fwrite(&LosGlobal->yaxis, sizeof(int), 1, fd);
	  my_fwrite(&LosGlobal->zaxis, sizeof(int), 1, fd);
	  my_fwrite(&LosGlobal->Xpos, sizeof(double), 1, fd);
	  my_fwrite(&LosGlobal->Ypos, sizeof(double), 1, fd);
	  my_fwrite(&LosGlobal->BoxSize, sizeof(double), 1, fd);
	  my_fwrite(&LosGlobal->Wmax, sizeof(double), 1, fd);
	  my_fwrite(&LosGlobal->Time, sizeof(double), 1, fd);
	}
    }

#ifdef HAVE_HDF5
  if(ThisTask == 0)
    {
      /* create groups for particle properties */
      sprintf(part_grp_name, "/LOS%d", num);
      hdf5_particle_grp = H5Gcreate(hdf5_file, part_grp_name, 0);

      /* write los projection properties */
      write_los_properties_in_hdf5(hdf5_particle_grp, count_total);
      /* write los properties */
      ndata = 0;

      /* positions */
      dims[0] = count_total;
      dims[1] = 3;
      rank = 2;
      hdf5_dataspace_in_file[ndata] = H5Screate_simple(rank, dims, NULL);
      hdf5_dset[ndata] = H5Dcreate(hdf5_particle_grp, "Positions", H5T_FLOAT,
				   hdf5_dataspace_in_file[ndata], H5P_DEFAULT);
      blocknr = IO_POS;
      write_attributes_in_hdf5(blocknr, hdf5_dset[ndata]);
      ndata += 1;


      /* velocity */
      dims[0] = count_total;
      dims[1] = 3;
      rank = 2;
      hdf5_dataspace_in_file[ndata] = H5Screate_simple(rank, dims, NULL);
      hdf5_dset[ndata] = H5Dcreate(hdf5_particle_grp, "Velocity", H5T_FLOAT,
				   hdf5_dataspace_in_file[ndata], H5P_DEFAULT);
      blocknr = IO_VEL;
      write_attributes_in_hdf5(blocknr, hdf5_dset[ndata]);
      ndata += 1;


      /* hsml */
      dims[0] = count_total;
      rank = 1;
      hdf5_dataspace_in_file[ndata] = H5Screate_simple(rank, dims, NULL);
      hdf5_dset[ndata] = H5Dcreate(hdf5_particle_grp, "SmoothingLength", H5T_FLOAT,
				   hdf5_dataspace_in_file[ndata], H5P_DEFAULT);
      blocknr = IO_HSML;
      write_attributes_in_hdf5(blocknr, hdf5_dset[ndata]);
      ndata += 1;


      /* thermal energy */
      dims[0] = count_total;
      rank = 1;
      hdf5_dataspace_in_file[ndata] = H5Screate_simple(rank, dims, NULL);
      hdf5_dset[ndata] = H5Dcreate(hdf5_particle_grp, "InternalEnergy", H5T_FLOAT,
				   hdf5_dataspace_in_file[ndata], H5P_DEFAULT);
      blocknr = IO_U;
      write_attributes_in_hdf5(blocknr, hdf5_dset[ndata]);
      ndata += 1;


      /* temperature */
      dims[0] = count_total;
      rank = 1;
      hdf5_dataspace_in_file[ndata] = H5Screate_simple(rank, dims, NULL);
      hdf5_dset[ndata] = H5Dcreate(hdf5_particle_grp, "Temperature", H5T_FLOAT,
				   hdf5_dataspace_in_file[ndata], H5P_DEFAULT);
      blocknr = IO_BG_TEMP;
      write_attributes_in_hdf5(blocknr, hdf5_dset[ndata]);
      ndata += 1;


      /* mass */
      dims[0] = count_total;
      rank = 1;
      hdf5_dataspace_in_file[ndata] = H5Screate_simple(rank, dims, NULL);
      hdf5_dset[ndata] = H5Dcreate(hdf5_particle_grp, "Mass", H5T_FLOAT,
				   hdf5_dataspace_in_file[ndata], H5P_DEFAULT);
      blocknr = IO_MASS;
      write_attributes_in_hdf5(blocknr, hdf5_dset[ndata]);
      ndata += 1;


      /* density */
      dims[0] = count_total;
      rank = 1;
      hdf5_dataspace_in_file[ndata] = H5Screate_simple(rank, dims, NULL);
      hdf5_dset[ndata] = H5Dcreate(hdf5_particle_grp, "Density", H5T_FLOAT,
				   hdf5_dataspace_in_file[ndata], H5P_DEFAULT);
      blocknr = IO_RHO;
      write_attributes_in_hdf5(blocknr, hdf5_dset[ndata]);
      ndata += 1;

      /* star formation rate */
      dims[0] = count_total;
      rank = 1;
      hdf5_dataspace_in_file[ndata] = H5Screate_simple(rank, dims, NULL);
      hdf5_dset[ndata] = H5Dcreate(hdf5_particle_grp, "StarFormationRate", H5T_FLOAT,
				   hdf5_dataspace_in_file[ndata], H5P_DEFAULT);
      blocknr = IO_SFR;
      write_attributes_in_hdf5(blocknr, hdf5_dset[ndata]);
      ndata += 1;


#ifdef BG_STELLAR_EVOLUTION
      /* Metallicity */
      dims[0] = count_total;
      rank = 1;
      hdf5_dataspace_in_file[ndata] = H5Screate_simple(rank, dims, NULL);
      hdf5_dset[ndata] = H5Dcreate(hdf5_particle_grp, "Metallicity", H5T_FLOAT,
				   hdf5_dataspace_in_file[ndata], H5P_DEFAULT);
      blocknr = IO_BG_METALLICITY;
      write_attributes_in_hdf5(blocknr, hdf5_dset[ndata]);
      ndata += 1;


#ifdef BG_SNIA_IRON
      dims[0] = count_total;
      rank = 1;
      hdf5_dataspace_in_file[ndata] = H5Screate_simple(rank, dims, NULL);
      hdf5_dset[ndata] = H5Dcreate(hdf5_particle_grp, "IronFromSNIa", H5T_FLOAT,
				   hdf5_dataspace_in_file[ndata], H5P_DEFAULT);
      blocknr = IO_BG_IRON_FROM_SNIA;
      write_attributes_in_hdf5(blocknr, hdf5_dset[ndata]);
      ndata += 1;

#endif

      /* data sets for metals */
      for(j = 0; j < BG_NELEMENTS; j++)
	{
	  element_dset[j] = -1;
	  selement_dset[j] = -1;
	}

      /* Metals */
      element_grp = H5Gcreate(hdf5_particle_grp, "ElementAbundance", 0);
      rank = 1;
      dims[0] = count_total;
      dims[1] = 1;

      blocknr = IO_BG_METALS;
      write_attributes_in_hdf5(blocknr, element_grp);

      for(j = 0; j < BG_NELEMENTS; j++)
	{
	  hdf5_dataspace_in_file[ndata] = H5Screate_simple(rank, dims, NULL);
	  hdf5_dset[ndata] = H5Dcreate(element_grp, ElementNames[j], H5T_FLOAT,
				       hdf5_dataspace_in_file[ndata], H5P_DEFAULT);
	  element_dset[j] = hdf5_dset[ndata];
	  ndata += 1;
	}

#else /* BG_STELLAR_EVOLUTION */
      /* Metallicity */
      dims[0] = count_total;
      rank = 1;
      hdf5_dataspace_in_file[ndata] = H5Screate_simple(rank, dims, NULL);
      hdf5_dset[ndata] = H5Dcreate(hdf5_particle_grp, "Metallicity", H5T_FLOAT,
				   hdf5_dataspace_in_file[ndata], H5P_DEFAULT);
      ndata += 1;

#endif /* BG_STELLAR_EVOLUTION */

      /* check whether allocation worked ok */
      for(i = 0; i < ndata; i++)
	if(hdf5_dset[i] < 0)
	  {
	    printf(" line of sight output: error allocating record number %d\n", i);
	    endrun(i);
	  }
      pcsum = 0;		/* cumulative number of particles written */
    }

#endif
  for(rep = 0, counttot = 0; rep < NTask; rep++)
    {
      extracted = 0;
      /* load particles */
      if(rep == 0 && extracted < count_local)
	extract_lineofsight_particles();


      /* Tasks 0 and rep exchange data */
      if(rep == ThisTask || ThisTask == 0)
	{
	  if(ThisTask > 0)
	    this_count = count_local;
	  else
	    this_count = countlist[rep];

	  nsent = 0;
	  while(nsent < this_count)
	    {
	      if(MaxLosBuffer > 0)
		{
		  if(extracted + MaxLosBuffer > this_count)
		    nthis_send = this_count - extracted;
		  else
		    nthis_send = MaxLosBuffer;
		}
	      else
		nthis_send = this_count;

	      extracted += nthis_send;
	      nsent += nthis_send;

	      if(ThisTask != 0 && rep == ThisTask && nthis_send > 0)
		MPI_Ssend(particles, sizeof(struct line_of_sight_particles) * nthis_send, MPI_BYTE, 0,
			  TAG_PDATA, MPI_COMM_WORLD);

	      if(ThisTask == 0)
		{
		  if(rep > 0 && nthis_send > 0)
		    MPI_Recv(particles, sizeof(struct line_of_sight_particles) * nthis_send,
			     MPI_BYTE, rep, TAG_PDATA, MPI_COMM_WORLD, &status);

#ifdef HAVE_HDF5
		  pc = nthis_send;
		  write_hdf5_lineofsight_particles(particles, pc, pcsum, count_total, hdf5_dset);
		  pcsum += pc;
#else
		  my_fwrite(particles, sizeof(struct line_of_sight_particles), nthis_send, fd);
#endif
		  counttot += nthis_send;
		}
	      if(ThisTask == rep && extracted < this_count)
		extract_lineofsight_particles();
	    }
	}
    }

  if(ThisTask == 0)
    {
#ifdef HAVE_HDF5
      /* close the data sets and data spaces */
      for(idata = 0; idata < ndata; idata++)
	{
	  H5Dclose(hdf5_dset[idata]);
	  H5Sclose(hdf5_dataspace_in_file[idata]);
	}

      /* close the groups */
#ifdef BG_STELLAR_EVOLUTION
      H5Gclose(element_grp);
#endif
/* #ifdef BG_METALSMOOTHING */
/*       H5Gclose(selement_grp); */
/* #endif */
      H5Gclose(hdf5_particle_grp);

      /* close the file if last LOS to write */
      if(num == N_LOS - 1)
	H5Fclose(hdf5_file);
#else
      fclose(fd);
      if(!(fd = fopen(fname, "r+")))
	{
	  printf("can't open file `%s'\n", fname);
	  endrun(113);
	}

      fseek(fd, 0, SEEK_CUR);
      my_fwrite(&counttot, sizeof(int), 1, fd);
      fclose(fd);
#endif
    }

#ifndef MyMaxLosBuffer
  myfree(particles);
#endif
  myfree(countlist);
}


#ifdef HAVE_HDF5
void write_los_properties_in_hdf5(hid_t hdf5_particle_grp, int count_total)
	  /* write global properties of this sightline to file */
{
  hid_t hdf5_dataspace = 0, hdf5_attribute = 0;

  /* number of particles */
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_particle_grp, "Number_of_part_this_los",
			     H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &count_total);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  /* determine whether sightline is parallel to x, y or z */
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_particle_grp, "x-axis", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &Los->xaxis);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_particle_grp, "y-axis", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &Los->yaxis);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_particle_grp, "z-axis", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &Los->zaxis);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  /* x, y-position of sightline  */
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_particle_grp, "x-position", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &Los->Xpos);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_particle_grp, "y-position", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &Los->Ypos);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
}

void write_hdf5_lineofsight_particles(struct line_of_sight_particles *particles, int pc, int pcsum,
				      int count_total, hid_t hdf5_dset[number_of_records])
{
  hid_t hdf5_dataspace_memory, hdf5_dataspace_metals = 0;

  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0;

  herr_t hdf5_status;

  hsize_t dims[2], count[2], start[2];

  /* */
  int rank = 0, isize, istart, idata;

  int i, j;

  float *this_metal;

  /* skip if nothing to write */
  if(pc > 0)
    {
      istart = 0;
      idata = 0;

      /* define data space in memory (structure "particles") */
      rank = 2;
      dims[0] = pc;
      dims[1] = 12;
#ifdef BG_STELLAR_EVOLUTION
      dims[1] += BG_NELEMENTS + 1;
#ifdef BG_SNIA_IRON
      dims[1] += 1;
#endif
#else
      dims[1] += 1;
#endif
      hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);

      rank = 1;
      dims[0] = pc;
      hdf5_dataspace_metals = H5Screate_simple(rank, dims, NULL);

      /* ---------------------------------------------------------------- */
      /* positions and velocity: 2D arrays!                                */
      /* ---------------------------------------------------------------- */

      isize = 3;

      /* define data space in file, and select hyperslab in file */
      rank = 2;
      dims[0] = count_total;
      dims[1] = 3;
      hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
      start[0] = pcsum;
      start[1] = 0;
      count[0] = pc;
      count[1] = 3;
      H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

      /* set the data type for float values */
      hdf5_datatype = H5Tcopy(H5T_FLOAT);

      /* select hyperslab in memory */
      start[0] = 0;
      start[1] = istart;
      istart += isize;
      count[0] = pc;
      count[1] = isize;

      /* positions */
      H5Sselect_hyperslab(hdf5_dataspace_memory, H5S_SELECT_SET, start, NULL, count, NULL);
      hdf5_status = H5Dwrite(hdf5_dset[idata], hdf5_datatype, hdf5_dataspace_memory,
			     hdf5_dataspace_in_file, H5P_DEFAULT, particles);
      idata += 1;

      /* velocity */
      start[1] = istart;
      istart += isize;
      H5Sselect_hyperslab(hdf5_dataspace_memory, H5S_SELECT_SET, start, NULL, count, NULL);
      hdf5_status = H5Dwrite(hdf5_dset[idata], hdf5_datatype, hdf5_dataspace_memory,
			     hdf5_dataspace_in_file, H5P_DEFAULT, particles);
      idata += 1;

      H5Sclose(hdf5_dataspace_in_file);

      /* ---------------------------------------------------------------- */
      /* all the other arrays are 1D                                      */
      /* ---------------------------------------------------------------- */

      isize = 1;

      /* define data space in file, and select hyperslab in file */
      rank = 1;
      dims[0] = count_total;
      hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
      start[0] = pcsum;
      count[0] = pc;
      H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

      /* select hyperslab in memory */
      start[0] = 0;
      count[0] = pc;
      count[1] = isize;

      /* hsml */
      start[1] = istart;
      istart += isize;
      H5Sselect_hyperslab(hdf5_dataspace_memory, H5S_SELECT_SET, start, NULL, count, NULL);
      hdf5_status = H5Dwrite(hdf5_dset[idata], hdf5_datatype, hdf5_dataspace_memory,
			     hdf5_dataspace_in_file, H5P_DEFAULT, particles);
      idata += 1;

      /* thermal energy */
      start[1] = istart;
      istart += isize;
      H5Sselect_hyperslab(hdf5_dataspace_memory, H5S_SELECT_SET, start, NULL, count, NULL);
      hdf5_status = H5Dwrite(hdf5_dset[idata], hdf5_datatype, hdf5_dataspace_memory,
			     hdf5_dataspace_in_file, H5P_DEFAULT, particles);
      idata += 1;

      /* temperature */
      start[1] = istart;
      istart += isize;
      H5Sselect_hyperslab(hdf5_dataspace_memory, H5S_SELECT_SET, start, NULL, count, NULL);
      hdf5_status = H5Dwrite(hdf5_dset[idata], hdf5_datatype, hdf5_dataspace_memory,
			     hdf5_dataspace_in_file, H5P_DEFAULT, particles);
      idata += 1;

      /* mass */
      start[1] = istart;
      istart += isize;
      H5Sselect_hyperslab(hdf5_dataspace_memory, H5S_SELECT_SET, start, NULL, count, NULL);
      hdf5_status = H5Dwrite(hdf5_dset[idata], hdf5_datatype, hdf5_dataspace_memory,
			     hdf5_dataspace_in_file, H5P_DEFAULT, particles);
      idata += 1;

      /* density */
      start[1] = istart;
      istart += isize;
      H5Sselect_hyperslab(hdf5_dataspace_memory, H5S_SELECT_SET, start, NULL, count, NULL);
      hdf5_status = H5Dwrite(hdf5_dset[idata], hdf5_datatype, hdf5_dataspace_memory,
			     hdf5_dataspace_in_file, H5P_DEFAULT, particles);
      idata += 1;

      /* star formation rate */
      start[1] = istart;
      istart += isize;
      H5Sselect_hyperslab(hdf5_dataspace_memory, H5S_SELECT_SET, start, NULL, count, NULL);
      hdf5_status = H5Dwrite(hdf5_dset[idata], hdf5_datatype, hdf5_dataspace_memory,
			     hdf5_dataspace_in_file, H5P_DEFAULT, particles);
      idata += 1;

#ifdef BG_STELLAR_EVOLUTION
      /* metallicity */
      start[1] = istart;
      istart += isize;
      H5Sselect_hyperslab(hdf5_dataspace_memory, H5S_SELECT_SET, start, NULL, count, NULL);
      hdf5_status = H5Dwrite(hdf5_dset[idata], hdf5_datatype, hdf5_dataspace_memory,
			     hdf5_dataspace_in_file, H5P_DEFAULT, particles);
      idata += 1;

      /* iron from SNIa */
#ifdef BG_SNIA_IRON
      start[1] = istart;
      istart += isize;
      H5Sselect_hyperslab(hdf5_dataspace_memory, H5S_SELECT_SET, start, NULL, count, NULL);
      hdf5_status = H5Dwrite(hdf5_dset[idata], hdf5_datatype, hdf5_dataspace_memory,
			     hdf5_dataspace_in_file, H5P_DEFAULT, particles);
      idata += 1;
#endif
      /* Element abundances + hydrogen */
      this_metal = (float *) mymalloc(sizeof(float) * count[0]);

      /* write all metals */
      for(j = 0; j < BG_NELEMENTS; j++)
	{
	  for(i = 0; i < (int) count[0]; i++)
	    this_metal[i] = particles[i].Metals[j];
	  hdf5_status = H5Dwrite(element_dset[j], hdf5_datatype, hdf5_dataspace_metals,
				 hdf5_dataspace_in_file, H5P_DEFAULT, this_metal);
	  idata += 1;
	}

      myfree(this_metal);
#else /* BG_STELLAR_EVOLUTION */
      /* metallicity */
      start[1] = istart;
      istart += isize;
      H5Sselect_hyperslab(hdf5_dataspace_memory, H5S_SELECT_SET, start, NULL, count, NULL);
      hdf5_status = H5Dwrite(hdf5_dset[idata], hdf5_datatype, hdf5_dataspace_memory,
			     hdf5_dataspace_in_file, H5P_DEFAULT, particles);
      idata += 1;
#endif /* BG_STELLAR_EVOLUTION */

      H5Tclose(hdf5_datatype);

      H5Sclose(hdf5_dataspace_metals);
      H5Sclose(hdf5_dataspace_in_file);
      H5Sclose(hdf5_dataspace_memory);
    }
}
#endif




void add_along_lines_of_sight(void)
{
#ifdef OUTPUTLINEOFSIGHT_SPECTRUM
  int n, bin, i, iz0, iz1, iz;

  double dx, dy, dz, r, r2, ne, nh0, nHeII, utherm, temp, meanWeight;

  double u, wk, weight, a3inv, h3inv;

  double z0, z1, dmax1, dmax2;

  for(i = 0; i < PIXELS; i++)
    {
      Los->Rho[i] = 0;
      Los->Vpec[i] = 0;
      Los->Temp[i] = 0;
      Los->Metallicity[i] = 0;

      Los->RhoHI[i] = 0;
      Los->NHI[i] = 0;
      Los->VpecHI[i] = 0;
      Los->TempHI[i] = 0;
      Los->TauHI[i] = 0;

      Los->RhoHeII[i] = 0;
      Los->NHeII[i] = 0;
      Los->VpecHeII[i] = 0;
      Los->TempHeII[i] = 0;
      Los->TauHeII[i] = 0;
    }

  a3inv = 1.0 / (All.Time * All.Time * All.Time);

  for(n = 0; n < N_gas; n++)
    {
      if(P[n].Type == 0)
	{
	  dx = los_periodic(P[n].Pos[Los->xaxis] - Los->Xpos);
	  dy = los_periodic(P[n].Pos[Los->yaxis] - Los->Ypos);

	  r2 = dx * dx + dy * dy;

	  if(r2 < PPP[n].Hsml * PPP[n].Hsml)
	    {
	      z0 = (P[n].Pos[Los->zaxis] - PPP[n].Hsml) / All.BoxSize * PIXELS;
	      z1 = (P[n].Pos[Los->zaxis] + PPP[n].Hsml) / All.BoxSize * PIXELS;
	      iz0 = (int) z0;
	      iz1 = (int) z1;
	      if(z0 < 0)
		iz0 -= 1;

	      for(iz = iz0; iz <= iz1; iz++)
		{
		  dz = los_periodic((iz + 0.5) / PIXELS * All.BoxSize - P[n].Pos[Los->zaxis]);
		  r = sqrt(r2 + dz * dz);

		  if(PPP[n].Hsml > All.BoxSize)
		    {
		      printf("Here:%d  n=%d %g\n", ThisTask, n, PPP[n].Hsml);
		      endrun(89);
		    }

		  if(r < PPP[n].Hsml)
		    {
		      u = r / PPP[n].Hsml;
		      h3inv = 1.0 / (PPP[n].Hsml * PPP[n].Hsml * PPP[n].Hsml);

		      if(u < 0.5)
			wk = h3inv * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		      else
			wk = h3inv * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

		      bin = iz;
		      while(bin >= PIXELS)
			bin -= PIXELS;
		      while(bin < 0)
			bin += PIXELS;

		      ne = SphP[n].Ne;
		      utherm = DMAX(All.MinEgySpec,
				    SphP[n].Entropy / GAMMA_MINUS1 * pow(SphP[n].d.Density *
									 a3inv, GAMMA_MINUS1));

		      AbundanceRatios(utherm, SphP[n].d.Density * a3inv, &ne, &nh0, &nHeII);

		      meanWeight = 4.0 / (3 * HYDROGEN_MASSFRAC + 1 + 4 * HYDROGEN_MASSFRAC * ne);

		      temp = meanWeight * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * utherm
			* All.UnitEnergy_in_cgs / All.UnitMass_in_g;

		      /* do total gas */
		      weight = P[n].Mass * wk;
		      Los->Rho[bin] += weight;
		      Los->Metallicity[bin] += P[n].Metallicity * weight;
		      Los->Temp[bin] += temp * weight;
		      Los->Vpec[bin] += P[n].Vel[Los->zaxis] * weight;

		      /* do neutral hydrogen */
		      weight = nh0 * HYDROGEN_MASSFRAC * P[n].Mass * wk;
		      Los->RhoHI[bin] += weight;
		      Los->TempHI[bin] += temp * weight;
		      Los->VpecHI[bin] += P[n].Vel[Los->zaxis] * weight;

		      /* do HeII */
		      weight = 4 * nHeII * HYDROGEN_MASSFRAC * P[n].Mass * wk;
		      Los->RhoHeII[bin] += weight;
		      Los->TempHeII[bin] += temp * weight;
		      Los->VpecHeII[bin] += P[n].Vel[Los->zaxis] * weight;
		    }
		}
	    }
	}
    }
#endif /* OUTPUTLINEOFSIGHT_SPECTRUM */
}


void sum_over_processors_and_normalize(void)
{
  int bin;

  MPI_Reduce(Los->Rho, LosGlobal->Rho, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Los->Metallicity, LosGlobal->Metallicity, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Los->Temp, LosGlobal->Temp, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Los->Vpec, LosGlobal->Vpec, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(Los->RhoHI, LosGlobal->RhoHI, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Los->TempHI, LosGlobal->TempHI, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Los->VpecHI, LosGlobal->VpecHI, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(Los->RhoHeII, LosGlobal->RhoHeII, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Los->TempHeII, LosGlobal->TempHeII, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Los->VpecHeII, LosGlobal->VpecHeII, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


  if(ThisTask == 0)
    {
      /* normalize results by the weights */
      for(bin = 0; bin < PIXELS; bin++)
	{
	  /* total gas density */
	  LosGlobal->Metallicity[bin] /= LosGlobal->Rho[bin];
	  LosGlobal->Temp[bin] /= LosGlobal->Rho[bin];
	  LosGlobal->Vpec[bin] /= (All.Time * LosGlobal->Rho[bin]);

	  /* neutral hydrogen quantities */
	  LosGlobal->VpecHI[bin] /= LosGlobal->RhoHI[bin];
	  LosGlobal->TempHI[bin] /= LosGlobal->RhoHI[bin];
	  LosGlobal->NHI[bin] = LosGlobal->RhoHI[bin] * (All.UnitMass_in_g / PROTONMASS);

	  /* HeII quantities */
	  LosGlobal->VpecHeII[bin] /= (All.Time * LosGlobal->RhoHeII[bin]);
	  LosGlobal->TempHeII[bin] /= LosGlobal->RhoHeII[bin];
	  LosGlobal->NHeII[bin] = LosGlobal->RhoHeII[bin] * (All.UnitMass_in_g / (4 * PROTONMASS));
	}
    }
}



void absorb_along_lines_of_sight(void)
{
  double dz, dv, b, fac, fac_HeII;

  int bin, k;


  if(ThisTask == 0)
    {
      dz = All.BoxSize / PIXELS;

      for(bin = 0; bin < PIXELS; bin++)
	{
	  LosGlobal->TauHI[bin] = 0;
	  LosGlobal->TauHeII[bin] = 0;

	  for(k = 0; k < PIXELS; k++)
	    {
	      dv = (k - bin);

	      while(dv < -PIXELS / 2)
		dv += PIXELS;
	      while(dv > PIXELS / 2)
		dv -= PIXELS;

	      dv = (dv * Wmax / PIXELS + LosGlobal->VpecHI[k]) * All.UnitVelocity_in_cm_per_s;

	      b = sqrt(2 * BOLTZMANN * LosGlobal->TempHI[k] / PROTONMASS);

	      LosGlobal->TauHI[bin] += LosGlobal->NHI[k] * exp(-dv * dv / (b * b)) / b * dz;


	      /* now HeII */
	      dv = (k - bin);

	      while(dv < -PIXELS / 2)
		dv += PIXELS;
	      while(dv > PIXELS / 2)
		dv -= PIXELS;

	      dv = (dv * Wmax / PIXELS + LosGlobal->VpecHeII[k]) * All.UnitVelocity_in_cm_per_s;

	      b = sqrt(2 * BOLTZMANN * LosGlobal->TempHeII[k] / (4 * PROTONMASS));

	      LosGlobal->TauHeII[bin] += LosGlobal->NHeII[k] * exp(-dv * dv / (b * b)) / b * dz;
	    }
	}


      /* multiply with correct prefactors */

      /*  to get things into cgs units */
      fac = 1 / pow(All.UnitLength_in_cm, 2);

      fac *= All.HubbleParam * All.HubbleParam;

      fac *= OSCILLATOR_STRENGTH * M_PI * LYMAN_ALPHA * sqrt(3 * THOMPSON / (8 * M_PI));	/* Ly-alpha cross section */

      fac *= C / (All.Time * All.Time) / sqrt(M_PI);

      /* Note: For HeII, the oscillator strength is equal to that of HI,
         and the Lyman-alpha wavelength is 4 times shorter */

      fac_HeII = fac * (OSCILLATOR_STRENGTH_HeII / OSCILLATOR_STRENGTH) * (LYMAN_ALPHA_HeII / LYMAN_ALPHA);

      for(bin = 0; bin < PIXELS; bin++)
	{
	  LosGlobal->TauHI[bin] *= fac;
	  LosGlobal->TauHeII[bin] *= fac_HeII;
	}

      LosGlobal->BoxSize = All.BoxSize;
      LosGlobal->Wmax = Wmax;
      LosGlobal->Time = All.Time;
    }

}



void output_lines_of_sight(int num)
{
  FILE *fd;

  int dummy;

  char fname[400];

  sprintf(fname, "%s/los/spec_los_z%05.3f_%03d.dat", All.OutputDir, 1 / All.Time - 1, num);

  if(!(fd = fopen(fname, "w")))
    {
      printf("can't open file `%s`\n", fname);
      exit(1);
    }

  dummy = PIXELS;
  my_fwrite(&dummy, sizeof(int), 1, fd);
  my_fwrite(&LosGlobal->BoxSize, sizeof(double), 1, fd);
  my_fwrite(&LosGlobal->Wmax, sizeof(double), 1, fd);
  my_fwrite(&LosGlobal->Time, sizeof(double), 1, fd);
  my_fwrite(&LosGlobal->Xpos, sizeof(double), 1, fd);
  my_fwrite(&LosGlobal->Ypos, sizeof(double), 1, fd);
  my_fwrite(&LosGlobal->xaxis, sizeof(int), 1, fd);
  my_fwrite(&LosGlobal->yaxis, sizeof(int), 1, fd);
  my_fwrite(&LosGlobal->zaxis, sizeof(int), 1, fd);

  my_fwrite(LosGlobal->TauHI, sizeof(double), PIXELS, fd);
  my_fwrite(LosGlobal->TempHI, sizeof(double), PIXELS, fd);
  my_fwrite(LosGlobal->VpecHI, sizeof(double), PIXELS, fd);
  my_fwrite(LosGlobal->NHI, sizeof(double), PIXELS, fd);

  my_fwrite(LosGlobal->TauHeII, sizeof(double), PIXELS, fd);
  my_fwrite(LosGlobal->TempHeII, sizeof(double), PIXELS, fd);
  my_fwrite(LosGlobal->VpecHeII, sizeof(double), PIXELS, fd);
  my_fwrite(LosGlobal->NHeII, sizeof(double), PIXELS, fd);

  my_fwrite(LosGlobal->Rho, sizeof(double), PIXELS, fd);
  my_fwrite(LosGlobal->Vpec, sizeof(double), PIXELS, fd);
  my_fwrite(LosGlobal->Temp, sizeof(double), PIXELS, fd);
  my_fwrite(LosGlobal->Metallicity, sizeof(double), PIXELS, fd);

  fclose(fd);
}


/* int find_next_lineofsighttime(int time0) */
/* { */
/*   double a1, a2, df1, df2, u1, u2; */
/*   double logTimeBegin, logTimeMax; */
/*   int i1, i2, im, time1; */

/*   logTimeBegin = log(All.TimeBegin); */
/*   logTimeMax = log(All.TimeMax); */

/*   a1 = logTimeBegin + time0 * All.Timebase_interval; */

/*   u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH; */
/*   i1 = (int) u1; */
/*   if(i1 >= DRIFT_TABLE_LENGTH) */
/*     i1 = DRIFT_TABLE_LENGTH - 1; */

/*   if(i1 <= 1) */
/*     df1 = u1 * GravKickTable[0]; */
/*   else */
/*     df1 = GravKickTable[i1 - 1] + (GravKickTable[i1] - GravKickTable[i1 - 1]) * (u1 - i1); */

/*   df2 = df1 + All.BoxSize / (C / All.UnitVelocity_in_cm_per_s); */

/*   i2 = DRIFT_TABLE_LENGTH - 1; */

/*   while(i2 - i1 > 0) */
/*     { */
/*       im = (i1 - 1 + i2) / 2; */

/*       if(GravKickTable[im] > df2) */
/* 	i2 = im; */
/*       else */
/* 	i1 = im + 1; */
/*     } */

/*   u2 = (df2 - GravKickTable[i2 - 1]) / (GravKickTable[i2] - GravKickTable[i2 - 1]) + i2; */

/*   a2 = u2 * (logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH + logTimeBegin; */

/*   time1 = (int) ((a2 - logTimeBegin) / All.Timebase_interval + 0.5); */

/*   return time1; */
/* } */

double los_periodic(double x)
{
  if(x >= 0.5 * All.BoxSize)
    x -= All.BoxSize;
  if(x < -0.5 * All.BoxSize)
    x += All.BoxSize;

  return x;
}

#endif
#endif
