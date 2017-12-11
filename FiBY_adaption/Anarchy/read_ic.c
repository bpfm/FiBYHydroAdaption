#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#ifdef BG_SFR
#include "bg_proto.h"
#include "bg_vars.h"
#endif

/* This function reads initial conditions that are in the default file format
 * of Gadget, i.e. snapshot files can be used as input files.  However, when a
 * snapshot file is used as input, not all the information in the header is
 * used: THE STARTING TIME NEEDS TO BE SET IN THE PARAMETERFILE.
 * Alternatively, the code can be started with restartflag==2, then snapshots
 * from the code can be used as initial conditions-files without having to
 * change the parameterfile.  For gas particles, only the internal energy is
 * read, the density and mean molecular weight will be recomputed by the code.
 * When InitGasU_ERG > 0 is given, the gas entropy will be initialzed the
 * corresponding value.
 */


#ifdef BG_SFR
static int relative_star_offset;
#endif


#ifdef AUTO_SWAP_ENDIAN_READIC
int swap_file = 8;
#endif

#ifdef AUTO_SWAP_ENDIAN_READIC
/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to swap ENDIAN -------------------------------*/
/*-------- char *data:    Pointer to the data ---------------------------------*/
/*-------- int n:         Number of elements to swap --------------------------*/
/*-------- int m:         Size of single element to swap ----------------------*/
/*--------                int,float = 4 ---------------------------------------*/
/*--------                double    = 8 ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
void swap_Nbyte(char *data, int n, int m)
{
  int i, j;

  char old_data[16];

  if(swap_file != 8)
    {
      for(j = 0; j < n; j++)
	{
	  memcpy(&old_data[0], &data[j * m], m);
	  for(i = 0; i < m; i++)
	    {
	      data[j * m + i] = old_data[m - i - 1];
	    }
	}
    }
}

/*------------------------------------------------------------------*/
/*----------- procedure to swap header if needed -------------------*/
/*------------------------------------------------------------------*/

void swap_header()
{
  swap_Nbyte((char *) &header.npart, 6, 4);
  swap_Nbyte((char *) &header.mass, 6, 8);
  swap_Nbyte((char *) &header.time, 1, 8);
  swap_Nbyte((char *) &header.redshift, 1, 8);
  swap_Nbyte((char *) &header.flag_sfr, 1, 4);
  swap_Nbyte((char *) &header.flag_feedback, 1, 4);
  swap_Nbyte((char *) &header.npartTotal, 6, 4);
  swap_Nbyte((char *) &header.flag_cooling, 1, 4);
  swap_Nbyte((char *) &header.num_files, 1, 4);
  swap_Nbyte((char *) &header.BoxSize, 1, 8);
  swap_Nbyte((char *) &header.Omega0, 1, 8);
  swap_Nbyte((char *) &header.OmegaLambda, 1, 8);
  swap_Nbyte((char *) &header.HubbleParam, 1, 8);
  swap_Nbyte((char *) &header.flag_stellarage, 1, 4);
  swap_Nbyte((char *) &header.flag_metals, 1, 4);
  swap_Nbyte((char *) &header.npartTotalHighWord, 6, 4);
  swap_Nbyte((char *) &header.flag_entropy_instead_u, 1, 4);
}

#endif


void read_ic(char *fname)
{
  int i, num_files, rest_files, ngroups, gr, filenr, masterTask, lastTask, groupMaster;

  double dmax1, dmax2;

  char buf[500];


#ifdef BG_SFR
  int generations;

  double original_gas_mass, mass, masstot;

  generations = All.Generations;

#ifdef BG_STELLAR_EVOLUTION
  int k;
#endif
#endif

#ifdef RESCALEVINI
  if(ThisTask == 0)
    {
      fprintf(stdout, "\nRescaling v_ini !\n\n");
      fflush(stdout);
    }
#endif


  CPU_Step[CPU_MISC] += measure_time();


  NumPart = 0;
  N_gas = 0;
#if defined(BG_SFR)
  N_star = 0;
#endif
  All.TotNumPart = 0;

  num_files = find_files(fname);

  rest_files = num_files;

  while(rest_files > NTask)
    {
      sprintf(buf, "%s.%d", fname, ThisTask + (rest_files - NTask));
      if(All.ICFormat == 3)
	sprintf(buf, "%s.%d.hdf5", fname, ThisTask + (rest_files - NTask));

      ngroups = NTask / All.NumFilesWrittenInParallel;
      if((NTask % All.NumFilesWrittenInParallel))
	ngroups++;
      groupMaster = (ThisTask / ngroups) * ngroups;

      for(gr = 0; gr < ngroups; gr++)
	{
	  if(ThisTask == (groupMaster + gr))	/* ok, it's this processor's turn */
	    read_file(buf, ThisTask, ThisTask);
	  MPI_Barrier(MPI_COMM_WORLD);
	}

      rest_files -= NTask;
    }


  if(rest_files > 0)
    {
      distribute_file(rest_files, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

      if(num_files > 1)
	{
	  sprintf(buf, "%s.%d", fname, filenr);
	  if(All.ICFormat == 3)
	    sprintf(buf, "%s.%d.hdf5", fname, filenr);
	}
      else
	{
	  sprintf(buf, "%s", fname);
	  if(All.ICFormat == 3)
	    sprintf(buf, "%s.hdf5", fname);
	}

      ngroups = rest_files / All.NumFilesWrittenInParallel;
      if((rest_files % All.NumFilesWrittenInParallel))
	ngroups++;

      for(gr = 0; gr < ngroups; gr++)
	{
	  if((filenr / All.NumFilesWrittenInParallel) == gr)	/* ok, it's this processor's turn */
	    read_file(buf, masterTask, lastTask);
	  MPI_Barrier(MPI_COMM_WORLD);
	}
    }

  myfree(CommBuffer);

  /* this makes sure that masses are initialized in the case that the mass-block
     is completely empty */
  for(i = 0; i < NumPart; i++)
    {
      if(All.MassTable[P[i].Type] != 0)
	P[i].Mass = All.MassTable[P[i].Type];
    }


#ifdef GENERATE_GAS_IN_ICS
  int count, j;

  double fac, d, a, b, rho;

  if(RestartFlag == 0)
    {
      for(i = 0, count = 0; i < NumPart; i++)
	if(P[i].Type == 1)
	  count++;

      memmove(P + count, P, sizeof(struct particle_data) * NumPart);

      NumPart += count;
      N_gas += count;

      fac = All.OmegaBaryon / All.Omega0;
      rho = All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

      for(i = count, j = 0; i < NumPart; i++)
	if(P[i].Type == 1)
	  {
	    P[j] = P[i];

	    d = pow(P[i].Mass / rho, 1.0 / 3);
	    a = 0.5 * All.OmegaBaryon / All.Omega0 * d;
	    b = 0.5 * (All.Omega0 - All.OmegaBaryon) / All.Omega0 * d;

	    P[j].Mass *= fac;
	    P[i].Mass *= (1 - fac);
	    P[j].Type = 0;
	    P[j].ID += 1000000000;

	    P[i].Pos[0] += a;
	    P[i].Pos[1] += a;
	    P[i].Pos[2] += a;
	    P[j].Pos[0] -= b;
	    P[j].Pos[1] -= b;
	    P[j].Pos[2] -= b;

	    j++;
	  }

      All.MassTable[0] *= fac;
      All.MassTable[1] *= (1 - fac);
    }
#endif

#if defined(BLACK_HOLES) && defined(SWALLOWGAS)
  if(RestartFlag == 0)
    {
      All.MassTable[5] = 0;
    }
#endif

#ifdef BG_SFR
  if(RestartFlag == 0)
    {
      if(All.MassTable[4] == 0 && All.MassTable[0] > 0)
	{
	  All.OrigGasMass = All.MassTable[0];

	  All.MassTable[0] = 0;
	  All.MassTable[4] = 0;
	}
      else
	{
	  for(i = 0, mass = 0; i < NumPart; i++)
	    {
	      if(P[i].Type == 0 || P[i].Type == 4)
		mass += P[i].Mass;
	    }

	  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  original_gas_mass = masstot / (All.TotN_gas + header.npartTotal[4]);


	  if(header.npartTotal[4])
	    {
	      for(i = 0, mass = 0; i < NumPart; i++)
		{
		  if(P[i].Type == 4)
		    mass += P[i].Mass;
		}

	      MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	      original_gas_mass = generations * masstot / header.npartTotal[4];
	    }

	  All.OrigGasMass = original_gas_mass;
	}

      if(ThisTask == 0)
	printf("original_gas_mass= %g\n", All.OrigGasMass);
    }

  if(RestartFlag == 2)
    {
      for(i = 0, mass = 0; i < NumPart; i++)
	{
#ifdef BLACK_HOLES
	  if(P[i].Type == 0 || P[i].Type == 4 || P[i].Type == 5)
#else
	  if(P[i].Type == 0 || P[i].Type == 4)
#endif
	    mass += P[i].Mass;
	}

      MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      original_gas_mass = masstot / (All.TotN_gas + header.npartTotal[4]);

      if(ThisTask == 0)
	printf("original_gas_mass= %g\n", original_gas_mass);

      All.OrigGasMass = original_gas_mass;
    }
#endif



#ifdef BG_SFR
  All.InitGasU = All.InitGasU_ERG * All.UnitMass_in_g / All.UnitEnergy_in_cgs;	/* unit conversion */
#endif


  if(RestartFlag == 0)
    {
#ifdef BG_SFR
      if(All.InitGasU > 0)
#endif
	{
	  for(i = 0; i < N_gas; i++)
	    {
	      if(ThisTask == 0 && i == 0 && SphP[i].Entropy == 0)
#ifdef BG_SFR
		printf("Initializing u from InitGasU !\n");
#endif
	      if(SphP[i].Entropy == 0)
		SphP[i].Entropy = All.InitGasU;

	      /* Note: the coversion to entropy will be done in the function init(),
	         after the densities have been computed */
	    }
	}
    }


  for(i = 0; i < N_gas; i++)
    SphP[i].Entropy = DMAX(All.MinEgySpec, SphP[i].Entropy);


/* convert metal abundances to metal masses */
#ifdef BG_STELLAR_EVOLUTION
  if(RestartFlag == 2)
    {
      if(ThisTask == 0)
	printf(" Converting element abundances to metal masses ... \n");

      for(i = 0; i < NumPart; i++)
	{
	  if(P[i].Type == 0)	/* gas particle */
	    {
	      for(k = 0; k < BG_NELEMENTS; k++)
		SphP[i].Metals[k] *= P[i].Mass;
#ifdef BG_SNIA_IRON
	      SphP[i].IronFromSNIa *= P[i].Mass;
#endif

#ifdef BG_METALSMOOTHING
	      for(k = 0; k < BG_NELEMENTS; k++)
		SphP[i].MetalsSmoothed[k] *= P[i].Mass;
#ifdef BG_SNIA_IRON
	      SphP[i].IronFromSNIaSmoothed *= P[i].Mass;
#endif
#endif
	    }

	  if(P[i].Type == 4)	/* star particle */
	    {
	      for(k = 0; k < BG_NELEMENTS; k++)
		StarP[P[i].StarID].Metals[k] *= P[i].Mass;
#ifdef BG_SNIA_IRON
	      StarP[P[i].StarID].IronFromSNIa *= StarP[P[i].StarID].InitialMass;
#endif
#ifdef BG_METALSMOOTHING
	      for(k = 0; k < BG_NELEMENTS; k++)
		StarP[P[i].StarID].MetalsSmoothed[k] *= StarP[P[i].StarID].InitialMass;
#ifdef BG_SNIA_IRON
	      StarP[P[i].StarID].IronFromSNIaSmoothed *= StarP[P[i].StarID].InitialMass;
#endif
#endif
	    }
	}
    }
#endif /* BG_STELLAR_EVOLUTION */


  MPI_Barrier(MPI_COMM_WORLD);


  if(ThisTask == 0)
    {
      printf("reading done.\n");
      fflush(stdout);
    }

  if(ThisTask == 0)
    {
      printf("Total number of particles :  %d%09d\n\n",
	     (int) (All.TotNumPart / 1000000000), (int) (All.TotNumPart % 1000000000));
      fflush(stdout);
    }

#ifdef USE_HDF5_FIX
  if(All.ICFormat == 3)
    {
      H5close();
      hdf5_memory_cleanup();
    }
#endif


  CPU_Step[CPU_SNAPSHOT] += measure_time();
}


/*! This function reads out the buffer that was filled with particle data.
 */
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type)
{
  int n, k;
  float *fp;
  MyIDType *ip;

#ifdef AUTO_SWAP_ENDIAN_READIC
  int vt, vpb;

  char *cp;
#endif

  fp = (float *) CommBuffer;
  ip = (MyIDType *) CommBuffer;

#ifdef AUTO_SWAP_ENDIAN_READIC
  cp = (char *) CommBuffer;
  vt = get_datatype_in_block(blocknr);
  vpb = get_values_per_blockelement(blocknr);
  if(vt == 2)
    swap_Nbyte(cp, pc * vpb, 8);
  else
    swap_Nbyte(cp, pc * vpb, 4);
#endif

  switch (blocknr)
    {
    case IO_POS:		/* positions */
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
	  P[offset + n].Pos[k] = *fp++;

      for(n = 0; n < pc; n++)
	P[offset + n].Type = type;	/* initialize type here as well */
      break;

    case IO_VEL:		/* velocities */
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
#ifdef RESCALEVINI
	  /* scaling v to use same IC's for different cosmologies */
	  P[offset + n].Vel[k] = (*fp++) * All.VelIniScale;
#else
	  P[offset + n].Vel[k] = *fp++;
#endif
      break;

    case IO_ID:		/* particle ID */
      for(n = 0; n < pc; n++)
	P[offset + n].ID = *ip++;
      break;

    case IO_MASS:		/* particle mass */
      for(n = 0; n < pc; n++)
	P[offset + n].Mass = *fp++;
      break;

    case IO_AO:		/* particle mass */
#ifdef ADAPTIVE_OUTPUT
      for(n = 0; n < pc; n++) {
	P[offset + n].StepsSinceLastOutput = *fp++;
      }
#endif
      break;


    case IO_U:			/* temperature */
      for(n = 0; n < pc; n++)
	SphP[offset + n].Entropy = *fp++;
      break;

    case IO_RHO:		/* density */
#ifndef BG_SFR
      for(n = 0; n < pc; n++)
	SphP[offset + n].d.Density = *fp++;
#else
      if(type == 4)
	for(n = 0; n < pc; n++)
	  StarP[offset + n - relative_star_offset].GasDensity = *fp++;
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset + n].d.Density = *fp++;
#endif
      break;

    case IO_WRHO:		/* weighted density */
#ifdef PRESSURE_ENTROPY_SPH
      for(n = 0; n < pc; n++)
	SphP[offset + n].cky.WeightedDensity = *fp++;
#endif
      break;

    case IO_HSML:		/* SPH smoothing length */
      for(n = 0; n < pc; n++)
	PPP[offset + n].Hsml = *fp++;
      break;

    case IO_AGE:		/* Age of stars */
      for(n = 0; n < pc; n++)
	{
#ifdef BG_SFR
	  P[offset + n].StarID = offset + n - relative_star_offset;
	  StarP[offset + n - relative_star_offset].PID = offset + n;
	  StarP[offset + n - relative_star_offset].StarBirthTime = *fp++;
#endif
	}
      break;


    case IO_Z:			/* Gas and star metallicity */
      break;

    case IO_BFLD:		/* Magnetic field */
      break;

    case IO_EGYPROM:
      break;

    case IO_EGYCOLD:
      break;

    case IO_CR_C0:		/* Adiabatic invariant for cosmic rays */
      break;

    case IO_CR_Q0:		/* Adiabatic invariant for cosmic rays */
      break;

    case IO_CR_P0:
      break;

    case IO_CR_E0:
      break;

    case IO_CR_n0:
      break;

    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
      break;

    case IO_BHMASS:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; n++)
	P[offset + n].BH_Mass = *fp++;
#endif
      break;

    case IO_BHMDOT:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; n++)
	P[offset + n].BH_Mdot = *fp++;
#endif
      break;

    case IO_BH_ENERGY:
#if defined(BLACK_HOLES) && defined(BH_THERMALFEEDBACK)
      for(n = 0; n < pc; n++)
	P[offset + n].BH_Energy = *fp++;
#endif
      break;

    case IO_BH_BIRTH_TIME:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; n++)
	P[offset + n].BH_BirthTime = *fp++;
#endif
      break;

    case IO_Zs:
      break;

    case IO_iMass:
      break;

    case IO_BG_SNII_KINETIC_FEEDBACK:
#ifdef BG_SNII_KINETIC_FEEDBACK
      if(type == 4)
	for(n = 0; n < pc; n++)
	  StarP[offset + n - relative_star_offset].WindFlag = *fp++;
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset + n].WindFlag = *fp++;
#endif
      break;

    case IO_BG_METALS:
#ifdef BG_STELLAR_EVOLUTION
      if(type == 4)
	for(n = 0; n < pc; n++)
	  for(k = 0; k < BG_NELEMENTS; k++)
	    StarP[offset + n - relative_star_offset].Metals[k] = *fp++;
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  for(k = 0; k < BG_NELEMENTS; k++)
	    SphP[offset + n].Metals[k] = *fp++;
#endif
      break;

    case IO_BG_METALS_SMOOTHED:
#ifdef BG_METALSMOOTHING
      if(type == 4)
	for(n = 0; n < pc; n++)
	  for(k = 0; k < BG_NELEMENTS; k++)
	    StarP[offset + n - relative_star_offset].MetalsSmoothed[k] = *fp++;
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  for(k = 0; k < BG_NELEMENTS; k++)
	    SphP[offset + n].MetalsSmoothed[k] = *fp++;
#endif
      break;

    case IO_BG_INITIAL_MASS:
#ifdef BG_STELLAR_EVOLUTION
      for(n = 0; n < pc; n++)
	StarP[offset + n - relative_star_offset].InitialMass = *fp++;
#endif
      break;

    case IO_BG_METALLICITY:
#ifdef BG_STELLAR_EVOLUTION
      if(type == 4)
	for(n = 0; n < pc; n++)
	  StarP[offset + n - relative_star_offset].Metallicity = *fp++;
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset + n].Metallicity = *fp++;
#endif
      break;

    case IO_BG_METALLICITY_SMOOTHED:
#ifdef BG_METALSMOOTHING
      if(type == 4)
	for(n = 0; n < pc; n++)
	  StarP[offset + n - relative_star_offset].MetallicitySmoothed = *fp++;
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset + n].MetallicitySmoothed = *fp++;
#endif
      break;

    case IO_BG_IRON_FROM_SNIA:
#ifdef BG_SNIA_IRON
      if(type == 4)
	for(n = 0; n < pc; n++)
	  StarP[offset + n - relative_star_offset].IronFromSNIa = *fp++;
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset + n].IronFromSNIa = *fp++;
#endif
      break;

    case IO_BG_IRON_FROM_SNIA_SMOOTHED:
#if defined(BG_SNIA_IRON) && defined(BG_METALSMOOTHING)
      if(type == 4)
	for(n = 0; n < pc; n++)
	  StarP[offset + n - relative_star_offset].IronFromSNIaSmoothed = *fp++;
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset + n].IronFromSNIaSmoothed = *fp++;
#endif
      break;

    case IO_BG_MAX_ENTROPY:
#ifdef BG_EXTRA_ARRAYS
      if(type == 4)
	for(n = 0; n < pc; n++)
	  StarP[offset + n - relative_star_offset].MaximumEntropy = *fp++;
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset + n].MaximumEntropy = *fp++;
#endif
      break;

    case IO_BG_MAX_TEMPERATURE:
#ifdef BG_EXTRA_ARRAYS
      if(type == 4)
	for(n = 0; n < pc; n++)
	  StarP[offset + n - relative_star_offset].MaximumTemperature = *fp++;
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset + n].MaximumTemperature = *fp++;
#endif
      break;

    case IO_BG_TIME_MAX_ENTROPY:
#ifdef BG_EXTRA_ARRAYS
      if(type == 4)
	for(n = 0; n < pc; n++)
	  StarP[offset + n - relative_star_offset].TimeMaximumEntropy = *fp++;
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset + n].TimeMaximumEntropy = *fp++;
#endif
      break;

    case IO_BG_TIME_MAX_TEMPERATURE:
#ifdef BG_EXTRA_ARRAYS
      if(type == 4)
	for(n = 0; n < pc; n++)
	  StarP[offset + n - relative_star_offset].TimeMaximumTemperature = *fp++;
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset + n].TimeMaximumTemperature = *fp++;
#endif
      break;

    case IO_BG_METALLICITY_WEIGHTED_REDSHIFT:
#if defined(BG_STELLAR_EVOLUTION) && defined(BG_Z_WEIGHTED_REDSHIFT)
      if(type == 4)
	{
	  for(n = 0; n < pc; n++)
	    StarP[offset + n - relative_star_offset].MetallicityWeightedRedshift = *fp++;
	}
      else if(type == 0)
	{
	  for(n = 0; n < pc; n++)
	    SphP[offset + n].MetallicityWeightedRedshift = *fp++;
	}
#endif
      break;

    case IO_BG_METALLICITY_WEIGHTED_POTENTIAL:
      break;

    case IO_SFR:
#ifdef BG_SFR
      for(n = 0; n < pc; n++)
	SphP[offset + n].Sfr = *fp++;
#endif
      break;

    case IO_BG_ONEOS:
#ifdef BG_SFR
      for(n = 0; n < pc; n++)
	SphP[offset + n].OnEOS = *fp++;
#endif
      break;

    case  IO_BG_e:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; n++)
	SphP[offset + n].n_e = *fp++;
#endif
	break;

    case  IO_BG_HII:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; n++)
	SphP[offset + n].x_Hp = *fp++;
#endif
	break;

    case  IO_BG_Hminus:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; n++)
	SphP[offset + n].x_Hm = *fp++;
#endif
	break;

    case  IO_BG_H2I:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; n++)
	SphP[offset + n].x_H2 = *fp++;
#endif
	break;

    case  IO_BG_H2II:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; n++)
	SphP[offset + n].x_H2p = *fp++;
#endif
	break;

    case  IO_BG_HeII:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; n++)
	SphP[offset + n].x_Hep = *fp++;
#endif
	break;

    case  IO_BG_HeIII:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; n++)
	SphP[offset + n].x_Hepp = *fp++;
#endif
	break;

    case  IO_BG_DII:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; n++)
	SphP[offset + n].x_Dp = *fp++;
#endif
	break;

    case  IO_BG_HD:
#ifdef BG_MOL_NETWORK
      for(n = 0; n < pc; n++)
	SphP[offset + n].x_HD = *fp++;
#endif
	break;

/*     case  IO_BG_DUST: */
/* #ifdef BG_DUST */
/*       for(n = 0; n < pc; n++) */
/*         SphP[offset + n].Dusticity = *fp++; */
/* #endif */
/*       break; */

/*     case  IO_BG_DUST_SMOOTHED: */
/* #ifdef BG_DUST */
/* #ifdef BG_METALSMOOTHING */
/*       for(n = 0; n < pc; n++) */
/*         SphP[offset + n].DusticitySmoothed = *fp++; */
/* #endif */
/* #endif */
/*       break; */

/*     case  IO_BG_DUST_SIZE: */
/* #ifdef BG_DUST_DESTRUCTION_SUBLIMATION */
/*       for(n = 0; n < pc; n++) */
/*         SphP[offset + n].DustGrainSize = *fp++; */
/* #endif */
/*       break; */

      /* the other input fields (if present) are not needed to define the 
         initial conditions of the code */

    case IO_POT:
    case IO_ACCEL:
    case IO_DTENTR:
    case IO_TSTP:
    case IO_DBDT:
    case IO_DIVB:
    case IO_ABVC:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_BSMTH:
    case IO_DENN:
    case IO_MACH:
    case IO_DTENERGY:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
    case IO_AMDC:
    case IO_PHI:
    case IO_CLDX:
      break;

#ifdef SUBFIND
    case IO_SUBFIND_MASS:
    case IO_SUBFIND_MASSTYPE:
    case IO_SUBFIND_POS:
    case IO_SUBFIND_CMPOS:
    case IO_SUBFIND_HALFMASS:
    case IO_SUBFIND_HALFMASSPROJ:
    case IO_SUBFIND_VMAXRAD:
    case IO_SUBFIND_SPIN:
    case IO_SUBFIND_VMAX:
    case IO_SUBFIND_CMVEL:
    case IO_SUBFIND_VELDISP:
    case IO_SUBFIND_STELLARVELDISP:
    case IO_SUBFIND_STELLARVELDISPHALFPROJ:
    case IO_SUBFIND_PARTPOS:
    case IO_SUBFIND_PARTVEL:
    case IO_SUBFIND_PARTMASS:
    case IO_M_MEAN:
    case IO_M_CRIT:
    case IO_M_TOPHAT:
    case IO_R_MEAN:
    case IO_R_CRIT:
    case IO_R_TOPHAT:
      break;
#endif

/*     case IO_LW_popII: */
/*     case IO_LW_popIII: */
    case IO_BG_STELLAR_AGE:
    case IO_BG_TEMP:
      break;

    case IO_LASTENTRY:
      endrun(220);
      break;
    }
}



/*! This function reads a snapshot file and distributes the data it contains
 *  to tasks 'readTask' to 'lastTask'.
 */
void read_file(char *fname, int readTask, int lastTask)
{
  int blockmaxlen;

  int i, n_in_file, n_for_this_task, ntask, pc, offset = 0, task;

  int blksize1, blksize2;

  MPI_Status status;

  FILE *fd = 0;

  int nall;

  int type, bnr;

  char label[4], expected_label[4];

  int nstart, bytes_per_blockelement, npart, nextblock, typelist[6];

  enum iofields blocknr;

  size_t bytes;

  double t0, t1;

#ifdef HAVE_HDF5
  char buf[500];

  int pcsum;

  hid_t hdf5_file = 0, hdf5_grp[6];

  hid_t this_grp = 0;

  hsize_t count[2], start[2], dims[2];

  hid_t hdf5_dataspace_in_file, hdf5_dataset;

  hid_t hdf5_datatype = 0, hdf5_dataspace_in_memory = 0;

  int rank;

#ifdef BG_STELLAR_EVOLUTION
  int k;

  char gbuf[500], mbuf[500];

  hid_t element_grp[6], selement_grp[6];

  hid_t element_dset[BG_NELEMENTS], hdf5_element_dataspace_in_file[BG_NELEMENTS];

  hsize_t myoffset[2];

  herr_t hdf5_status;
#endif
#endif

#define SKIP  {my_fread(&blksize1,sizeof(int),1,fd);}
#define SKIP2  {my_fread(&blksize2,sizeof(int),1,fd);}


  if(ThisTask == readTask)
    {
      if(All.ICFormat == 1 || All.ICFormat == 2)
	{
	  if(!(fd = fopen(fname, "r")))
	    {
	      printf("can't open file `%s' for reading initial conditions.\n", fname);
	      endrun(123);
	    }

	  if(All.ICFormat == 2)
	    {
	      SKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_file = blksize1;
#endif
	      my_fread(&label, sizeof(char), 4, fd);
	      my_fread(&nextblock, sizeof(int), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_Nbyte((char *) &nextblock, 1, 4);
#endif
	      printf("Reading header => '%c%c%c%c' (%d byte)\n", label[0], label[1], label[2], label[3],
		     nextblock);
	      SKIP2;
	    }

	  SKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
	  if(All.ICFormat == 1)
	    {
	      if(blksize1 != 256)
		swap_file = 1;
	    }
#endif
	  my_fread(&header, sizeof(header), 1, fd);
	  SKIP2;

	  /* NEW */
	  /*
	  if(All.ICFormat == 1)
	    {
	      header.npartTotalHighWord[0] = 1;
	      header.npartTotalHighWord[1] = 1;
	    }
	  */
	  /* NEW */

#ifdef AUTO_SWAP_ENDIAN_READIC
	  swap_Nbyte((char *) &blksize1, 1, 4);
	  swap_Nbyte((char *) &blksize2, 1, 4);
#endif

	  if(blksize1 != 256 || blksize2 != 256)
	    {
	      printf("incorrect header format\n");
	      fflush(stdout);
	      endrun(890);
	    }
#ifdef AUTO_SWAP_ENDIAN_READIC
	  swap_header();
#endif
	}


#ifdef HAVE_HDF5
      if(All.ICFormat == 3)
	{
	  for(type = 0; type < 6; type++)
	    {
	      hdf5_grp[type] = -1;
#ifdef BG_STELLAR_EVOLUTION
	      element_grp[type] = -1;
	      selement_grp[type] = -1;
#endif
	    }

	  read_header_attributes_in_hdf5(fname);

	  hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

	  for(type = 0; type < 6; type++)
	    {
	      if(header.npart[type] > 0)
		{
		  sprintf(buf, "/PartType%d", type);
		  hdf5_grp[type] = H5Gopen(hdf5_file, buf);
#ifdef BG_STELLAR_EVOLUTION
		  if(type == 0 || type == 4)
		    element_grp[type] = H5Gopen(hdf5_grp[type], "ElementAbundance");
		  /* group for smoothed metallicities */
		  if(type == 0 || type == 4)
		    selement_grp[type] = H5Gopen(hdf5_grp[type], "SmoothedElementAbundance");
#endif
		}
	    }
	}
#endif

      for(task = readTask + 1; task <= lastTask; task++)
	MPI_Ssend(&header, sizeof(header), MPI_BYTE, task, TAG_HEADER, MPI_COMM_WORLD);
    }
  else
    MPI_Recv(&header, sizeof(header), MPI_BYTE, readTask, TAG_HEADER, MPI_COMM_WORLD, &status);

#ifdef AUTO_SWAP_ENDIAN_READIC
  MPI_Bcast(&swap_file, sizeof(swap_file), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  if(All.TotNumPart == 0)
    {
      if(header.num_files <= 1)
	for(i = 0; i < 6; i++)
	  {
	    header.npartTotal[i] = header.npart[i];
#ifdef BG_SFR
	    header.npartTotalHighWord[i] = 0;
#endif
	  }

      All.TotN_gas = header.npartTotal[0] + (((long long) header.npartTotalHighWord[0]) << 32);
#if defined(BG_SFR)
      All.TotN_star = header.npartTotal[4];
#endif

      for(i = 0, All.TotNumPart = 0; i < 6; i++)
	{
	  All.TotNumPart += header.npartTotal[i];
	  All.TotNumPart += (((long long) header.npartTotalHighWord[i]) << 32);
	}

#ifdef GENERATE_GAS_IN_ICS
      if(RestartFlag == 0)
	{
	  All.TotN_gas += header.npartTotal[1];
	  All.TotNumPart += header.npartTotal[1];
	}
#endif

      for(i = 0; i < 6; i++)
	All.MassTable[i] = header.mass[i];

      All.MaxPart = (int) (All.PartAllocFactor * (All.TotNumPart / NTask));	/* sets the maximum number of particles that may */
      All.MaxPartSph = (int) (All.PartAllocFactor * (All.TotN_gas / NTask));	/* sets the maximum number of particles that may 
										   reside on a processor */
#ifdef INHOMOG_GASDISTR_HINT
      All.MaxPartSph = All.MaxPart;
#endif

#ifdef BG_SFR
      if(All.TotN_star == 0)
	All.MaxPartStar =
	  (int) (All.PartAllocFactor * (All.TotN_gas / NTask) * All.Generations * All.StarAllocFactor);
      else
	All.MaxPartStar =
	  (int) (All.PartAllocFactor * (All.TotN_star / NTask + All.TotN_gas / NTask) * All.Generations *
		 All.StarAllocFactor);
#endif

      allocate_memory();

      if(!(CommBuffer = mymalloc(bytes = All.BufferSize * 1024 * 1024)))
	{
	  printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(2);
	}

      if(RestartFlag == 2)
	All.Time = All.TimeBegin = header.time;
    }

  if(ThisTask == readTask)
    {
      for(i = 0, n_in_file = 0; i < 6; i++)
	n_in_file += header.npart[i];

      printf("\nreading file `%s' on task=%d (contains %d particles.)\n"
	     "distributing this file to tasks %d-%d\n"
	     "Type 0 (gas):   %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 1 (halo):  %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 2 (disk):  %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 3 (bulge): %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 4 (stars): %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 5 (bndry): %8d  (tot=%6d%09d) masstab=%g\n\n", fname, ThisTask, n_in_file, readTask,
	     lastTask, header.npart[0], (int) (header.npartTotal[0] / 1000000000),
	     (int) (header.npartTotal[0] % 1000000000), All.MassTable[0], header.npart[1],
	     (int) (header.npartTotal[1] / 1000000000), (int) (header.npartTotal[1] % 1000000000),
	     All.MassTable[1], header.npart[2], (int) (header.npartTotal[2] / 1000000000),
	     (int) (header.npartTotal[2] % 1000000000), All.MassTable[2], header.npart[3],
	     (int) (header.npartTotal[3] / 1000000000), (int) (header.npartTotal[3] % 1000000000),
	     All.MassTable[3], header.npart[4], (int) (header.npartTotal[4] / 1000000000),
	     (int) (header.npartTotal[4] % 1000000000), All.MassTable[4], header.npart[5],
	     (int) (header.npartTotal[5] / 1000000000), (int) (header.npartTotal[5] % 1000000000),
	     All.MassTable[5]);
      fflush(stdout);
    }


  ntask = lastTask - readTask + 1;


  /* to collect the gas particles all at the beginning (in case several
     snapshot files are read on the current CPU) we move the collisionless
     particles such that a gap of the right size is created */

#ifdef BG_SFR
  relative_star_offset = N_gas - N_star;
#endif

  for(type = 0, nall = 0; type < 6; type++)
    {
      n_in_file = header.npart[type];

      n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
	n_for_this_task++;


      if(type == 0)
	{
	  if(N_gas + n_for_this_task > All.MaxPartSph)
	    {
	      printf("Not enough space on task=%d for SPH particles (space for %d, need at least %d)\n",
		     ThisTask, All.MaxPartSph, N_gas + n_for_this_task);
	      fflush(stdout);
	      endrun(172);
	    }
	}

#ifdef BG_SFR
      if(type < 4)
	relative_star_offset += n_for_this_task;
#endif

      nall += n_for_this_task;
    }

  if(NumPart + nall > All.MaxPart)
    {
      printf("Not enough space on task=%d (space for %d, need at least %d)\n",
	     ThisTask, All.MaxPart, NumPart + nall);
      fflush(stdout);
      endrun(173);
    }

  memmove(&P[N_gas + nall], &P[N_gas], (NumPart - N_gas) * sizeof(struct particle_data));
  nstart = N_gas;

#ifdef BG_SFR
  for(i = 0; i < N_star; i++)
    StarP[i].PID += nall;
#endif


  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum iofields) bnr;

      if(blocknr == IO_LASTENTRY)
	break;

      if(blockpresent(blocknr))
	{
#ifndef READ_HSML
	  if(RestartFlag == 0 && blocknr > IO_U && blocknr != IO_BFLD)
#else
	  if(RestartFlag == 0 && blocknr > IO_U && blocknr != IO_BFLD && blocknr != IO_HSML)
#endif
	    continue;		/* ignore all other blocks in initial conditions */


#ifdef BINISET
	  if(RestartFlag == 0 && blocknr == IO_BFLD)
	    continue;
#endif
	  t0 = second();

	  if(ThisTask == readTask)
	    {
	      get_dataset_name(blocknr, buf);
	      printf("reading block %d (%s)... ", blocknr, buf);
	      fflush(stdout);
	    }

	  bytes_per_blockelement = get_bytes_per_blockelement(blocknr);

	  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / bytes_per_blockelement;

	  npart = get_particles_in_block(blocknr, &typelist[0]);

	  if(npart > 0)
	    {
	      if(ThisTask == readTask)
		{
		  if(All.ICFormat == 2)
		    {
		      SKIP;
		      my_fread(&label, sizeof(char), 4, fd);
		      my_fread(&nextblock, sizeof(int), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
		      swap_Nbyte((char *) &nextblock, 1, 4);
#endif
		      printf("Reading header => '%c%c%c%c' (%d byte)\n", label[0], label[1], label[2],
			     label[3], nextblock);
		      SKIP2;

		      get_Tab_IO_Label(blocknr, expected_label);
		      if(strncmp(label, expected_label, 4) != 0)
			{
			  printf("incorrect block-structure!\n");
			  printf("expected '%c%c%c%c' but found '%c%c%c%c'\n",
				 label[0], label[1], label[2], label[3],
				 expected_label[0], expected_label[1], expected_label[2], expected_label[3]);
			  fflush(stdout);
			  endrun(1890);
			}
		    }
		  if(All.ICFormat == 1 || All.ICFormat == 2)
		    SKIP;
		}

	      for(type = 0, offset = 0; type < 6; type++)
		{
		  n_in_file = header.npart[type];
#ifdef HAVE_HDF5
		  pcsum = 0;
#endif
		  if(typelist[type] == 0)
		    {
		      n_for_this_task = n_in_file / ntask;
		      if((ThisTask - readTask) < (n_in_file % ntask))
			n_for_this_task++;

		      offset += n_for_this_task;
		    }
		  else
		    {
		      for(task = readTask; task <= lastTask; task++)
			{
			  n_for_this_task = n_in_file / ntask;
			  if((task - readTask) < (n_in_file % ntask))
			    n_for_this_task++;

			  if(task == ThisTask)
			    if(NumPart + n_for_this_task > All.MaxPart)
			      {
				printf("too many particles. %d %d %d\n", NumPart, n_for_this_task,
				       All.MaxPart);
				endrun(1313);
			      }


			  do
			    {
			      pc = n_for_this_task;

			      if(pc > blockmaxlen)
				pc = blockmaxlen;

			      if(ThisTask == readTask)
				{
				  if(All.ICFormat == 1 || All.ICFormat == 2)
				    my_fread(CommBuffer, bytes_per_blockelement, pc, fd);
#ifdef HAVE_HDF5
				  if(All.ICFormat == 3 && pc > 0)
				    {
				      this_grp = hdf5_grp[type];
#ifdef BG_STELLAR_EVOLUTION
				      /* group name */
				      get_group_name(blocknr, gbuf);
				      if(strcmp(gbuf, "ElementAbundance") == 0)
					this_grp = element_grp[type];
				      if(strcmp(gbuf, "SmoothedElementAbundance") == 0)
					this_grp = selement_grp[type];
#endif
				      /* dataset name */
				      get_dataset_name(blocknr, buf);

				      /* data type */
				      switch (get_datatype_in_block(blocknr))
					{
					case 0:
					  hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
					  break;
					case 1:
					  hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
					  break;
					case 2:
					  hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
					  break;
					}

#ifdef BG_STELLAR_EVOLUTION
				      /* special cases first */
				      if(strcmp(gbuf, "ElementAbundance") == 0
					 || strcmp(gbuf, "SmoothedElementAbundance") == 0)
					{
					  rank = 1;
					  dims[0] = header.npart[type];
					  dims[1] = 1;
					  for(k = 0; k < BG_NELEMENTS; k++)
					    {
					      strcpy(mbuf, ElementNames[k]);
					      hdf5_element_dataspace_in_file[k] =
						H5Screate_simple(rank, dims, NULL);
					      element_dset[k] = H5Dopen(this_grp, mbuf);
					    }
					  /* create data space in memory */
					  rank = 2;
					  dims[0] = pc;
					  dims[1] = BG_NELEMENTS;
					  hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);

					  /* read each metal */
					  for(k = 0; k < BG_NELEMENTS; k++)
					    {
					      /* select hyperslab in memory */
					      myoffset[0] = 0;
					      myoffset[1] = k;
					      count[0] = pc;
					      count[1] = 1;
					      hdf5_status =
						H5Sselect_hyperslab(hdf5_dataspace_in_memory, H5S_SELECT_SET,
								    myoffset, NULL, count, NULL);

					      /* select hyperslab in file */
					      myoffset[0] = pcsum;
					      count[0] = pc;
					      hdf5_status =
						H5Sselect_hyperslab(hdf5_element_dataspace_in_file[k],
								    H5S_SELECT_SET, myoffset, NULL, count,
								    NULL);

					      /* read */
					      hdf5_status =
						H5Dread(element_dset[k], hdf5_datatype,
							hdf5_dataspace_in_memory,
							hdf5_element_dataspace_in_file[k], H5P_DEFAULT,
							CommBuffer);

					    }

					  /* */
					  pcsum += pc;
					  H5Sclose(hdf5_dataspace_in_memory);

					  for(k = 0; k < BG_NELEMENTS; k++)
					    {
					      H5Sclose(hdf5_element_dataspace_in_file[k]);
					      H5Dclose(element_dset[k]);
					    }
					}
				      else
					{
#endif
					  /* array size in file */
					  dims[0] = header.npart[type];
					  dims[1] = get_values_per_blockelement(blocknr);
					  if(dims[1] == 1)
					    rank = 1;
					  else
					    rank = 2;

					  hdf5_dataset = H5Dopen(this_grp, buf);
					  hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);


					  /* data space in memory */
					  dims[0] = pc;
					  hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);

					  start[0] = pcsum;
					  start[1] = 0;

					  count[0] = pc;
					  count[1] = get_values_per_blockelement(blocknr);
					  pcsum += pc;

					  H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET,
							      start, NULL, count, NULL);

					  H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory,
						  hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);

					  H5Sclose(hdf5_dataspace_in_memory);
					  H5Sclose(hdf5_dataspace_in_file);
					  H5Dclose(hdf5_dataset);
#ifdef BG_STELLAR_EVOLUTION
					}
#endif
				      H5Tclose(hdf5_datatype);
				    }
#endif
				}

			      if(ThisTask == readTask && task != readTask && pc > 0)
				MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA,
					  MPI_COMM_WORLD);

			      if(ThisTask != readTask && task == ThisTask && pc > 0)
				MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, readTask,
					 TAG_PDATA, MPI_COMM_WORLD, &status);

			      if(ThisTask == task)
				{
				  empty_read_buffer(blocknr, nstart + offset, pc, type);

				  offset += pc;
				}

			      n_for_this_task -= pc;
			    }
			  while(n_for_this_task > 0);
			}
		    }
		}


	      if(ThisTask == readTask)
		{
		  if(All.ICFormat == 1 || All.ICFormat == 2)
		    {
		      SKIP2;
#ifdef AUTO_SWAP_ENDIAN_READIC
		      swap_Nbyte((char *) &blksize1, 1, 4);
		      swap_Nbyte((char *) &blksize2, 1, 4);
#endif
		      if(blksize1 != blksize2)
			{
			  printf("incorrect block-sizes detected!\n");
			  printf("Task=%d   blocknr=%d  blksize1=%d  blksize2=%d\n", ThisTask,
				 blocknr, blksize1, blksize2);
			  if(blocknr == IO_ID)
			    {
			      printf
				("Possible mismatch of 32bit and 64bit ID's in IC file and GADGET compilation !\n");
			    }
			  fflush(stdout);
			  endrun(1889);
			}
		    }
		}
	    }


	  t1 = second();

	  if(ThisTask == readTask)
	    {
	      printf("(took %g sec)\n", timediff(t0, t1));
	      fflush(stdout);
	    }
	}
    }


  for(type = 0; type < 6; type++)
    {
      n_in_file = header.npart[type];

      n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
	n_for_this_task++;

      NumPart += n_for_this_task;

      if(type == 0)
	N_gas += n_for_this_task;
#if defined(BG_SFR)
      if(type == 4)
	N_star += n_for_this_task;
#endif
    }

  if(ThisTask == readTask)
    {
      if(All.ICFormat == 1 || All.ICFormat == 2)
	fclose(fd);
#ifdef HAVE_HDF5
      if(All.ICFormat == 3)
	{
	  for(type = 5; type >= 0; type--)
	    if(header.npart[type] > 0)
	      {
		H5Gclose(hdf5_grp[type]);
#ifdef BG_STELLAR_EVOLUTION
		if(element_grp[type] > 0)
		  H5Gclose(element_grp[type]);

		if(selement_grp[type] > 0)
		  H5Gclose(selement_grp[type]);
#endif
	      }
	  H5Fclose(hdf5_file);
	}
#endif
    }

  if(ThisTask == 0)
    printf("Done variables input\n\n");
}



/*! This function determines on how many files a given snapshot is distributed.
 */
int find_files(char *fname)
{
  FILE *fd;

  char buf[200], buf1[200];

  int dummy, fnumber = 0;

  sprintf(buf, "%s.%d", fname, fnumber);
  sprintf(buf1, "%s", fname);

  if(All.ICFormat == 3)
    {
      sprintf(buf, "%s.%d.hdf5", fname, fnumber);
      sprintf(buf1, "%s.hdf5", fname);
    }

#ifndef  HAVE_HDF5
  if(All.ICFormat == 3)
    {
      if(ThisTask == 0)
	printf("Code wasn't compiled with HDF5 support enabled!\n");
      endrun(0);
    }
#endif

  header.num_files = 0;

  if(ThisTask == 0)
    {
      if((fd = fopen(buf, "r")))
	{
	  if(All.ICFormat == 1 || All.ICFormat == 2)
	    {
	      if(All.ICFormat == 2)
		{
		  my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
		  swap_file = dummy;
#endif
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		}

	      my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      if(All.SnapFormat == 1)
		{
		  if(dummy == 256)
		    swap_file = 8;
		  else
		    swap_file = dummy;
		}
#endif
	      my_fread(&header, sizeof(header), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_header();
#endif
	      my_fread(&dummy, sizeof(dummy), 1, fd);
	    }
	  fclose(fd);

#ifdef HAVE_HDF5
	  if(All.ICFormat == 3)
	    read_header_attributes_in_hdf5(buf);
#endif
	}
    }

#ifdef AUTO_SWAP_ENDIAN_READIC
  MPI_Bcast(&swap_file, sizeof(swap_file), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
  MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(header.num_files > 0)
    return header.num_files;

  if(ThisTask == 0)
    {
      if((fd = fopen(buf1, "r")))
	{
	  if(All.ICFormat == 1 || All.ICFormat == 2)
	    {
	      if(All.ICFormat == 2)
		{
		  my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
		  swap_file = dummy;
#endif
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		}

	      my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      if(All.SnapFormat == 1)
		{
		  if(dummy == 256)
		    swap_file = 8;
		  else
		    swap_file = dummy;
		}
#endif
	      my_fread(&header, sizeof(header), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_header();
#endif
	      my_fread(&dummy, sizeof(dummy), 1, fd);
	    }
	  fclose(fd);

#ifdef HAVE_HDF5
	  if(All.ICFormat == 3)
	    read_header_attributes_in_hdf5(buf1);
#endif
	  header.num_files = 1;
	}
    }

#ifdef AUTO_SWAP_ENDIAN_READIC
  MPI_Bcast(&swap_file, sizeof(swap_file), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
  MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(header.num_files > 0)
    return header.num_files;

  if(ThisTask == 0)
    {
      printf("\nCan't find initial conditions file.");
      printf("neither as '%s'\nnor as '%s'\n", buf, buf1);
      fflush(stdout);
    }

  endrun(0);
  return 0;
}



/*! This function assigns a certain number of files to processors, such that
 *  each processor is exactly assigned to one file, and the number of cpus per
 *  file is as homogenous as possible. The number of files may at most be
 *  equal to the number of processors.
 */
void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask,
		     int *filenr, int *master, int *last)
{
  int ntask, filesleft, filesright, tasksleft, tasksright;

  if(nfiles > 1)
    {
      ntask = lasttask - firsttask + 1;

      filesleft = (int) ((((double) (ntask / 2)) / ntask) * nfiles);
      if(filesleft <= 0)
	filesleft = 1;
      if(filesleft >= nfiles)
	filesleft = nfiles - 1;

      filesright = nfiles - filesleft;

      tasksleft = ntask / 2;
      tasksright = ntask - tasksleft;

      distribute_file(filesleft, firstfile, firsttask, firsttask + tasksleft - 1, filenr, master, last);
      distribute_file(filesright, firstfile + filesleft, firsttask + tasksleft,
		      lasttask, filenr, master, last);
    }
  else
    {
      if(ThisTask >= firsttask && ThisTask <= lasttask)
	{
	  *filenr = firstfile;
	  *master = firsttask;
	  *last = lasttask;
	}
    }
}



#ifdef HAVE_HDF5
void read_header_attributes_in_hdf5(char *fname)
{
  hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;

  hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hdf5_headergrp = H5Gopen(hdf5_file, "/Header");

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, header.npart);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total_HighWord");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "MassTable");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "ExpansionFactor");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "BoxSize");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumFilesPerSnapshot");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
  H5Aclose(hdf5_attribute);

  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);
}
#endif
