#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef BG_SFR

#include "allvars.h"
#include "proto.h"
#include "bg_proto.h"
#include "bg_vars.h"
#include "bg_yields.h"
#ifdef BG_DUST
#include "bg_dust.h"
#endif

#define EXTRASPACEFAC 2

static struct metaldata_in
{
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyFloat Hsml;

#ifdef BG_SNII_KINETIC_FEEDBACK
  MyFloat StarMass;
  MyDouble NgbMassSum;
#endif

#ifdef BG_DOUBLE_IMF
  MyFloat BirthDensity;
#endif

#ifdef BG_STELLAR_EVOLUTION
  MyFloat Metals[BG_NELEMENTS];
  MyFloat MetalMass;
#ifdef BG_DUST
  MyFloat DustMass;
#endif
  MyFloat SolidAngleWeightSum;
  MyFloat Energy;
#ifdef BG_SNIA_IRON
  MyFloat IronFromSNIa;
#endif
#endif

  unsigned int ID;
  int NodeList[NODELISTLENGTH];
} *MetalDataIn, *MetalDataGet;


static struct solidangledata_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int NodeList[NODELISTLENGTH];
}
 *SolidAngleDataIn, *SolidAngleDataGet;


static struct solidangledata_out
{
  MyDouble NgbMassSum;
  MyFloat SolidAngleWeightSum;
}
 *SolidAngleDataResult, *SolidAngleDataOut;




/*
 * These routines do metal-enrichment and wind formation for
 * the STELLA project on BlueGene.
 */

#ifdef BG_STELLAR_EVOLUTION
MyFloat MetalsReleased[BG_NELEMENTS];

MyFloat MetalMassReleased;
MyFloat DustMassReleased;

MyFloat EnergyReleased;		/* Energy from SNIa */

#ifdef BG_SNIA_IRON
MyFloat IronFromSNIa;		/* Iron from SNIa */
#endif
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK
MyFloat StarMass;
#endif


int bg_enrich_ngb_treefind(MyDouble searchcenter[3], MyFloat hsml, int target,
			   int *startnode, int mode, int *nexport, int *nsend_local)
{
  int numngb, no, p, task, nexport_save;

  struct NODE *current;

  MyDouble dx, dy, dz, dist;

#ifdef PERIODIC
  MyDouble xtmp;
#endif


  nexport_save = *nexport;

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(P[p].Type > 0)
	    continue;

	  dist = hsml;
	  dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  Ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(123127);

	      if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
		{
		  Exportflag[task] = target;
		  Exportnodecount[task] = NODELISTLENGTH;
		}

	      if(Exportnodecount[task] == NODELISTLENGTH)
		{
		  if(*nexport >= All.BunchSize)
		    {
		      /* out if buffer space. Need to discard work for this particle and interrupt */
		      *nexport = nexport_save;
		      if(nexport_save == 0)
			endrun(14003);	/* in this case, the buffer is too small to process even a single particle */
		      for(task = 0; task < NTask; task++)
			nsend_local[task] = 0;
		      for(no = 0; no < nexport_save; no++)
			nsend_local[DataIndexTable[no].Task]++;
		      return -1;
		    }
		  Exportnodecount[task] = 0;
		  Exportindex[task] = *nexport;
		  DataIndexTable[*nexport].Task = task;
		  DataIndexTable[*nexport].Index = target;
		  DataIndexTable[*nexport].IndexGet = *nexport;

#ifdef BG_STELLAR_EVOLUTION
		  MetalDataGet[*nexport].Energy = EnergyReleased;
		  MetalDataGet[*nexport].MetalMass = MetalMassReleased;
#ifdef BG_DUST
                  MetalDataGet[*nexport].DustMass = DustMassReleased;
#endif

		  memcpy(MetalDataGet[*nexport].Metals, MetalsReleased, sizeof(MyFloat) * BG_NELEMENTS);
#ifdef BG_SNIA_IRON
		  MetalDataGet[*nexport].IronFromSNIa = IronFromSNIa;
#endif
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK
		  MetalDataGet[*nexport].StarMass = StarMass;
#endif
		  *nexport = *nexport + 1;
		  nsend_local[task]++;
		}

	      DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
		DomainNodeIndex[no - (All.MaxPart + MaxNodes)];

	      if(Exportnodecount[task] < NODELISTLENGTH)
		DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;

	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
		  return numngb;
		}
	    }

	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;;
	  dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return numngb;
}


/* Main driver routine that deals with neighbour finding and parallelization 
   (export to other CPUs if needed) */
void bg_enrich(void)
{
  int i, j, ndone, ndone_flag, dummy, iter;

  int ngrp, sendTask, recvTask, nexport, nimport;

  int place, getplace;

  int imax1, imax2;

  double hubble_a, time_hubble_a, a3inv, tstart, tend, timestellarevol = 0;

  double dt, dtime, age_of_star_in_Gyr, age_of_star_in_Gyr_begstep, dtime_in_Gyr;

  MPI_Status status;

  /*
#ifdef BG_STELLAR_EVOLUTION
  double NumberOfSNIa;
#endif
  */

#if defined(BG_STELLAR_EVOLUTION) && defined(BG_VERBOSE)
  double TotalMetalsReleased[BG_NELEMENTS], GlobalMetalsReleased[BG_NELEMENTS];

  double TotalEnergyReleased, TotalMetalMassReleased;

  double GlobalEnergyReleased, GlobalMetalMassReleased;

  for(i = 0; i < BG_NELEMENTS; i++)
    TotalMetalsReleased[i] = 0;

  TotalEnergyReleased = 0;
  TotalMetalMassReleased = 0;
#endif

#ifdef BG_SNIA_IRON
  IronFromSNIa = 0;
#endif

  if(ThisTask == 0)
    {
      printf("Start enrichment...\n");
      fflush(stdout);
    }

#ifdef BG_STELLAR_EVOLUTION
  for(i = 0; i < BG_NELEMENTS; i++)
    MetalsReleased[i] = 0;

  EnergyReleased = MetalMassReleased = DustMassReleased = 0;
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK
  ExpectedMassLoading = ActualMassLoading = 0;
#endif

  if(All.ComovingIntegrationOn)	/* Factors for comoving integration of hydro */
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = All.Hubble * sqrt(All.Omega0 /
				   (All.Time * All.Time * All.Time)
				   + (1 - All.Omega0 - All.OmegaLambda) /
				   (All.Time * All.Time) + All.OmegaLambda);
      time_hubble_a = All.Time * hubble_a;
    }
  else
    a3inv = time_hubble_a = 1;

#ifdef BG_STELLAR_EVOLUTION
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0)
	SphP[i].OldMass = P[i].Mass;
    }
#endif

  /* compute the solid angle weight sum */
  bg_enrich_determine_weights();

  if(ThisTask == 0)
    {
      printf("Weight are updated, distributing metals...\n");
      fflush(stdout);
    }

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     (1 + EXTRASPACEFAC) * sizeof(struct metaldata_in)));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));

  iter = 0;
  i = FirstActiveParticle;	/* begin with this index */

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */

      MetalDataGet =
	(struct metaldata_in *) mymalloc(EXTRASPACEFAC * All.BunchSize * sizeof(struct metaldata_in));

      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	if(P[i].Type == 4)
	  {
	    /* the actual time-step, in Gadget units */
	    dt = TISTEP(P[i].TimeBin) * All.Timebase_interval;

	    if(All.ComovingIntegrationOn)
	      dtime = All.Time * dt / time_hubble_a;
	    else
	      dtime = dt;

	    dtime *= All.UnitTime_in_s / All.HubbleParam;	/* convert to  seconds */
	    dtime_in_Gyr = dtime / SEC_PER_MEGAYEAR / 1000;

	    age_of_star_in_Gyr = bg_get_elapsed_time(StarP[P[i].StarID].StarBirthTime, All.Time, 1);	/* note: this is valid for a flat universe! */
	    age_of_star_in_Gyr_begstep = age_of_star_in_Gyr - dtime_in_Gyr;

	    /* this function encapsulates the stellar evolution. It returns how much
	       mass in each element is returned during stellar evolution of the stellar
	       population between [age_of_star_in_Myr, age_of_star_in_Myr + dtime].
	       It also returns the amount of energy released by SNIa */

	    tstart = second();

	    /* global variables set to zero for each star particle */
	    NumberOfSNIa = NumberOfSNII = DustMassReleased = EnergyReleased = 0;

#ifdef BG_STELLAR_EVOLUTION
#ifdef BG_SNIA_IRON
	    bg_stellar_evolution(age_of_star_in_Gyr_begstep, dtime_in_Gyr, i, MetalsReleased, &MetalMassReleased,
				 &DustMassReleased, &IronFromSNIa, &EnergyReleased);
#else
	    bg_stellar_evolution(age_of_star_in_Gyr_begstep, dtime_in_Gyr, i, MetalsReleased, &MetalMassReleased,
				 &DustMassReleased, &EnergyReleased);
#endif

	    tend = second();
	    timestellarevol += timediff(tstart, tend);

	    if(dtime_in_Gyr > 0)
	      StarP[P[i].StarID].SNIaRate = NumberOfSNIa / All.HubbleParam / (dtime_in_Gyr * 1.e9);
#endif

	    /* now we check whether the star is eligible to perform feedback */
#ifdef BG_SNII_KINETIC_FEEDBACK
	    if(age_of_star_in_Gyr_begstep <= All.SNII_WindDelay_YR * 1.0e-9 &&
	       age_of_star_in_Gyr > All.SNII_WindDelay_YR * 1.0e-9)
	      {
		StarMass = StarP[P[i].StarID].InitialMass;	/* all SNII went off, the star can perform feedback */
	      }
	    else
	      StarMass = 0;	/* this signals that the star is not allowed to pruduce wind */
#endif

/* 	    /\* DEBUG *\/ */
/* 	    printf("[TASK %d] ID = %d, Type = %d, NumberOfSNIa = %f, NumberOfSNII = %f\n", */
/* 		    ThisTask, P[i].ID, P[i].Type, NumberOfSNIa, NumberOfSNII); */
/* 	    /\* DEBUG *\/ */

	    /* distribute ejected mass and SNIa energy and possibly launch the wind */
	    if(bg_enrich_evaluate(i, 0, &nexport, Send_count) < 0)
	      break;

#if defined(BG_STELLAR_EVOLUTION) && defined(BG_VERBOSE)
	    int iel;

	    for(iel = 0; iel < BG_NELEMENTS; iel++)
	      TotalMetalsReleased[iel] += MetalsReleased[iel];

	    TotalEnergyReleased += EnergyReleased;
	    TotalMetalMassReleased += MetalMassReleased;
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK
	    if(age_of_star_in_Gyr_begstep <= All.SNII_WindDelay_YR * 1.0e-9 &&
	       age_of_star_in_Gyr > All.SNII_WindDelay_YR * 1.0e-9)
	      {
#ifdef BG_DOUBLE_IMF
		if(StarP[P[i].StarID].GasDensity * a3inv > All.IMF_PhysDensThresh)
		  ExpectedMassLoading += StarMass * All.SNII_WindMassLoading1;
		else
		  ExpectedMassLoading += StarMass * All.SNII_WindMassLoading;
#else
		ExpectedMassLoading += StarMass * All.SNII_WindMassLoading;
#endif
	      }
#endif
	  }

#ifdef MYSORT
      mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#else
      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#endif


      MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);


      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      if(nimport > EXTRASPACEFAC * All.BunchSize)
	{
	  printf("On task=%d, we need to do a realloc unfortunately... nimport=%d  fac*Bunchsize=%d\n",
		 ThisTask, nimport, EXTRASPACEFAC * All.BunchSize);
	  fflush(stdout);
	  MetalDataGet =
	    (struct metaldata_in *) myrealloc(MetalDataGet,
					      IMAX(nimport,
						   EXTRASPACEFAC * All.BunchSize) *
					      sizeof(struct metaldata_in));
	}

      MetalDataIn = (struct metaldata_in *) mymalloc(nexport * sizeof(struct metaldata_in));


      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  getplace = DataIndexTable[j].IndexGet;

	  MetalDataIn[j].Pos[0] = P[place].Pos[0];
	  MetalDataIn[j].Pos[1] = P[place].Pos[1];
	  MetalDataIn[j].Pos[2] = P[place].Pos[2];

	  MetalDataIn[j].Vel[0] = P[place].Vel[0];
	  MetalDataIn[j].Vel[1] = P[place].Vel[1];
	  MetalDataIn[j].Vel[2] = P[place].Vel[2];

	  MetalDataIn[j].Hsml = PPP[place].Hsml;
	  MetalDataIn[j].ID = P[place].ID;

#ifdef BG_STELLAR_EVOLUTION
	  MetalDataIn[j].Energy = MetalDataGet[getplace].Energy;
	  MetalDataIn[j].SolidAngleWeightSum = StarP[P[place].StarID].SolidAngleWeightSum;
#ifdef BG_SNIA_IRON
	  MetalDataIn[j].IronFromSNIa = MetalDataGet[getplace].IronFromSNIa;
#endif
	  memcpy(MetalDataIn[j].Metals, MetalDataGet[getplace].Metals, sizeof(MyFloat) * (BG_NELEMENTS));

	  MetalDataIn[j].MetalMass = MetalDataGet[getplace].MetalMass;
#ifdef BG_DUST
          MetalDataIn[j].DustMass = MetalDataGet[getplace].DustMass;
#endif
#endif /* BG_STELLAR_EVOLUTION */

#ifdef BG_SNII_KINETIC_FEEDBACK
	  MetalDataIn[j].StarMass = MetalDataGet[getplace].StarMass;
	  MetalDataIn[j].NgbMassSum = StarP[P[place].StarID].NgbMassSum;
#endif

#ifdef BG_DOUBLE_IMF
	  MetalDataIn[j].BirthDensity = StarP[P[place].StarID].GasDensity;
#endif
	  memcpy(MetalDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	}


      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&MetalDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct metaldata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &MetalDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct metaldata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
		}
	    }
	}


      myfree(MetalDataIn);

      for(j = 0; j < nimport; j++)
	bg_enrich_evaluate(j, 1, &dummy, &dummy);

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      myfree(MetalDataGet);

      iter++;
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);


  CPU_Step[CPU_ENRICHSTEVOL] += timestellarevol;


#ifdef BG_SNII_KINETIC_FEEDBACK
  ActualMassLoading2ndLoop = bg_snii_kinetic_feedback_second_loop();
#endif

#if defined(BG_STELLAR_EVOLUTION) && defined(BG_VERBOSE)
  MPI_Reduce(TotalMetalsReleased, GlobalMetalsReleased, BG_NELEMENTS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&TotalMetalMassReleased, &GlobalMetalMassReleased, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&TotalEnergyReleased, &GlobalEnergyReleased, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("Metals produced:\n");
      for(i = 0; i < BG_NELEMENTS; i++)
	{
	  GlobalMetalsReleased[i] *= All.UnitMass_in_g / SOLAR_MASS / All.HubbleParam;
	  printf(" element = %s total amount [M_sun] = %g\n", ElementNames[i], GlobalMetalsReleased[i]);
	}

      /* print global metal produced in a more convenient format */
      printf("--Yields--[Msun]-- %f", All.Time);
      for(i = 0; i < BG_NELEMENTS; i++)
	printf(" %g", GlobalMetalsReleased[i]);
      printf("\n");

      printf("Metal mass released [Msun]: %e\n",
	     GlobalMetalMassReleased * All.UnitMass_in_g / SOLAR_MASS / All.HubbleParam);
      printf("SNIa energy released [erg]: %e\n",
	     GlobalEnergyReleased * All.UnitEnergy_in_cgs / All.HubbleParam);
    }
#endif

  if(TimeBinActive[16])		/* we only do full logs on sufficiently coarse steps */
    {
#ifdef BG_STELLAR_EVOLUTION
      /* log metallicity data */
      output_metallicity_data();

      /* log SNIa data */
      output_snia_data();
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK
      /* log wind data */
      output_wind_data();
#endif
    }
}


/*! This function represents the core of the enrichment calculation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */

int bg_enrich_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int startnode;

  unsigned int id;

  double h, h2;

  double ascale, a3inv, redshift;

#if defined(BG_STELLAR_EVOLUTION) || defined(BG_SNII_KINETIC_FEEDBACK)
  int j, k, n;

  int numngb_inbox, listindex = 0;

  double dx, dy, dz, r, r2;
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK
  double mass_star, ngb_mass;
#endif

#ifdef BG_STELLAR_EVOLUTION
  double dmZPart, mZPart, omega_frac, dmPart;

  MyFloat *metals, energy, metal_mass;
#ifdef BG_DUST
  double dmDPart, mDPart;

  MyFloat dust_mass;
#endif

  double solidAngleWeightSum;

  double OldEnergy, NewEnergy, dmax1, dmax2;

#ifdef BG_SNIA_IRON
  double iron_from_snia;
#endif

#ifdef BG_DOUBLE_IMF
  MyFloat birth_density;
#endif
#endif /* BG_STELLAR_EVOLUTION */

  MyDouble *pos, *vel;

  int static first_call = 0, HydrogenIndex, HeliumIndex;

  if(first_call == 0)
    {
      HydrogenIndex = element_index("Hydrogen");
      HeliumIndex = element_index("Helium");

      first_call = 1;
    }

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = P[target].Vel;
      h = PPP[target].Hsml;
#ifdef BG_STELLAR_EVOLUTION
      solidAngleWeightSum = StarP[P[target].StarID].SolidAngleWeightSum;
      energy = EnergyReleased;
      metals = MetalsReleased;
      metal_mass = MetalMassReleased;
#ifdef BG_DUST
      dust_mass = DustMassReleased;
#endif
#ifdef BG_DOUBLE_IMF
      birth_density = StarP[P[target].StarID].GasDensity;
#endif
#ifdef BG_SNIA_IRON
      iron_from_snia = IronFromSNIa;
#endif
#endif /* BG_STELLAR_EVOLUTION */

#ifdef BG_SNII_KINETIC_FEEDBACK
      mass_star = StarMass;
      ngb_mass = StarP[P[target].StarID].NgbMassSum;
#endif

      id = P[target].ID;
    }
  else
    {
      pos = MetalDataGet[target].Pos;
      vel = MetalDataGet[target].Vel;
      h = MetalDataGet[target].Hsml;
#ifdef BG_STELLAR_EVOLUTION
      solidAngleWeightSum = MetalDataGet[target].SolidAngleWeightSum;
      energy = MetalDataGet[target].Energy;
      metals = MetalDataGet[target].Metals;
      metal_mass = MetalDataGet[target].MetalMass;
#ifdef BG_DUST
      dust_mass = MetalDataGet[target].DustMass;
#endif
#ifdef BG_DOUBLE_IMF
      birth_density = MetalDataGet[target].BirthDensity;
#endif
#ifdef BG_SNIA_IRON
      iron_from_snia = MetalDataGet[target].IronFromSNIa;
#endif
#endif /* BG_STELLAR_EVOLUTION */

#ifdef BG_SNII_KINETIC_FEEDBACK
      mass_star = MetalDataGet[target].StarMass;
      ngb_mass = MetalDataGet[target].NgbMassSum;
#endif

      id = MetalDataGet[target].ID;
    }

  if(All.ComovingIntegrationOn)
    {
      ascale = All.Time;
      a3inv = 1 / (All.Time * All.Time * All.Time);
      redshift = 1 / All.Time - 1;	/* Joop fixed this. */
    }
  else
    {
      ascale = a3inv = 1;
      redshift = 0;
    }

  h2 = h * h;

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = MetalDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }


#if defined(BG_STELLAR_EVOLUTION) || defined(BG_SNII_KINETIC_FEEDBACK)
  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb_inbox = bg_enrich_ngb_treefind(&pos[0], h, target, &startnode, mode, nexport, nsend_local);

	  if(numngb_inbox < 0)
	    return -1;

	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = Ngblist[n];

	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	      if(dx > boxHalf_X)
		dx -= boxSize_X;
	      if(dx < -boxHalf_X)
		dx += boxSize_X;
	      if(dy > boxHalf_Y)
		dy -= boxSize_Y;
	      if(dy < -boxHalf_Y)
		dy += boxSize_Y;
	      if(dz > boxHalf_Z)
		dz -= boxSize_Z;
	      if(dz < -boxHalf_Z)
		dz += boxSize_Z;
#endif
	      r2 = dx * dx + dy * dy + dz * dz;

	      if(r2 < h2)
		{
		  r = sqrt(r2);

#ifdef BG_STELLAR_EVOLUTION
		  /* compute solid angle fraction */
		  omega_frac = bg_evaluate_solid_angle(j, r, h) / solidAngleWeightSum;

#ifdef BG_SNIA_IRON
		  SphP[j].IronFromSNIa += iron_from_snia * omega_frac;
#endif

		  /* original mass in metals */
		  mZPart = P[j].Mass * SphP[j].Metallicity;

		  /* metal mass variation */
		  dmZPart = metal_mass * omega_frac;

		  /* update total metal mass */
		  SphP[j].Metallicity = mZPart + dmZPart;

#ifdef BG_DUST
		  /* dust destruction */
		  if(SphP[j].Dusticity > 0)
		    dust_destruction(j);

                  /* original mass in dust */
                  mDPart = P[j].Mass * SphP[j].Dusticity;

                  /* dust mass variation */
                  dmDPart = dust_mass * omega_frac;

                  /* update total dust mass */
                  SphP[j].Dusticity = mDPart + dmDPart;

#ifdef BG_DUST_DESTRUCTION_SUBLIMATION
		  /* update average dust grain size */
		  SphP[j].DustGrainSize = (SphP[j].DustGrainSize * mDPart + All.InitDustGrainSize_nM * dmDPart) / (mDPart + dmDPart);
#endif
#endif

		  if(SphP[j].Metallicity < 0)
		    printf("wrong Z (1): %e %e %e %e %e %d\n", SphP[j].Metallicity, P[j].Mass, metal_mass,
			   dmZPart, omega_frac, mode);

#ifdef BG_Z_WEIGHTED_REDSHIFT
		  if(dmZPart > 0)
		    SphP[j].MetallicityWeightedRedshift =
		      (redshift * dmZPart +
		       SphP[j].MetallicityWeightedRedshift * mZPart) / (dmZPart + mZPart);
#endif

		  /* add H and He to metal mass released to find the total mass variation */
		  dmPart = metal_mass + metals[HydrogenIndex] + metals[HeliumIndex];


		  /* update particle momentum */
		  for(k = 0; k < 3; k++)
		    {
		      P[j].Vel[k] = (P[j].Mass * P[j].Vel[k] + dmPart * omega_frac * vel[k]) / (P[j].Mass + dmPart * omega_frac);

		      SphP[j].VelPred[k] = (P[j].Mass * SphP[j].VelPred[k] + dmPart * omega_frac * vel[k]) / (P[j].Mass + dmPart * omega_frac);
		    }

		  /* update element masses */
		  for(k = 0; k < BG_NELEMENTS; k++)
		    SphP[j].Metals[k] += metals[k] * omega_frac;

		  /* update particle mass */
		  P[j].Mass += dmPart * omega_frac;

		  /* convert back to metallicity */
		  SphP[j].Metallicity /= P[j].Mass;
#ifdef BG_DUST
                  /* convert back to metallicity */
                  SphP[j].Dusticity /= P[j].Mass;
#endif

		  /*
		   * Add energy relased to old energy to obtain new energy.
		   * Then compute the corresponding new entropy
		   */
		  if(All.SNIa_EnergyTransferOn == 1)
		    {
		      OldEnergy = DMAX(All.MinEgySpec, SphP[j].Entropy / GAMMA_MINUS1 * pow(SphP[j].d.Density * a3inv, GAMMA_MINUS1));	/* E/m */
		      NewEnergy = OldEnergy * SphP[j].OldMass + energy * omega_frac;	/* E */
		      SphP[j].Entropy *= (NewEnergy / SphP[j].OldMass) / OldEnergy;	/* update entropy */
		    }
#endif /* BG_STELLAR_EVOLUTION */

#ifdef BG_SNII_KINETIC_FEEDBACK
#ifdef BG_DOUBLE_IMF
		  bg_snii_kinetic_feedback(j, id, &pos[0], mass_star, ngb_mass, birth_density);
#else
		  bg_snii_kinetic_feedback(j, id, &pos[0], mass_star, ngb_mass);
#endif
#endif
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = MetalDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }
#endif /* BG_STELLAR_EVOLUTION || BG_KINETIC_FEEDBACK */

  return 0;
}





/*
 * ----------------------------------------------------------------------
 * Output the amount of metals in gas particles, star forming particles,
 * and star particles.
 * ----------------------------------------------------------------------
 */

#ifdef BG_STELLAR_EVOLUTION
void output_metallicity_data(void)
{
  int i, j;

  double el_mass_gas[BG_NELEMENTS + 1], el_mass_stars[BG_NELEMENTS + 1],
    el_mass_sf[BG_NELEMENTS + 1], el_mass_total[BG_NELEMENTS + 1], global_el_mass[BG_NELEMENTS + 1];

  double a3inv;


  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  for(i = 0; i < BG_NELEMENTS + 1; i++)
    {
      el_mass_gas[i] = 0.0;
      el_mass_stars[i] = 0.0;
      el_mass_sf[i] = 0.0;
    }

  for(j = 0; j < NumPart; j++)
    {
      if(P[j].Type == 0)	/* gas particle */
	{
	  for(i = 0; i < BG_NELEMENTS; i++)
	    el_mass_gas[i] += SphP[j].Metals[i];

	  el_mass_gas[BG_NELEMENTS] += P[j].Mass;

	  if(SphP[j].Sfr > 0)	/* star forming particle */
	    {
	      for(i = 0; i < BG_NELEMENTS; i++)
		el_mass_sf[i] += SphP[j].Metals[i];

	      el_mass_sf[BG_NELEMENTS] += P[j].Mass;
	    }
	}
      else if(P[j].Type == 4)	/* star particle */
	{
	  for(i = 0; i < BG_NELEMENTS; i++)
	    el_mass_stars[i] += StarP[P[j].StarID].Metals[i];

	  el_mass_stars[BG_NELEMENTS] += P[j].Mass;
	}
    }

  for(i = 0; i < BG_NELEMENTS; i++)
    el_mass_total[i] = el_mass_gas[i] + el_mass_stars[i];

  el_mass_total[BG_NELEMENTS] = el_mass_gas[BG_NELEMENTS] + el_mass_stars[BG_NELEMENTS];

  /* output total metal mass in gas particles */
  MPI_Reduce(el_mass_gas, global_el_mass, BG_NELEMENTS + 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdMetGas, "%e", All.Time);
      /* convert to solar masses */
      for(i = 0; i < BG_NELEMENTS + 1; i++)
	{
	  global_el_mass[i] *= All.UnitMass_in_g / SOLAR_MASS;
	  fprintf(FdMetGas, " %e", global_el_mass[i]);
	}
      fprintf(FdMetGas, "\n");
      fflush(FdMetGas);
    }

  /* output total metal mass in star particles */
  MPI_Reduce(el_mass_stars, global_el_mass, BG_NELEMENTS + 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdMetStars, "%e", All.Time);
      /* convert to solar masses */
      for(i = 0; i < BG_NELEMENTS + 1; i++)
	{
	  global_el_mass[i] *= All.UnitMass_in_g / SOLAR_MASS;
	  fprintf(FdMetStars, " %e", global_el_mass[i]);
	}
      fprintf(FdMetStars, "\n");
      fflush(FdMetStars);
    }

  /* output total metal mass in star forming particles */
  MPI_Reduce(el_mass_sf, global_el_mass, BG_NELEMENTS + 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdMetSF, "%e", All.Time);
      /* convert to solar masses */
      for(i = 0; i < BG_NELEMENTS + 1; i++)
	{
	  global_el_mass[i] *= All.UnitMass_in_g / SOLAR_MASS;
	  fprintf(FdMetSF, " %e", global_el_mass[i]);
	}
      fprintf(FdMetSF, "\n");
      fflush(FdMetSF);
    }

  /* output total metal mass in all particles */
  MPI_Reduce(el_mass_total, global_el_mass, BG_NELEMENTS + 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdMetTot, "%e", All.Time);
      /* convert to solar masses */
      for(i = 0; i < BG_NELEMENTS + 1; i++)
	{
	  global_el_mass[i] *= All.UnitMass_in_g / SOLAR_MASS;
	  fprintf(FdMetTot, " %e", global_el_mass[i]);
	}
      fprintf(FdMetTot, "\n");
      fflush(FdMetTot);
    }
}
#endif

/*
 * ----------------------------------------------------------------------
 * Output wind expected and actual mass loading
 * ----------------------------------------------------------------------
 */

#ifdef BG_SNII_KINETIC_FEEDBACK
void output_wind_data(void)
{
  double expected_mass_loading, actual_mass_loading, actual_mass_loading_2nd_loop;

  double global_expected_mass_loading, global_actual_mass_loading, global_actual_mass_loading_2nd_loop;

  expected_mass_loading = ExpectedMassLoading;
  actual_mass_loading = ActualMassLoading;
  actual_mass_loading_2nd_loop = ActualMassLoading2ndLoop;

  MPI_Reduce(&expected_mass_loading, &global_expected_mass_loading, 1, MPI_DOUBLE, MPI_SUM, 0,
	     MPI_COMM_WORLD);
  MPI_Reduce(&actual_mass_loading, &global_actual_mass_loading, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&actual_mass_loading_2nd_loop, &global_actual_mass_loading_2nd_loop, 1, MPI_DOUBLE, MPI_SUM, 0,
	     MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdKineticFeedback, "%e %e %e %e\n", All.Time, global_expected_mass_loading,
	      global_actual_mass_loading, global_actual_mass_loading + global_actual_mass_loading_2nd_loop);
      fflush(FdKineticFeedback);
    }
}
#endif

/*
 * ----------------------------------------------------------------------
 * Output SNIa number and rate
 * ----------------------------------------------------------------------
 */

#ifdef BG_STELLAR_EVOLUTION
void output_snia_data()
{
  int i;

  double TotalSNIaRate, GlobalSNIaRate;

  for(i = 0, TotalSNIaRate = 0; i < N_star; i++)
    TotalSNIaRate += StarP[i].SNIaRate;

  MPI_Reduce(&TotalSNIaRate, &GlobalSNIaRate, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdSNIa, "%e %e\n", All.Time, GlobalSNIaRate);
      fflush(FdSNIa);
    }
}
#endif

void bg_enrich_determine_weights(void)
{
  int i, j, ndone, ndone_flag, dummy;

  int ngrp, sendTask, recvTask, nexport, nimport, place;

  MPI_Status status;


  /* allocate buffers to arrange communication */
  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct solidangledata_in) +
					     sizeof(struct solidangledata_out) +
					     sizemax(sizeof(struct solidangledata_in),
						     sizeof(struct solidangledata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));


  i = FirstActiveParticle;	/* begin with this index */

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */

      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	if(P[i].Type == 4)
	  {
	    if(bg_enrich_evaluate_weightsum(i, 0, &nexport, Send_count) < 0)
	      break;
	  }

      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);

      MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      SolidAngleDataGet = (struct solidangledata_in *) mymalloc(nimport * sizeof(struct solidangledata_in));
      SolidAngleDataIn = (struct solidangledata_in *) mymalloc(nexport * sizeof(struct solidangledata_in));

      /* prepare particle data for export */

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  SolidAngleDataIn[j].Pos[0] = P[place].Pos[0];
	  SolidAngleDataIn[j].Pos[1] = P[place].Pos[1];
	  SolidAngleDataIn[j].Pos[2] = P[place].Pos[2];
	  SolidAngleDataIn[j].Hsml = PPP[place].Hsml;

	  memcpy(SolidAngleDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	}

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&SolidAngleDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct solidangledata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &SolidAngleDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct solidangledata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
		}
	    }
	}

      myfree(SolidAngleDataIn);

      SolidAngleDataResult =
	(struct solidangledata_out *) mymalloc(nimport * sizeof(struct solidangledata_out));
      SolidAngleDataOut = (struct solidangledata_out *) mymalloc(nexport * sizeof(struct solidangledata_out));

      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
	bg_enrich_evaluate_weightsum(j, 1, &dummy, &dummy);

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      /* get the result */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&SolidAngleDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct solidangledata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B,
			       &SolidAngleDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct solidangledata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
		}
	    }

	}

      /* add the result to the local particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

#ifdef BG_STELLAR_EVOLUTION
	  StarP[P[place].StarID].SolidAngleWeightSum += SolidAngleDataOut[j].SolidAngleWeightSum;
#endif
#ifdef BG_SNII_KINETIC_FEEDBACK
	  StarP[P[place].StarID].NgbMassSum += SolidAngleDataOut[j].NgbMassSum;
#endif
	}

      myfree(SolidAngleDataOut);
      myfree(SolidAngleDataResult);
      myfree(SolidAngleDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
}




int bg_enrich_evaluate_weightsum(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n;

  int startnode, numngb_inbox, listindex = 0;

  double h, h2, omega, ngb_mass;

  double dx, dy, dz, r, r2;

#ifdef BG_SNII_KINETIC_FEEDBACK
  double time_in_wind_in_Gyr;
#endif

  MyDouble *pos;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
    }
  else
    {
      pos = SolidAngleDataGet[target].Pos;
      h = SolidAngleDataGet[target].Hsml;
    }

  h2 = h * h;

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = SolidAngleDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  omega = ngb_mass = 0.0;

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb_inbox = ngb_treefind_variable(&pos[0], h, target, &startnode, mode, nexport, nsend_local);

	  if(numngb_inbox < 0)
	    return -1;

	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = Ngblist[n];

	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	      if(dx > boxHalf_X)
		dx -= boxSize_X;
	      if(dx < -boxHalf_X)
		dx += boxSize_X;
	      if(dy > boxHalf_Y)
		dy -= boxSize_Y;
	      if(dy < -boxHalf_Y)
		dy += boxSize_Y;
	      if(dz > boxHalf_Z)
		dz -= boxSize_Z;
	      if(dz < -boxHalf_Z)
		dz += boxSize_Z;
#endif
	      r2 = dx * dx + dy * dy + dz * dz;

	      if(r2 < h2)
		{
		  r = sqrt(r2);

#ifdef BG_STELLAR_EVOLUTION
		  omega += bg_evaluate_solid_angle(j, r, h);
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK
		  if(SphP[j].WindFlag < 0)
		    {
		      /*
		       * check whether the wind particle is still active computing
		       * the elapsed time since the particle was put in the wind
		       *
		       * note: this is valid for a flat universe (mode = 1)!
		       * use mode = 0 for any other cosmology
		       */
		      time_in_wind_in_Gyr = bg_get_elapsed_time(-(double) SphP[j].WindFlag, All.Time, 1);
		    }
		  else
		    time_in_wind_in_Gyr = 0;

		  /* if the particle is currently in the wind we skip it */
#ifdef BG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING
		  if(SphP[j].WindFlag == 0 || time_in_wind_in_Gyr > All.SNII_WindDecouplingTime_YR * 1.0e-9)
#else
		  if(SphP[j].WindFlag == 0)
#endif
		    {
#endif
		      ngb_mass += P[j].Mass;
#ifdef BG_SNII_KINETIC_FEEDBACK
		    }
#endif
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = SolidAngleDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  if(mode == 0)
    {
#ifdef BG_STELLAR_EVOLUTION
      StarP[P[target].StarID].SolidAngleWeightSum = omega;
#endif
#if defined(BG_SNII_KINETIC_FEEDBACK) || defined(BG_SNII_THERMAL_FEEDBACK)
      StarP[P[target].StarID].NgbMassSum = ngb_mass;
#endif
    }
  else
    {
#ifdef BG_STELLAR_EVOLUTION
      SolidAngleDataResult[target].SolidAngleWeightSum = omega;
#endif
#if defined(BG_SNII_KINETIC_FEEDBACK) || defined(BG_SNII_THERMAL_FEEDBACK)
      SolidAngleDataResult[target].NgbMassSum = ngb_mass;
#endif
    }
  return 0;
}


#ifdef BG_STELLAR_EVOLUTION
double bg_evaluate_solid_angle(int j, double r, double h)
{
  double u;

  double wk;

  double hinv, hinv3;

  hinv = 1 / h;
  hinv3 = hinv * hinv * hinv;

  u = r * hinv;

  if(u < 0.5)
    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
  else
    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

  return wk * SphP[j].OldMass / SphP[j].d.Density;
}
#endif

#endif
