#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "bg_proto.h"
#include "bg_vars.h"


#if defined(BLACK_HOLES) && defined(BH_THERMALFEEDBACK) 

#define EXTRASPACEFAC 2

#ifdef BG_THERMAL_FEEDBACK_MAKE_PARTICLES_ACTIVE
#define EXTRALENGTHFAC 8
#else
#define EXTRALENGTHFAC 1
#endif

static struct feedbackdata_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
#ifdef BH_THERMALFEEDBACK
  MyDouble NgbMassSum;
  int NumNgb;
  MyDouble BH_Energy;
#endif
  unsigned int ID;
  int NodeList[NODELISTLENGTH];
} *FeedbackDataIn, *FeedbackDataGet;


static struct ngbmassdata_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int NodeList[NODELISTLENGTH];
}
 *NgbMassDataIn, *NgbMassDataGet;


static struct ngbmassdata_out
{
  MyDouble NgbMassSum;
  int NumNgb;
}
 *NgbMassDataResult, *NgbMassDataOut;

#ifdef BH_DEBUG
double total_energy, Total_energy;
double total_energy_ready, Total_energy_ready;
double total_energy_used, Total_energy_used;
double total_energy_removed, Total_energy_removed;
int    min_numngb, Min_numngb;
int    max_numngb, Max_numngb;
double min_ngbmass, Min_ngbmass;
double max_ngbmass, Max_ngbmass;
double nactive, Nactive;
#endif

void bg_bh_ngbmass(void);

int bg_bh_ngbmass_evaluate(int target, int mode, int *nexport, int *nsend_local);

int bg_bh_thermal_feedback_ngb_treefind(MyDouble searchcenter[3], MyFloat hsml, int target,
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
		      /* out of buffer space. Need to discard work for this particle and interrupt */
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


void bg_bh_thermal_feedback(void)
{
  int i, j, ndone, ndone_flag, dummy, iter;

  int ngrp, sendTask, recvTask, nexport, nimport;

  int place, getplace;

  int imax1, imax2;

  double hubble_a, time_hubble_a, a3inv, critical_energy, u_to_temp_fac;

  MPI_Status status;


  if(All.ComovingIntegrationOn)
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

  /* Conversion between internal energy per unit mass and temperature for a fully ionized plasma
     of primordial composition */
  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1
    * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  /* compute the neighbour masses */
  bg_bh_ngbmass();

  if(ThisTask == 0)
    {
      printf("Ngb masses are updated, applying AGN feedback...\n");
      fflush(stdout);
    }

#ifdef BH_DEBUG
  total_energy         = 0.0;
  total_energy_used    = 0.0;
  total_energy_ready   = 0.0;
  total_energy_removed = 0.0;
  min_numngb           = 32000;
  max_numngb           = -32000;
  min_ngbmass          = 1e10;
  max_ngbmass          = -1e10;
  nactive              = 0;
#endif

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     (1 + EXTRASPACEFAC) * sizeof(struct feedbackdata_in)));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));

#ifdef BH_DEBUG
  i = FirstActiveParticle;	/* begin with this index */
  for(nexport = 0; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 5  && P[i].NgbMassSum > 0) {
      critical_energy = All.BlackHoleMinHeatTemp * (P[i].NgbMassSum/(float)P[i].NumNgb) * All.BlackHoleNumberOfNeighboursToHeat / u_to_temp_fac;
      total_energy += P[i].BH_Energy;
      if (P[i].BH_Energy > critical_energy) total_energy_ready += P[i].BH_Energy;
      if (P[i].NumNgb < min_numngb) min_numngb = P[i].NumNgb; 
      if (P[i].NumNgb > max_numngb) max_numngb = P[i].NumNgb; 
      if (P[i].NgbMassSum < min_ngbmass) min_ngbmass = P[i].NgbMassSum; 
      if (P[i].NgbMassSum > max_ngbmass) max_ngbmass = P[i].NgbMassSum; 
      nactive++;
    }

  MPI_Allreduce(&min_numngb, &Min_numngb, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&max_numngb, &Max_numngb, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&min_ngbmass, &Min_ngbmass, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&max_ngbmass, &Max_ngbmass, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&nactive, &Nactive, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (ThisTask == 0 && Nactive > 0) {
    printf("Min and max NumNgb  = %d %d\n",Min_numngb,Max_numngb);
    printf("Min and max NgbMass = %e %e\n",Min_ngbmass,Max_ngbmass);
    printf("Nactive = %f\n",Nactive);
  }
#endif

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

      FeedbackDataGet =
	(struct feedbackdata_in *) mymalloc(EXTRASPACEFAC * All.BunchSize * sizeof(struct feedbackdata_in));

      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	if(P[i].Type == 5 && P[i].NgbMassSum > 0)
	  {
	    /* apply AGN feedback */
	    if(bg_bh_thermal_feedback_evaluate(i, 0, &nexport, Send_count) < 0)
	      break;
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
	  FeedbackDataGet =
	    (struct feedbackdata_in *) myrealloc(FeedbackDataGet,
					      IMAX(nimport,
						   EXTRASPACEFAC * All.BunchSize) *
					      sizeof(struct feedbackdata_in));
	}


      FeedbackDataIn = (struct feedbackdata_in *) mymalloc(nexport * sizeof(struct feedbackdata_in));

      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  getplace = DataIndexTable[j].IndexGet;

	  FeedbackDataIn[j].Pos[0]     = P[place].Pos[0];
	  FeedbackDataIn[j].Pos[1]     = P[place].Pos[1];
	  FeedbackDataIn[j].Pos[2]     = P[place].Pos[2];
	  FeedbackDataIn[j].Hsml       = PPP[place].Hsml;
	  FeedbackDataIn[j].ID         = P[place].ID;
	  FeedbackDataIn[j].NgbMassSum = P[place].NgbMassSum;
	  FeedbackDataIn[j].NumNgb     = P[place].NumNgb;
	  FeedbackDataIn[j].BH_Energy  = P[place].BH_Energy;
 
	  memcpy(FeedbackDataIn[j].NodeList,
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
		  MPI_Sendrecv(&FeedbackDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct feedbackdata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &FeedbackDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct feedbackdata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
		}
	    }
	}


      myfree(FeedbackDataIn);

      for(j = 0; j < nimport; j++)
	bg_bh_thermal_feedback_evaluate(j, 1, &dummy, &dummy);

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      myfree(FeedbackDataGet);

      iter++;
    }
  while(ndone < NTask);

  /* Subtract used energy */
  i = FirstActiveParticle;	/* begin with this index */
  for(nexport = 0; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 5 && P[i].NgbMassSum > 0 && P[i].NumNgb > 0)
      {
	critical_energy = All.BlackHoleMinHeatTemp * (P[i].NgbMassSum/(float)P[i].NumNgb) * All.BlackHoleNumberOfNeighboursToHeat / u_to_temp_fac;
	if (P[i].BH_Energy > critical_energy) {
#ifdef BH_DEBUG
	  total_energy_removed += P[i].BH_Energy;
#endif
	  P[i].BH_Energy = 0.0;
	}
      }
      
  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

#ifdef BH_DEBUG
  MPI_Allreduce(&total_energy, &Total_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&total_energy_ready, &Total_energy_ready, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&total_energy_used, &Total_energy_used, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&total_energy_removed, &Total_energy_removed, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (ThisTask == 0) {
    if (Nactive > 0) {
      printf("BH feedback total energy in simulation           : %f\n",Total_energy);
      printf("BH feedback total energy ready to be distributed : %f\n",Total_energy_ready);
      printf("BH feedback total energy actually used           : %f\n",Total_energy_used);
      printf("BH feedback total energy removed from reservoir  : %f\n",Total_energy_removed);
    }
      if (Total_energy_ready > 0) {
	fprintf(FdBlackHolesEnergy,"%f %e %e %e %e\n",All.Time,Total_energy,Total_energy_ready,Total_energy_used,Total_energy_removed);
	fflush(FdBlackHolesEnergy);
    }
  }


#endif

  /* fflush(FdBlackHolesFeedback); */

}


int bg_bh_thermal_feedback_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n;

  unsigned int id;

  int startnode, numngb_inbox, listindex = 0;

  double h, h2;

  double dx, dy, dz, r2, ascale, a3inv;

  double OldEnergy, NewEnergy, dmax1, dmax2;

#ifdef VERBOSE_LOGFILES
  double OldTemperature;
#endif

#ifdef BG_THERMAL_FEEDBACK_MAKE_PARTICLES_ACTIVE
  double csnd, vel;
#endif

  double ngb_mass, bh_energy;

  int num_ngb;

  MyDouble *pos;

  double u_to_temp_fac, p;

  double critical_energy, critical_energy_per_unit_mass, nheat;


  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
      ngb_mass = P[target].NgbMassSum;
      num_ngb  = P[target].NumNgb;
      bh_energy = P[target].BH_Energy;
      id = P[target].ID;
    }
  else
    {
      pos = FeedbackDataGet[target].Pos;
      h = FeedbackDataGet[target].Hsml;
      ngb_mass = FeedbackDataGet[target].NgbMassSum;
      num_ngb  = FeedbackDataGet[target].NumNgb;
      bh_energy = FeedbackDataGet[target].BH_Energy;
      id = FeedbackDataGet[target].ID;
    }

  if(All.ComovingIntegrationOn)
    {
      ascale = All.Time;
      a3inv = 1 / (All.Time * All.Time * All.Time);
    }
  else
    {
      ascale = a3inv = 1;
    }


  if (num_ngb == 0) return 0;

  h2 = h * h;

  /* Conversion between internal energy per unit mass and temperature for a fully ionized plasma
     of primordial composition */
  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1
    * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = FeedbackDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }


  critical_energy_per_unit_mass = All.BlackHoleMinHeatTemp * All.BlackHoleNumberOfNeighboursToHeat / u_to_temp_fac;

  critical_energy = critical_energy_per_unit_mass * ngb_mass / num_ngb;

  if (bh_energy > critical_energy)
    {
      nheat = bh_energy / critical_energy;
      p = nheat / num_ngb;
    }
  else
    return 0;

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{

	  numngb_inbox = bg_bh_thermal_feedback_ngb_treefind(&pos[0], EXTRALENGTHFAC * h, target, &startnode, mode, nexport, nsend_local);

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
		  if (get_random_number(P[j].ID + id) < p)
		    {
#ifdef VERBOSE_LOGFILES
		      /* get old temperature for logging */
		      OldTemperature = bg_get_temperature(j);
#endif
		      /* current energy per unit mass */
		      OldEnergy = DMAX(All.MinEgySpec, SphP[j].Entropy / GAMMA_MINUS1 *
				       pow(SphP[j].d.Density * a3inv, GAMMA_MINUS1));
		      
		      /* updated energy per unit mass */
		      /* NewEnergy = OldEnergy + critical_energy / (ngb_mass / (float) num_ngb); */
		      NewEnergy = OldEnergy + critical_energy_per_unit_mass;
#ifdef BH_DEBUG
		      total_energy_used += critical_energy;
#endif
		      /* updated entropy */
		      SphP[j].Entropy *= NewEnergy / OldEnergy;

		      /* store energy */
		      /* SphP[j].i.Injected_Energy += critical_energy/(ngb_mass/(float)num_ngb) * P[j].Mass; */

#ifdef VERBOSE_LOGFILES
		      fprintf(FdBlackHolesFeedback, "%f %13d %13d %e %e %e %e %e %f %f %f %d\n",
			      All.Time, id, P[j].ID, OldEnergy * P[j].Mass,
			      NewEnergy * P[j].Mass, OldTemperature,
			      bg_get_temperature(j), P[j].Mass, critical_energy, nheat, p, mode);

		      fflush(FdBlackHolesFeedback);
#endif

#ifdef BG_THERMAL_FEEDBACK_MAKE_PARTICLES_ACTIVE
		      /* use a trick to shorten the timestep in the next call to timestep.c
			 change the MaxSignalVel to twice the sound speed given by the
			 entropy change of heated particles */
		      csnd = sqrt(GAMMA * SphP[j].Entropy * pow(SphP[j].d.Density, GAMMA_MINUS1));

		      vel = sqrt(SphP[j].VelPred[0] * SphP[j].VelPred[0] +
				 SphP[j].VelPred[1] * SphP[j].VelPred[1] +
				 SphP[j].VelPred[2] * SphP[j].VelPred[2]);

/* #ifdef ALTERNATIVE_VISCOUS_TIMESTEP */
/* 		      double dmax1, dmax2, dt_courant; */

/* 		      if(All.ComovingIntegrationOn) */
/* 			dt_courant = All.CourantFac * All.Time * DMAX(PPP[j].Hsml, All.SofteningTable[0]) / (pow(All.Time, 3 * (1 - GAMMA) / 2.0) * (2 * csnd + 3 * vel)); */
/* 		      else */
/* 			dt_courant = All.CourantFac * DMAX(PPP[j].Hsml, All.SofteningTable[0]) / (2 * csnd + 3 * vel); */

/* 		      if(dt_courant < 2 * All.CourantFac * SphP[j].MinViscousDt) */
/* 			SphP[j].MinViscousDt = dt_courant / (2 * All.CourantFac); */
/* #else */
		      if(2 * csnd + 3 * vel > SphP[j].MaxSignalVel) /* if the particle has already been touched skip it */
			SphP[j].MaxSignalVel = 2 * csnd + 3 * vel;
/* #endif */
		      if(!TimeBinActive[P[j].TimeBin])
			synchronize_particle(j);
#endif
		    }
		}
            }
        }

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = FeedbackDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  return 0;
}


void bg_bh_ngbmass(void)
{
  int i, j, ndone, ndone_flag, dummy;

  int ngrp, sendTask, recvTask, nexport, nimport, place;

  MPI_Status status;


  /* allocate buffers to arrange communication */
  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct ngbmassdata_in) +
					     sizeof(struct ngbmassdata_out) +
					     sizemax(sizeof(struct ngbmassdata_in),
						     sizeof(struct ngbmassdata_out))));
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
         if(P[i].Type == 5)
	  {
	    if(bg_bh_ngbmass_evaluate(i, 0, &nexport, Send_count) < 0)
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

      NgbMassDataGet = (struct ngbmassdata_in *) mymalloc(nimport * sizeof(struct ngbmassdata_in));
      NgbMassDataIn = (struct ngbmassdata_in *) mymalloc(nexport * sizeof(struct ngbmassdata_in));

      /* prepare particle data for export */

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  NgbMassDataIn[j].Pos[0] = P[place].Pos[0];
	  NgbMassDataIn[j].Pos[1] = P[place].Pos[1];
	  NgbMassDataIn[j].Pos[2] = P[place].Pos[2];
	  NgbMassDataIn[j].Hsml = PPP[place].Hsml;

	  memcpy(NgbMassDataIn[j].NodeList,
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
		  MPI_Sendrecv(&NgbMassDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct ngbmassdata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &NgbMassDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct ngbmassdata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
		}
	    }
	}

      myfree(NgbMassDataIn);

      NgbMassDataResult =
	(struct ngbmassdata_out *) mymalloc(nimport * sizeof(struct ngbmassdata_out));
      NgbMassDataOut = (struct ngbmassdata_out *) mymalloc(nexport * sizeof(struct ngbmassdata_out));

      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
	bg_bh_ngbmass_evaluate(j, 1, &dummy, &dummy);

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
		  MPI_Sendrecv(&NgbMassDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct ngbmassdata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B,
			       &NgbMassDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct ngbmassdata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
		}
	    }

	}

      /* add the result to the local particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  P[place].NgbMassSum += NgbMassDataOut[j].NgbMassSum;
	  P[place].NumNgb += NgbMassDataOut[j].NumNgb;
	}

      myfree(NgbMassDataOut);
      myfree(NgbMassDataResult);
      myfree(NgbMassDataGet);
    }
  while(ndone < NTask);


#ifdef BH_DEBUG
  i = FirstActiveParticle;
  for(nexport = 0; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 5 && P[i].NgbMassSum <= 0) {
      printf("BH Warning! Particle %u zero NgbMass Hsml=%f NumNgb=%d Mass=%f\n",P[i].ID,PPP[i].Hsml, P[i].NumNgb,P[i].Mass);
    }
#endif

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
}


int bg_bh_ngbmass_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n;

  int startnode, numngb_inbox, listindex = 0;

  double h, h2, ngb_mass;

  int num_ngb;

  double dx, dy, dz, r2;

  MyDouble *pos;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
    }
  else
    {
      pos = NgbMassDataGet[target].Pos;
      h = NgbMassDataGet[target].Hsml;
    }

  h2 = h * h;

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = NgbMassDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  ngb_mass = 0.0;
  num_ngb = 0;

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

	      if(r2 < h2) {
		ngb_mass += P[j].Mass;
		num_ngb++;
	      }
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = NgbMassDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  if(mode == 0) {
    P[target].NgbMassSum = ngb_mass;
    P[target].NumNgb = num_ngb;
  } else {
    NgbMassDataResult[target].NgbMassSum = ngb_mass;
    NgbMassDataResult[target].NumNgb = num_ngb;
  }

  return 0;
}
#endif

