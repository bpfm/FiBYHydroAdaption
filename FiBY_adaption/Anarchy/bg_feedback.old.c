#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "bg_proto.h"
#include "bg_vars.h"


#define EXTRASPACEFAC 2

#define EXTRALENGTHFAC 16

static struct feedbackdata_in
{
  MyDouble Pos[3];
  MyFloat Hsml;

  MyFloat StarMass;
  MyDouble NgbMassSum;

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
}
 *NgbMassDataResult, *NgbMassDataOut;


MyFloat StarMass;

void synchronize_particle(int j);


#if defined(BG_SNII_THERMAL_FEEDBACK) || (defined(BG_POPIII) && defined(BG_POPIII_THERMAL_FEEDBACK))

void bg_thermal_feedback(void)
{
  /* compute neighbours total mass */
  bg_ngbmass();

  if(ThisTask == 0)
    {
      printf("Ngb masses are updated\n");
      fflush(stdout);
    }

#ifdef BG_SNII_THERMAL_FEEDBACK
  if(ThisTask == 0)
    {
      printf("Applying SNII feedback...\n");
      fflush(stdout);
    }

  bg_snii_thermal_feedback();
#endif

#if defined(BG_POPIII) && defined(BG_POPIII_THERMAL_FEEDBACK)
  if(ThisTask == 0)
    {
      printf("Applying (PI)SNII feedback...\n");
      fflush(stdout);
    }

  bg_popiii_thermal_feedback();
#endif

  reconstruct_timebins();
}


void synchronize_particle(int i)
{
  int j, no, ti_step_old, ti_step_new, ti_prev_for_bin, ti_next_for_bin;
  int bin = 0, binold, dt_bin;

  double dt_entr, dt_gravkick, dt_hydrokick;
  double time0, time1;

  MyFloat dp[3], dp_dm[3];
  MyDouble dv[3];

  /* the particle is on an active time bin, the new time step and proper kick will be computed in timestep.c */
  if(TimeBinActive[P[i].TimeBin])
    return;

  /* find the largest active bin with particles in it */
  for(j = 0; j < TIMEBINS; j++)
    if(TimeBinActive[j])
      if(j > bin)
	bin = j;

  /* actual particle time bin */
  binold = P[i].TimeBin;

  if(binold > 0)
    {
      dt_bin = (1 << binold);
      ti_prev_for_bin = (All.Ti_Current / dt_bin) * dt_bin;
      ti_next_for_bin = (All.Ti_Current / dt_bin) * dt_bin + dt_bin;
    }
  else
    {
      printf("[SYNC %2d] bin=%d (!!!!!)\n", ThisTask, bin);
      endrun(5172);
    }

  /* DEBUG */
  /*
  printf("[TASK %2d] PID=%d, binold=%d, bin=%d, ti_begstep=%d (%d - %d), ti_current=%d, ti_old_step=%d, ti_new_step=%d, Ti_Current=%d, %d\n",
	 ThisTask,
	 P[i].ID,
	 binold,
	 bin,
	 P[i].Ti_begstep,
	 ti_prev_for_bin,
	 ti_next_for_bin,
	 P[i].Ti_current,
	 (1 << binold),
	 (1 << bin),
	 All.Ti_Current,
	 (int) (log(All.Time / All.TimeBegin) / All.Timebase_interval));
  */
  /* DEBUG */

  ti_step_old = binold ? (1 << binold) : 0;
  ti_step_new = bin ? (1 << bin) : 0;

  /*
  printf("[SYNC %2d] PID=%d, ti_step_new=%d, %d\n", ThisTask, P[i].ID, ti_step_new, (1 << bin));
  */

  time0 = P[i].Ti_begstep + ti_step_old / 2;
  time1 = All.Ti_Current - ti_step_new / 2;

  if(All.ComovingIntegrationOn)
    {
      dt_entr = (time1 - time0) * All.Timebase_interval;
      dt_gravkick = get_gravkick_factor(time0, time1);
      dt_hydrokick = get_hydrokick_factor(time0, time1);
    }
  else
    dt_entr = dt_gravkick = dt_hydrokick = (time1 - time0) * All.Timebase_interval;

  /* update the entropy at half step */
  dt_entr = time1 - time0; //(All.Ti_Current - ti_step_new / 2) - (P[i].Ti_begstep + ti_step_old / 2);
  dt_entr *= All.Timebase_interval;

  /*  if(SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr < 0)
    printf("[SYNC %2d] PID=%d, entropy=%g, dtentropy=%g, dt_entr=%g\n",
	   ThisTask, P[i].ID, SphP[i].Entropy, SphP[i].e.DtEntropy, dt_entr);
  */

  SphP[i].Entropy += SphP[i].e.DtEntropy * dt_entr;

  /* update velocity at half step */
#ifdef PMGRID
  for(j = 0; j < 3; j++)
    dv[j] =
      (P[i].g.GravAccel[j] + P[i].GravPM[j]) * dt_gravkick + SphP[i].a.HydroAccel[j] * dt_hydrokick;
#else
  for(j = 0; j < 3; j++)
    dv[j] = P[i].g.GravAccel[j] * dt_gravkick + SphP[i].a.HydroAccel[j] * dt_hydrokick;
#endif

  /* update node momentum. if not done, node would miss some momentum because
     of the artificial drifting of the particle velocity in it */
  for(j = 0; j < 3; j++)
    {
      P[i].Vel[j] += dv[j];

      dp[j] = P[i].Mass * dv[j];
#ifdef SCALARFIELD
      if(P[i].Type != 0)
        dp_dm[j] = P[i].Mass * dv[j];
      else
        dp_dm[j] = 0;
#endif 
    }

  no = Father[i];

  while(no >= 0)
    {
      force_drift_node(no, All.Ti_Current);

      for(j = 0; j < 3; j++)
        {
          Extnodes[no].dp[j] += dp[j];
#ifdef SCALARFIELD
          Extnodes[no].dp_dm[j] += dp_dm[j];
#endif
        }

      /* these will be updated in the next kick */
      /*
      if(Extnodes[no].vmax < vmax)
	Extnodes[no].vmax = vmax;

      Nodes[no].u.d.bitflags |= (1 << BITFLAG_NODEHASBEENKICKED);

      Extnodes[no].Ti_lastkicked = All.Ti_Current;

      if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))
        {
          if(Extnodes[no].Flag != GlobFlag)
            {
              Extnodes[no].Flag = GlobFlag;
              DomainList[DomainNumChanged++] = no;
            }
          break;
        }
      */

      no = Nodes[no].u.d.father;
    }

  /* new bin for the particle */
  P[i].TimeBin = bin;

  /* new ti_current and ti_begstep for the particle */
  P[i].Ti_current = All.Ti_Current;
  P[i].Ti_begstep = All.Ti_Current - ti_step_new;

  /* the particle is made active... update the following */
  NumForceUpdate++;
}


void bg_ngbmass(void)
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
	if(P[i].Type == 4)
	  {
	    if(bg_ngbmass_evaluate(i, 0, &nexport, Send_count) < 0)
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
	bg_ngbmass_evaluate(j, 1, &dummy, &dummy);

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
	  StarP[P[place].StarID].NgbMassSum += NgbMassDataOut[j].NgbMassSum;
	}

      myfree(NgbMassDataOut);
      myfree(NgbMassDataResult);
      myfree(NgbMassDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
}


int bg_ngbmass_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n;

  int startnode, numngb_inbox, listindex = 0;

  double h, h2, ngb_mass;

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
		ngb_mass += P[j].Mass;
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

  if(mode == 0)
    StarP[P[target].StarID].NgbMassSum = ngb_mass;
  else
    NgbMassDataResult[target].NgbMassSum = ngb_mass;

  return 0;
}


int bg_thermal_feedback_ngb_treefind(MyDouble searchcenter[3], MyFloat hsml, int target,
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

		  FeedbackDataGet[*nexport].StarMass = StarMass;

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
#endif




#ifdef BG_SNII_THERMAL_FEEDBACK

void bg_snii_thermal_feedback(void)
{
  int i, j, ndone, ndone_flag, dummy, iter;

  int ngrp, sendTask, recvTask, nexport, nimport;

  int place, getplace;

  int imax1, imax2;

  double hubble_a, time_hubble_a, a3inv;

  double time_begstep, dt, dtime, age_of_star_in_Gyr, age_of_star_in_Gyr_endstep, dtime_in_Gyr;

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


  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     (1 + EXTRASPACEFAC) * sizeof(struct feedbackdata_in)));
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

      FeedbackDataGet =
	(struct feedbackdata_in *) mymalloc(EXTRASPACEFAC * All.BunchSize * sizeof(struct feedbackdata_in));

      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
#if defined(BG_POPIII) && defined(BG_POPIII_THERMAL_FEEDBACK)
#ifdef BG_METALSMOOTHING
	if(P[i].Type == 4 && StarP[P[i].StarID].MetallicitySmoothed >= All.POPIII_MetallicityThreshold)
#else
	if(P[i].Type == 4 && StarP[P[i].StarID].Metallicity >= All.POPIII_MetallicityThreshold)
#endif
#else
	if(P[i].Type == 4)
#endif
	  {
	    if(All.ComovingIntegrationOn)
	      time_begstep = All.TimeBegin * exp(P[i].Ti_begstep * All.Timebase_interval);
	    else
	      time_begstep = All.TimeBegin + P[i].Ti_begstep * All.Timebase_interval;

	    /* the actual time-step, in Gadget units */
	    dt = TISTEP(P[i].TimeBin) * All.Timebase_interval;

	    if(All.ComovingIntegrationOn)
	      dtime = All.Time * dt / time_hubble_a;
	    else
	      dtime = dt;

	    dtime *= All.UnitTime_in_s / All.HubbleParam;	/* convert to  seconds */
	    dtime_in_Gyr = dtime / SEC_PER_MEGAYEAR / 1000;

	    age_of_star_in_Gyr = bg_get_elapsed_time(StarP[P[i].StarID].StarBirthTime, time_begstep, 1);	/* note: this is valid for a flat universe! */
	    age_of_star_in_Gyr_endstep = age_of_star_in_Gyr + dtime_in_Gyr;

	    /* DEBUG */
	    /*
	    if(age_of_star_in_Gyr < All.SNII_WindDelay_YR * 1e-9 &&
	       age_of_star_in_Gyr_endstep >= All.SNII_WindDelay_YR * 1e-9)
	      {
		printf("[TASK %d] ID = %d, age_0 = %f, age_1 = %f, All.SNII_WindDelay_YR = %f\n",
		       ThisTask, P[i].ID, age_of_star_in_Gyr, age_of_star_in_Gyr_endstep, All.SNII_WindDelay_YR);
		fflush(stdout);
	      }
	    */
	    /* DEBUG */

	    /* now we check whether the star is eligible to spawn a wind */
	    if(age_of_star_in_Gyr <= All.SNII_WindDelay_YR * 1.0e-9 &&
	       age_of_star_in_Gyr_endstep > All.SNII_WindDelay_YR * 1.0e-9)
	      StarMass = StarP[P[i].StarID].InitialMass;	/* all SNII went off, the star can produce wind */
	    else
	      StarMass = 0;	/* this signals that the star is not allowed to produce wind */

	    /* DEBUG */
	    /*
	    if(StarMass > 0)
	      {
		printf("[TASK %d] ID = %d, performing POPII thermal feedback\n", ThisTask, P[i].ID);
		fflush(stdout);
	      }
	    */
	    /* DEBUG */

	    /* apply SNII feedback */
	    if(StarMass > 0)
	      if(bg_snii_thermal_feedback_evaluate(i, 0, &nexport, Send_count) < 0)
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

	  FeedbackDataIn[j].Pos[0] = P[place].Pos[0];
	  FeedbackDataIn[j].Pos[1] = P[place].Pos[1];
	  FeedbackDataIn[j].Pos[2] = P[place].Pos[2];

	  FeedbackDataIn[j].Hsml = PPP[place].Hsml;

	  FeedbackDataIn[j].ID = P[place].ID;

	  FeedbackDataIn[j].StarMass = FeedbackDataGet[getplace].StarMass;
	  FeedbackDataIn[j].NgbMassSum = StarP[P[place].StarID].NgbMassSum;

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
	if(FeedbackDataGet[j].StarMass > 0)
	  bg_snii_thermal_feedback_evaluate(j, 1, &dummy, &dummy);

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      myfree(FeedbackDataGet);

      iter++;
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
}


int bg_snii_thermal_feedback_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n, ngb_heat = 0;

  unsigned int id;

  int startnode, numngb_inbox, listindex = 0;

  double h, h2;

  double dx, dy, dz, r2, ascale, a3inv;
#ifdef BG_VERBOSE_LOGFILES
  int idummy = 0;
  double ddummy = 0;
#endif
  double OldEnergy, NewEnergy, OldTemperature, dmax1, dmax2;

  double mass_star, ngb_mass;

  MyDouble *pos;


  /* NEW */
  double fac3, hubble_a, step_fac;
  double csnd, entropy, dt_courant;

  int ti_step;

  if(All.ComovingIntegrationOn)
    {
      fac3 = pow(All.Time, 3 * (1 - GAMMA) / 2.0);
      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
        + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time)
#ifdef DARKENERGY
        + DarkEnergy_a(All.Time);
#else
      + All.OmegaLambda;
#endif
      hubble_a = All.Hubble * sqrt(hubble_a);
    }
  else
    fac3 = hubble_a = 1;
  /* NEW */

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
      mass_star = StarMass;
      ngb_mass = StarP[P[target].StarID].NgbMassSum;
      id = P[target].ID;
    }
  else
    {
      pos = FeedbackDataGet[target].Pos;
      h = FeedbackDataGet[target].Hsml;
      mass_star = FeedbackDataGet[target].StarMass;
      ngb_mass = FeedbackDataGet[target].NgbMassSum;
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


  h2 = h * h;

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = FeedbackDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }


  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb_inbox = bg_thermal_feedback_ngb_treefind(&pos[0], EXTRALENGTHFAC * h, target, &startnode, mode, nexport, nsend_local);

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

	      if(r2 < EXTRALENGTHFAC * EXTRALENGTHFAC * h2)
		{
		  if(All.SNII_EnergyTransferOn == 1 && mass_star > 0)
		    {
		      /* NEW */
                      /* note that some particles may not be neighbours of any active particle
			 a drifting to the current time is then required */
		      if(P[j].Ti_current != All.Ti_Current)
			drift_particle(j, All.Ti_Current);

		      /* use a trick to shorten the timestep in the next call to timestep.c
			 change the MaxSignalVel to twice the sound speed given by the
			 entropy change of heated particles */
		      entropy = All.SNII_EnergyPerUnitMass * GAMMA_MINUS1 / pow(SphP[j].d.Density * a3inv, GAMMA_MINUS1);
		      csnd = sqrt(GAMMA * entropy * pow(SphP[j].d.Density, GAMMA_MINUS1));

		      step_fac = sqrt(r2 / (4 * h2));

		      if(step_fac < 1) /* put particles within 2*h to the timestep of heated neighbours */
			step_fac = 1;  /* or propagate to other particles a timestep incresing with distance */

		      step_fac = 1 / step_fac;

		      /* DEBUG */
		      //		      printf("[TASK %2d] PID=%d, csnd=%g, step_fac=%g\n", ThisTask, P[j].ID, csnd, step_fac);
		      /* DEBUG */

		      if(2 * step_fac * csnd > SphP[j].MaxSignalVel) /* if the particle has already been touched skip it */
			{
			  SphP[j].MaxSignalVel = 2 * step_fac * csnd;

			  /* now make the particle active so it will be considered in the call to timestep.c */
			  /*
			  if(All.ComovingIntegrationOn)
			    dt_courant = 2 * All.CourantFac * All.Time * PPP[j].Hsml / (fac3 * csnd);
			  else
			    dt_courant = 2 * All.CourantFac * PPP[j].Hsml / csnd;

			  ti_step = (int) (dt_courant * hubble_a / All.Timebase_interval);
			  */
			  synchronize_particle(j);
			}
		      /* NEW */

		      /* if particle distance is less than h, perform stochastic feedback */
		      if(r2 < h2)
			{
			  if(get_random_number(P[j].ID + id) < All.SNII_EnergyFraction *
			     All.SNII_AvailableEnergyPerUnitMass / All.SNII_EnergyPerUnitMass *
			     mass_star / ngb_mass)
			    {
			      /* get old temperature for logging */
			      OldTemperature = bg_get_temperature(j);

			      /* current energy per unit mass */
			      OldEnergy = DMAX(All.MinEgySpec, SphP[j].Entropy / GAMMA_MINUS1 *
					       pow(SphP[j].d.Density * a3inv, GAMMA_MINUS1));

			      /* updated energy per unit mass */
			      NewEnergy = OldEnergy + All.SNII_EnergyPerUnitMass;

			      /* updated entropy */
			      SphP[j].Entropy *= NewEnergy / OldEnergy;

			      /* store energy */
			      //SphP[j].i.Injected_Energy += All.SNII_EnergyPerUnitMass * P[j].Mass;
			      /* DEBUG */
			      /*
			      printf("[TASK %2d] PID=%d, SNII_energy=%g, inj_energy=%g\n",
				     ThisTask, P[j].ID, All.SNII_EnergyPerUnitMass * P[j].Mass, SphP[j].i.Injected_Energy);
			      */
			      /* DEBUG */
#ifdef BG_VERBOSE_LOGFILES
			      fprintf(FdThermalFeedback, "%e %13d %13d %e %e %e %e %e %e\n",
				      All.Time, id, P[j].ID, OldEnergy * P[j].Mass,
				      NewEnergy * P[j].Mass, OldTemperature,
				      bg_get_temperature(j), mass_star, P[j].Mass);

			      fflush(FdThermalFeedback);
#endif
			      ngb_heat++;
			    }
			}
		    }
		}
	    }
	}

#ifdef BG_VERBOSE_LOGFILES
      if(ngb_heat == 0 && mass_star > 0)
	fprintf(FdThermalFeedback, "%e %13d %13d %e %e %e %e %e %e\n",
		All.Time, id, idummy, ddummy, ddummy, ddummy, ddummy, ddummy, ddummy);
#endif

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

#endif /* BG_SNII_THERMAL_FEEDBACK */







#if defined(BG_POPIII) && defined(BG_POPIII_THERMAL_FEEDBACK)

void bg_popiii_thermal_feedback(void)
{
  int i, j, ndone, ndone_flag, dummy, iter;

  int ngrp, sendTask, recvTask, nexport, nimport;

  int place, getplace;

  int imax1, imax2;

  double hubble_a, time_hubble_a, a3inv, metallicity;

  double time_begstep, dt, dtime, age_of_star_in_Gyr, age_of_star_in_Gyr_endstep, dtime_in_Gyr;

  double dying_mass, dying_mass_endstep;

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


  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     (1 + EXTRASPACEFAC) * sizeof(struct feedbackdata_in)));
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

      FeedbackDataGet =
	(struct feedbackdata_in *) mymalloc(EXTRASPACEFAC * All.BunchSize * sizeof(struct feedbackdata_in));

      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
#ifdef BG_METALSMOOTHING
	if(P[i].Type == 4 && StarP[P[i].StarID].MetallicitySmoothed < All.POPIII_MetallicityThreshold)
#else
	if(P[i].Type == 4 && StarP[P[i].StarID].Metallicity < All.POPIII_MetallicityThreshold)
#endif
	  {
	    if(All.ComovingIntegrationOn)
	      time_begstep = All.TimeBegin * exp(P[i].Ti_begstep * All.Timebase_interval);
	    else
	      time_begstep = All.TimeBegin + P[i].Ti_begstep * All.Timebase_interval;

	    /* the actual time-step, in Gadget units */
	    dt = TISTEP(P[i].TimeBin) * All.Timebase_interval;

	    if(All.ComovingIntegrationOn)
	      dtime = All.Time * dt / time_hubble_a;
	    else
	      dtime = dt;

	    dtime *= All.UnitTime_in_s / All.HubbleParam;	/* convert to  seconds */
	    dtime_in_Gyr = dtime / SEC_PER_MEGAYEAR / 1000;

	    age_of_star_in_Gyr = bg_get_elapsed_time(StarP[P[i].StarID].StarBirthTime, time_begstep, 1);	/* note: this is valid for a flat universe! */
	    age_of_star_in_Gyr_endstep = age_of_star_in_Gyr + dtime_in_Gyr;

#ifdef BG_METALSMOOTHING
	    metallicity = StarP[P[i].StarID].MetallicitySmoothed;
#else
	    metallicity = StarP[P[i].StarID].Metallicity;
#endif

	    dying_mass = dying_mass_msun(age_of_star_in_Gyr, metallicity);
	    dying_mass_endstep = dying_mass_msun(age_of_star_in_Gyr_endstep, metallicity);

	    /* DEBUG */
	    /*
	    if(dying_mass > All.POPIII_MinMass_MSUN && dying_mass_endstep <= All.POPIII_MinMass_MSUN)
	      {
		printf("[TASK %d] ID = %d, dying_mass_0 = %f, dying_mass_1 = %f, All.POPIII_MinMass_MSUN = %f\n",
		       ThisTask, P[i].ID, dying_mass, dying_mass_endstep, All.POPIII_MinMass_MSUN);
		fflush(stdout);
	      }
	    */
	    /* DEBUG */

	    /* now we check whether the star is eligible to spawn a wind */
	    if(dying_mass > All.POPIII_MinMass_MSUN && dying_mass_endstep <= All.POPIII_MinMass_MSUN)
	      StarMass = StarP[P[i].StarID].InitialMass;	/* all SNII/PISN went off, the star can produce wind */
	    else
	      StarMass = 0;	/* this signals that the star is not allowed to produce wind */

	    /* DEBUG */
	    /*
	    if(StarMass > 0)
	      {
		printf("[TASK %d] ID = %d, performing POPIII thermal feedback\n", ThisTask, P[i].ID);
		fflush(stdout);
	      }
	    */
	    /* DEBUG */

	    /* apply POPIII feedback */
	    if(bg_popiii_thermal_feedback_evaluate(i, 0, &nexport, Send_count) < 0)
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

	  FeedbackDataIn[j].Pos[0] = P[place].Pos[0];
	  FeedbackDataIn[j].Pos[1] = P[place].Pos[1];
	  FeedbackDataIn[j].Pos[2] = P[place].Pos[2];

	  FeedbackDataIn[j].Hsml = PPP[place].Hsml;

	  FeedbackDataIn[j].ID = P[place].ID;

	  FeedbackDataIn[j].StarMass = FeedbackDataGet[getplace].StarMass;
	  FeedbackDataIn[j].NgbMassSum = StarP[P[place].StarID].NgbMassSum;

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
	bg_popiii_thermal_feedback_evaluate(j, 1, &dummy, &dummy);

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      myfree(FeedbackDataGet);

      iter++;
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
}


int bg_popiii_thermal_feedback_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n, ngb_heat = 0;

  unsigned int id;

  int startnode, numngb_inbox, listindex = 0;

  double h, h2;

  double dx, dy, dz, r2, ascale, a3inv;
#ifdef BG_VERBOSE_LOGFILES
  int idummy = 0;
  double ddummy = 0;
#endif
  double OldEnergy, NewEnergy, OldTemperature, dmax1, dmax2;

  double mass_star, ngb_mass;

  MyDouble *pos;


  /* NEW */
  double fac3, hubble_a, step_fac;
  double csnd, entropy, dt_courant;

  int ti_step;

  if(All.ComovingIntegrationOn)
    {
      fac3 = pow(All.Time, 3 * (1 - GAMMA) / 2.0);
      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
        + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time)
#ifdef DARKENERGY
        + DarkEnergy_a(All.Time);
#else
      + All.OmegaLambda;
#endif
      hubble_a = All.Hubble * sqrt(hubble_a);
    }
  else
    fac3 = hubble_a = 1;
  /* NEW */


  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
      mass_star = StarMass;
      ngb_mass = StarP[P[target].StarID].NgbMassSum;
      id = P[target].ID;
    }
  else
    {
      pos = FeedbackDataGet[target].Pos;
      h = FeedbackDataGet[target].Hsml;
      mass_star = FeedbackDataGet[target].StarMass;
      ngb_mass = FeedbackDataGet[target].NgbMassSum;
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


  h2 = h * h;

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = FeedbackDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }



  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb_inbox = bg_thermal_feedback_ngb_treefind(&pos[0], EXTRALENGTHFAC * h, target, &startnode, mode, nexport, nsend_local);

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

	      if(r2 < EXTRALENGTHFAC * EXTRALENGTHFAC * h2)
		{
		  if(All.POPIII_EnergyTransferOn == 1 && mass_star > 0)
		    {
		      /* NEW */
		      /* note that some particles may not be neighbours of any active particle
			 a drifting to the current time is then required */
		      if(P[j].Ti_current != All.Ti_Current)
			drift_particle(j, All.Ti_Current);

		      /* use a trick to shorten the timestep in the next call to timestep.c
			 change the MaxSignalVel to twice the sound speed given by the
			 entropy change of heated particles */
		      entropy = All.POPIII_EnergyPerUnitMass * GAMMA_MINUS1 / pow(SphP[j].d.Density * a3inv, GAMMA_MINUS1);
		      csnd = sqrt(GAMMA * entropy * pow(SphP[j].d.Density, GAMMA_MINUS1));

		      step_fac = sqrt(r2 / h2);

		      if(step_fac < 1) /* put particles within 2*h to the timestep of heated neighbours */
			step_fac = 1;  /* or propagate to other particles a timestep increasing with distance */

		      step_fac = 1 / step_fac;

/* 		      /\* DEBUG *\/ */
/* 		      printf("[TASK %2d] PID=%d, csnd=%g, step_fac=%g\n", ThisTask, P[j].ID, csnd, step_fac); */
/* 		      /\* DEBUG *\/ */

		      if(2 * step_fac * csnd > SphP[j].MaxSignalVel) /* if the particle has already been touched skip it */
			{
			  SphP[j].MaxSignalVel = 2 * step_fac * csnd;

			  /* now make the particle active so it will be considered in the call to timestep.c */
			  /*
			  if(All.ComovingIntegrationOn)
			    dt_courant = 2 * All.CourantFac * All.Time * PPP[j].Hsml / (fac3 * csnd);
			  else
			    dt_courant = 2 * All.CourantFac * PPP[j].Hsml / csnd;

			  ti_step = (int) (dt_courant * hubble_a / All.Timebase_interval);
			  */
			  synchronize_particle(j);
			}
		      /* NEW */

		      /* if particle distance is less than h, perform stochastic feedback */
		      if(r2 < h2)
			{
			  if(get_random_number(P[j].ID + id) < All.POPIII_EnergyFraction *
			     All.POPIII_AvailableEnergyPerUnitMass / All.POPIII_EnergyPerUnitMass *
			     mass_star / ngb_mass)
			    {
			      /* get old temperature for logging */
			      OldTemperature = bg_get_temperature(j);

			      /* current energy per unit mass */
			      OldEnergy = DMAX(All.MinEgySpec, SphP[j].Entropy / GAMMA_MINUS1 *
					       pow(SphP[j].d.Density * a3inv, GAMMA_MINUS1));

			      /* updated energy per unit mass */
			      NewEnergy = OldEnergy + All.POPIII_EnergyPerUnitMass;

			      /* updated entropy */
			      SphP[j].Entropy *= NewEnergy / OldEnergy;

			      /* store energy */
			      //SphP[j].i.Injected_Energy += All.POPIII_EnergyPerUnitMass * P[j].Mass;
			      /* DEBUG */
			      /*
			      printf("[TASK %2d] PID=%d, SNII_energy=%g, inj_energy=%g\n",
				     ThisTask, P[j].ID, All.SNII_EnergyPerUnitMass * P[j].Mass, SphP[j].i.Injected_Energy);
			      */
			      /* DEBUG */
#ifdef BG_VERBOSE_LOGFILES
			      fprintf(FdThermalFeedback, "%e %13d %13d %e %e %e %e %e %e\n",
				      All.Time, id, P[j].ID, OldEnergy * P[j].Mass,
				      NewEnergy * P[j].Mass, OldTemperature,
				      bg_get_temperature(j), mass_star, P[j].Mass);

			      fflush(FdThermalFeedback);
#endif
			      ngb_heat++;
			    }
			}
		    }
		}
	    }
	}

#ifdef BG_VERBOSE_LOGFILES
      if(ngb_heat == 0 && mass_star > 0)
	fprintf(FdThermalFeedback, "%e %13d %13d %e %e %e %e %e %e\n",
		All.Time, id, idummy, ddummy, ddummy, ddummy, ddummy, ddummy, ddummy);
#endif

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
#endif /* BG_POPIII_THERMAL_FEEDBACK */










#ifdef BG_SNII_KINETIC_FEEDBACK
#ifdef BG_DOUBLE_IMF
void bg_snii_kinetic_feedback(int j, unsigned int id, MyDouble pos[3], double mass_star, double ngb_mass, double birth_density)
#else
void bg_snii_kinetic_feedback(int j, unsigned int id, MyDouble pos[3], double mass_star, double ngb_mass)
#endif
{
  int k;

  double prob, dir[3], theta, phi, norm, ascale;

#ifdef BG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING
  double time_begstep, time_in_wind_in_Gyr;
#endif


  if(All.ComovingIntegrationOn)
    ascale = All.Time;
  else
    ascale = 1;


  /* compute probability for the wind */
#ifdef BG_DOUBLE_IMF
  if(birth_density * a3inv > All.IMF_PhysDensThresh)
    prob = All.SNII_WindMassLoading1 * mass_star / ngb_mass;	/* Second IMF */
  else
    prob = All.SNII_WindMassLoading * mass_star / ngb_mass;	/* First IMF */
#else
  prob = All.SNII_WindMassLoading * mass_star / ngb_mass;
#endif


  if(All.SNII_WindOn == 1 && mass_star > 0)
    {
#ifdef BG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING
      if(SphP[j].WindFlag < 0)
	{
	  /*
	   * check whether the wind particle is still active computing
	   * the elapsed time since the particle was put in the wind
	   */
	  if(All.ComovingIntegrationOn)
	    time_begstep = All.TimeBegin * exp(P[j].Ti_begstep * All.Timebase_interval);
	  else
	    time_begstep = All.TimeBegin + P[j].Ti_begstep * All.Timebase_interval;
	  
	  /*
	   * note: this is valid for a flat universe (mode = 1)!
	   * use mode = 0 for any other cosmology
	   */
	  time_in_wind_in_Gyr =
	    bg_get_elapsed_time(-(double) SphP[j].WindFlag, time_begstep, 1);
	}
      else
	time_in_wind_in_Gyr = 0;

      /* if the particle is currently in the wind we skip it */
      if(SphP[j].WindFlag == 0
	 || time_in_wind_in_Gyr > All.SNII_WindDecouplingTime_YR * 1.0e-9)
#else
      if(SphP[j].WindFlag == 0)
#endif
	{
	  if(get_random_number(P[j].ID + id) < prob)	/* put this gas particle into the wind */
	    {
	      SphP[j].Sfr = 0;
	      
	      if(SphP[j].OnEOS == 1)
		SphP[j].OnEOS = -(MyFloat) All.Time;
		
	      ActualMassLoading += P[j].Mass;

	      if(All.SNII_WindIsotropicOn == 0)
		{
		  /* kick the particle in the star-gas particle direction */
		  for(k = 0, norm = 0; k < 3; k++)
		    {
		      dir[k] = P[j].Pos[k] - pos[k];
		      norm += dir[k] * dir[k];
		    }

		  if(norm > 0)
		    {
		      norm = sqrt(norm);
		      for(k = 0; k < 3; k++)
			dir[k] /= norm;
		    }
		  else
		    {
		      theta = acos(2 * get_random_number(P[j].ID + 1 + id) - 1);
		      phi = 2 * M_PI * get_random_number(P[j].ID + 2 + id);
		      
		      dir[0] = sin(theta) * cos(phi);
		      dir[1] = sin(theta) * sin(phi);
		      dir[2] = cos(theta);
		    }
		}
	      else
		{
		  /* kick the particle in a random direction */
		  theta = acos(2 * get_random_number(P[j].ID + 1 + id) - 1);
		  phi = 2 * M_PI * get_random_number(P[j].ID + 2 + id);

		  dir[0] = sin(theta) * cos(phi);
		  dir[1] = sin(theta) * sin(phi);
		  dir[2] = cos(theta);
		}

	      for(k = 0; k < 3; k++)
		{
		  P[j].Vel[k] += All.SNII_WindSpeed_KMpS * ascale * dir[k];
		  SphP[j].VelPred[k] += All.SNII_WindSpeed_KMpS * ascale * dir[k];
		}

	      SphP[j].WindFlag = -(MyFloat) All.Time;	/* flag the particle with the time it was put in the wind */
	    }
	}
    }
}

double bg_snii_kinetic_feedback_second_loop(void)
{
  int i, k;

  double MassSF, GlobalMassSF;

  double prob = 0, theta, phi, dir[3], ascale;

  double star_forming_mass, global_star_forming_mass;

  double global_expected_mass_loading, global_actual_mass_loading, global_actual_mass_loading_2nd_loop;

  double expected_mass_loading, actual_mass_loading, actual_mass_loading_2nd_loop;

  static double extra_mass = 0;

#ifdef BG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING
  double time_begstep, time_in_wind_in_Gyr;
#endif


  if(All.ComovingIntegrationOn)
    ascale = All.Time;
  else
    ascale = 1;


  actual_mass_loading_2nd_loop = MassSF = 0;

  if(All.SNII_WindOn == 1)
    {
      expected_mass_loading = ExpectedMassLoading;
      actual_mass_loading = ActualMassLoading;

      MPI_Allreduce(&expected_mass_loading, &global_expected_mass_loading, 1, MPI_DOUBLE, MPI_SUM,
		    MPI_COMM_WORLD);
      MPI_Allreduce(&actual_mass_loading, &global_actual_mass_loading, 1, MPI_DOUBLE, MPI_SUM,
		    MPI_COMM_WORLD);

      if(ThisTask == 0)
	printf("global_expected_mass_loading = %f, global_actual_mass_loading = %f\n",
	       global_expected_mass_loading, global_actual_mass_loading);

      extra_mass = extra_mass + global_expected_mass_loading - global_actual_mass_loading;

      if(extra_mass <= 0)
	return 0;

      for(i = 0; i < N_gas; i++)
	{
#ifdef BG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING
	  if(SphP[i].WindFlag < 0)
	    {
	      /*
	       * check whether the wind particle is still active computing
	       * the elapsed time since the particle was put in the wind
	       */
	      if(All.ComovingIntegrationOn)
		time_begstep = All.TimeBegin * exp(P[i].Ti_begstep * All.Timebase_interval);
	      else
		time_begstep = All.TimeBegin + P[i].Ti_begstep * All.Timebase_interval;

	      /*
	       * note: this is valid for a flat universe (mode = 1)!
	       * use mode = 0 for any other cosmology
	       */
	      time_in_wind_in_Gyr = bg_get_elapsed_time(-(double) SphP[i].WindFlag, time_begstep, 1);
	    }
	  else
	    time_in_wind_in_Gyr = 0;

	  /* if the particle is currently in the wind we skip it */
	  if(SphP[i].WindFlag == 0 || time_in_wind_in_Gyr > All.SNII_WindDecouplingTime_YR * 1.0e-9)
#else
	  if(SphP[i].WindFlag == 0)
#endif
	    if(SphP[i].Sfr > 0)
	      MassSF += P[i].Mass;
	}

      star_forming_mass = MassSF;

      MPI_Allreduce(&star_forming_mass, &global_star_forming_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      GlobalMassSF = global_star_forming_mass;

      if(ThisTask == 0)
	printf("GlobalMassSF = %e\n", GlobalMassSF);

      if(GlobalMassSF == 0)
	return 0;

      for(i = 0; i < N_gas; i++)
	{
#ifdef BG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING
	  if(SphP[i].WindFlag < 0)
	    {
	      /*
	       * check whether the wind particle is still active computing
	       * the elapsed time since the particle was put in the wind
	       */
	      if(All.ComovingIntegrationOn)
		time_begstep = All.TimeBegin * exp(P[i].Ti_begstep * All.Timebase_interval);
	      else
		time_begstep = All.TimeBegin + P[i].Ti_begstep * All.Timebase_interval;

	      /*
	       * note: this is valid for a flat universe (mode = 1)!
	       * use mode = 0 for any other cosmology
	       */
	      time_in_wind_in_Gyr = bg_get_elapsed_time(-(double) SphP[i].WindFlag, time_begstep, 1);
	    }
	  else
	    time_in_wind_in_Gyr = 0;

	  /* if the particle is currently in the wind we skip it */
	  if(SphP[i].WindFlag == 0 || time_in_wind_in_Gyr > All.SNII_WindDecouplingTime_YR * 1.0e-9)
#else
	  if(SphP[i].WindFlag == 0)
#endif
	    {
	      if(GlobalMassSF > 0 && SphP[i].Sfr > 0)
		prob = extra_mass / GlobalMassSF;
	      else
		prob = -1;
	    }

	  if(get_random_number(P[i].ID) < prob)	/* put this gas particle into the wind */
	    {
	      SphP[i].Sfr = 0;

	      if(SphP[i].OnEOS == 1)
		SphP[i].OnEOS = -(MyFloat) All.Time;

	      actual_mass_loading_2nd_loop += P[i].Mass;

	      /* kick the particle in a random direction */
	      theta = acos(2 * get_random_number(P[i].ID + 1) - 1);
	      phi = 2 * M_PI * get_random_number(P[i].ID + 2);

	      dir[0] = sin(theta) * cos(phi);
	      dir[1] = sin(theta) * sin(phi);
	      dir[2] = cos(theta);

	      for(k = 0; k < 3; k++)
		{
		  P[i].Vel[k] += All.SNII_WindSpeed_KMpS * ascale * dir[k];
		  SphP[i].VelPred[k] += All.SNII_WindSpeed_KMpS * ascale * dir[k];
		}

	      SphP[i].WindFlag = -(MyFloat) All.Time;	/* flag the particle with the time it was put in the wind */
	    }
	}
    }

  MPI_Allreduce(&actual_mass_loading_2nd_loop, &global_actual_mass_loading_2nd_loop, 1, MPI_DOUBLE, MPI_SUM,
		MPI_COMM_WORLD);

  extra_mass = extra_mass - global_actual_mass_loading_2nd_loop;

  return actual_mass_loading_2nd_loop;
}
#endif
