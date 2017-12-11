#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#ifdef SUBFIND

#include "allvars.h"
#include "proto.h"
#include "domain.h"
#include "fof.h"
#include "subfind.h"


static int *toGo, *toGoSph;

static int *local_toGo, *local_toGoSph;

static int *list_NumPart;

static int *list_N_gas;

#ifdef BG_SFR
static int *list_N_star, *local_toGoStar, *toGoStar;
#endif


void subfind_countToGo(int mode);


void subfind_distribute_groups(void)
{
  int i, nexport = 0, nimport = 0, target, ngrp, sendTask, recvTask;

  struct group_properties *send_Group;

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Ngroups; i++)
    {
      target = (Group[i].GrNr - 1) % NTask;
      if(target != ThisTask)
	Send_count[target]++;
    }
  MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, Recv_offset[0] = Send_offset[0] = 0; i < NTask; i++)
    {
      Recv_count[i] = Sendcount_matrix[i * NTask + ThisTask];

      nimport += Recv_count[i];
      nexport += Send_count[i];

      if(i > 0)
	{
	  Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
	  Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
	}
    }

  send_Group = (struct group_properties *) mymalloc(nexport * sizeof(struct group_properties));

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Ngroups; i++)
    {
      target = (Group[i].GrNr - 1) % NTask;
      if(target != ThisTask)
	{
	  send_Group[Send_offset[target] + Send_count[target]] = Group[i];
	  Send_count[target]++;

	  Group[i] = Group[Ngroups - 1];
	  Ngroups--;
	  i--;
	}
    }

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the group info */
	      MPI_Sendrecv(&send_Group[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct group_properties), MPI_BYTE,
			   recvTask, TAG_DENS_A,
			   &Group[Ngroups + Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct group_properties), MPI_BYTE,
			   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  Ngroups += nimport;

  myfree(send_Group);
}










void subfind_distribute_particles(int mode)
{
  int count_togo = 0, count_togo_sph = 0, count_get = 0, count_get_sph = 0;

  int *count, *count_sph, *offset, *offset_sph;

  int *count_recv, *count_recv_sph, *offset_recv, *offset_recv_sph;

#ifdef BG_SFR
  int count_togo_star = 0, count_get_star = 0, *count_star, *offset_star;

  int *count_recv_star, *offset_recv_star;
#endif
  int i, n, ngrp, target, n_requests;

  struct particle_data *partBuf;

  struct sph_particle_data *sphBuf;

#ifdef BG_SFR
  struct star_particle_data *starBuf;
#endif
  MPI_Request *requests;

  toGo = (int *) mymalloc((sizeof(int) * NTask * NTask));
  toGoSph = (int *) mymalloc((sizeof(int) * NTask * NTask));
#ifdef BG_SFR
  toGoStar = (int *) mymalloc((sizeof(int) * NTask * NTask));
#endif
  local_toGo = (int *) mymalloc((sizeof(int) * NTask));
  local_toGoSph = (int *) mymalloc((sizeof(int) * NTask));
#ifdef BG_SFR
  local_toGoStar = (int *) mymalloc((sizeof(int) * NTask));
#endif
  list_NumPart = (int *) mymalloc((sizeof(int) * NTask));
  list_N_gas = (int *) mymalloc((sizeof(int) * NTask));
#ifdef BG_SFR
  list_N_star = (int *) mymalloc((sizeof(int) * NTask));
#endif

  requests = (MPI_Request *) mymalloc(16 * NTask * sizeof(MPI_Request));

  count = (int *) mymalloc(NTask * sizeof(int));
  count_sph = (int *) mymalloc(NTask * sizeof(int));
#ifdef BG_SFR
  count_star = (int *) mymalloc(NTask * sizeof(int));
#endif
  offset = (int *) mymalloc(NTask * sizeof(int));
  offset_sph = (int *) mymalloc(NTask * sizeof(int));
#ifdef BG_SFR
  offset_star = (int *) mymalloc(NTask * sizeof(int));
#endif

  count_recv = (int *) mymalloc(NTask * sizeof(int));
  count_recv_sph = (int *) mymalloc(NTask * sizeof(int));
#ifdef BG_SFR
  count_recv_star = (int *) mymalloc(NTask * sizeof(int));
#endif
  offset_recv = (int *) mymalloc(NTask * sizeof(int));
  offset_recv_sph = (int *) mymalloc(NTask * sizeof(int));
#ifdef BG_SFR
  offset_recv_star = (int *) mymalloc(NTask * sizeof(int));
#endif


  if(ThisTask == 0)
    printf("Start counting particles to send (mode = %1d)... ", mode);

  subfind_countToGo(mode);

  if(ThisTask == 0)
    printf("done\n");


  //printf("ThisTask = %d, pos 1, Ncollective = %d\n", ThisTask, Ncollective);
  //fflush(stdout);


  for(i = 1, offset_sph[0] = 0; i < NTask; i++)
    offset_sph[i] = offset_sph[i - 1] + toGoSph[ThisTask * NTask + i - 1];

#ifdef BG_SFR
  offset_star[0] = offset_sph[NTask - 1] + toGoSph[ThisTask * NTask + NTask - 1];

  for(i = 1; i < NTask; i++)
    offset_star[i] = offset_star[i - 1] + toGoStar[ThisTask * NTask + i - 1];

  offset[0] = offset_star[NTask - 1] + toGoStar[ThisTask * NTask + NTask - 1];

  for(i = 1; i < NTask; i++)
    offset[i] =
      offset[i - 1] + (toGo[ThisTask * NTask + i - 1] - toGoSph[ThisTask * NTask + i - 1] -
		       toGoStar[ThisTask * NTask + i - 1]);
#else
  offset[0] = offset_sph[NTask - 1] + toGoSph[ThisTask * NTask + NTask - 1];

  for(i = 1; i < NTask; i++)
    offset[i] = offset[i - 1] + (toGo[ThisTask * NTask + i - 1] - toGoSph[ThisTask * NTask + i - 1]);
#endif


  //printf("ThisTask = %d, pos 2\n", ThisTask);
  //fflush(stdout);


  for(i = 0; i < NTask; i++)
    {
      count_togo += toGo[ThisTask * NTask + i];
      count_togo_sph += toGoSph[ThisTask * NTask + i];
#ifdef BG_SFR
      count_togo_star += toGoStar[ThisTask * NTask + i];
#endif

      count_get += toGo[i * NTask + ThisTask];
      count_get_sph += toGoSph[i * NTask + ThisTask];
#ifdef BG_SFR
      count_get_star += toGoStar[i * NTask + ThisTask];
#endif
    }


  //printf("ThisTask = %d, pos 3\n", ThisTask);
  //fflush(stdout);


  partBuf = (struct particle_data *) mymalloc_msg(count_togo * sizeof(struct particle_data), "partBuf");
  sphBuf =
    (struct sph_particle_data *) mymalloc_msg(count_togo_sph * sizeof(struct sph_particle_data), "sphBuf");
#ifdef BG_SFR
  starBuf =
    (struct star_particle_data *) mymalloc_msg(count_togo_star * sizeof(struct star_particle_data),
					       "starBuf");
#endif


  //printf("ThisTask = %d, pos 4\n", ThisTask);
  //fflush(stdout);


  for(i = 0; i < NTask; i++)
    {
      count[i] = count_sph[i] = 0;
#ifdef BG_SFR
      count_star[i] = 0;
#endif
    }


  //printf("ThisTask = %d, pos 5\n", ThisTask);
  //fflush(stdout);


  for(n = 0; n < NumPart; n++)
    {
      if(P[n].Type & 16)
	{
	  P[n].Type &= 15;

	  target = -1000;

	  if(mode == 0)
	    {
	      //              if(P[n].GrNr > Ncollective && P[n].GrNr <= TotNgroups)    /* particle is in small group */
	      if(!(P[n].GrNr <= TotNgroups))	/* particle is in small group */
		continue;

	      target = (P[n].GrNr - 1) % NTask;
	    }
	  else if(mode == 1)
	    target = P[n].targettask;
	  else if(mode == 2)
	    target = P[n].origintask;


	  if(target < 0)
	    printf("ThisTask = %d, target = %d, P[n].GrNr = %d, mode = %d, TotNgroups = %d\n", ThisTask,
		   target, P[n].GrNr, mode, TotNgroups);

/* 	  printf("ThisTask = %d, target = %d, part n = %d\n", ThisTask, target, n); */
/* 	  fflush(stdout); */

	  if(P[n].Type == 0)
	    {
	      partBuf[offset_sph[target] + count_sph[target]] = P[n];
	      sphBuf[offset_sph[target] + count_sph[target]] = SphP[n];
	      count_sph[target]++;
	    }
#ifdef BG_SFR
	  else if(P[n].Type == 4)
	    {
	      partBuf[offset_star[target] + count_star[target]] = P[n];
	      starBuf[offset_star[target] - offset_star[0] + count_star[target]] = StarP[P[n].StarID];
	      count_star[target]++;
	    }
#endif
	  else
	    {
	      partBuf[offset[target] + count[target]] = P[n];
	      count[target]++;
	    }


	  if(P[n].Type == 0)
	    {
	      P[n] = P[N_gas - 1];
	      P[N_gas - 1] = P[NumPart - 1];
#ifdef BG_SFR
	      if((P[n].Type & 15) == 4)
		StarP[P[n].StarID].PID = n;
	      if((P[N_gas - 1].Type & 15) == 4)
		StarP[P[N_gas - 1].StarID].PID = N_gas - 1;
#endif
	      SphP[n] = SphP[N_gas - 1];

	      NumPart--;
	      N_gas--;
	      n--;
	    }
#ifdef BG_SFR
	  else if(P[n].Type == 4)
	    {
	      StarP[P[n].StarID] = StarP[N_star - 1];
	      P[StarP[N_star - 1].PID].StarID = P[n].StarID;

	      P[n] = P[NumPart - 1];

	      if((P[n].Type & 15) == 4 && n != (NumPart - 1))
		StarP[P[n].StarID].PID = n;

	      NumPart--;
	      N_star--;
	      n--;
	    }
#endif
	  else
	    {
	      P[n] = P[NumPart - 1];
#ifdef BG_SFR
	      if((P[n].Type & 15) == 4)
		StarP[P[n].StarID].PID = n;
#endif
	      NumPart--;
	      n--;
	    }
	}
    }


  //printf("ThisTask = %d, pos 6\n", ThisTask);
  //fflush(stdout);


  if(count_get_sph)
    {
      memmove(P + N_gas + count_get_sph, P + N_gas, (NumPart - N_gas) * sizeof(struct particle_data));

#ifdef BG_SFR
      for(i = 0; i < N_star; i++)
	if(StarP[i].PID >= (unsigned int) N_gas)
	  StarP[i].PID += count_get_sph;
#endif
    }


  //printf("ThisTask = %d, pos 7\n", ThisTask);
  //fflush(stdout);


  for(i = 0; i < NTask; i++)
    {
      count_recv_sph[i] = toGoSph[i * NTask + ThisTask];
#ifndef BG_SFR
      count_recv[i] = toGo[i * NTask + ThisTask] - toGoSph[i * NTask + ThisTask];
#else
      count_recv_star[i] = toGoStar[i * NTask + ThisTask];
      count_recv[i] =
	toGo[i * NTask + ThisTask] - toGoSph[i * NTask + ThisTask] - toGoStar[i * NTask + ThisTask];
#endif
    }


  //printf("ThisTask = %d, pos 8\n", ThisTask);
  //fflush(stdout);


  for(i = 1, offset_recv_sph[0] = N_gas; i < NTask; i++)
    offset_recv_sph[i] = offset_recv_sph[i - 1] + count_recv_sph[i - 1];


#ifdef BG_SFR
  for(i = 1, offset_recv_star[0] = N_star; i < NTask; i++)
    offset_recv_star[i] = offset_recv_star[i - 1] + count_recv_star[i - 1];
#endif

#ifndef BG_SFR
  offset_recv[0] = NumPart + count_get_sph;
#else
  offset_recv[0] = NumPart + count_get_sph + count_get_star;
#endif


  //printf("ThisTask = %d, pos 9\n", ThisTask);
  //fflush(stdout);


  for(i = 1; i < NTask; i++)
    offset_recv[i] = offset_recv[i - 1] + count_recv[i - 1];


  //printf("ThisTask = %d, pos 10\n", ThisTask);
  //fflush(stdout);


  n_requests = 0;

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
	{
	  if(count_recv_sph[target] > 0)
	    {
	      MPI_Irecv(P + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct particle_data),
			MPI_BYTE, target, TAG_PDATA_SPH, MPI_COMM_WORLD, &requests[n_requests++]);

	      MPI_Irecv(SphP + offset_recv_sph[target],
			count_recv_sph[target] * sizeof(struct sph_particle_data), MPI_BYTE, target,
			TAG_SPHDATA, MPI_COMM_WORLD, &requests[n_requests++]);
	    }

#ifdef BG_SFR
	  if(count_recv_star[target] > 0)
	    {
	      MPI_Irecv(StarP + offset_recv_star[target],
			count_recv_star[target] * sizeof(struct star_particle_data), MPI_BYTE, target,
			TAG_STARDATA, MPI_COMM_WORLD, &requests[n_requests++]);

	      MPI_Irecv(P + offset_recv_star[target] + NumPart + count_get_sph - N_star,
			count_recv_star[target] * sizeof(struct particle_data), MPI_BYTE, target,
			TAG_PDATA_STAR, MPI_COMM_WORLD, &requests[n_requests++]);
	    }
#endif

	  if(count_recv[target] > 0)
	    {
	      MPI_Irecv(P + offset_recv[target], count_recv[target] * sizeof(struct particle_data),
			MPI_BYTE, target, TAG_PDATA, MPI_COMM_WORLD, &requests[n_requests++]);
	    }
	}
    }


  //printf("ThisTask = %d, pos 11\n", ThisTask);
  //fflush(stdout);


  MPI_Barrier(MPI_COMM_WORLD);	/* not really necessary, but this will guarantee that all receives are
				   posted before the sends, which helps the stability of MPI on 
				   bluegene, and perhaps some mpich1-clusters */

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
	{
	  if(count_sph[target] > 0)
	    {
	      MPI_Isend(partBuf + offset_sph[target], count_sph[target] * sizeof(struct particle_data),
			MPI_BYTE, target, TAG_PDATA_SPH, MPI_COMM_WORLD, &requests[n_requests++]);

	      MPI_Isend(sphBuf + offset_sph[target], count_sph[target] * sizeof(struct sph_particle_data),
			MPI_BYTE, target, TAG_SPHDATA, MPI_COMM_WORLD, &requests[n_requests++]);
	    }

#ifdef BG_SFR
	  if(count_star[target] > 0)
	    {
	      MPI_Isend(partBuf + offset_star[target], count_star[target] * sizeof(struct particle_data),
			MPI_BYTE, target, TAG_PDATA_STAR, MPI_COMM_WORLD, &requests[n_requests++]);

	      MPI_Isend(starBuf + offset_star[target] - offset_star[0],
			count_star[target] * sizeof(struct star_particle_data), MPI_BYTE, target,
			TAG_STARDATA, MPI_COMM_WORLD, &requests[n_requests++]);
	    }
#endif

	  if(count[target] > 0)
	    {
	      MPI_Isend(partBuf + offset[target], count[target] * sizeof(struct particle_data),
			MPI_BYTE, target, TAG_PDATA, MPI_COMM_WORLD, &requests[n_requests++]);
	    }
	}
    }

  MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);

  MPI_Barrier(MPI_COMM_WORLD);

  //printf("ThisTask = %d, pos 12\n", ThisTask);
  //fflush(stdout);


#ifdef BG_SFR
  for(i = N_star; i < N_star + count_get_star; i++)
    {
      StarP[i].PID = NumPart + count_get_sph - N_star + i;
      P[NumPart + count_get_sph - N_star + i].StarID = i;
    }
#endif

  NumPart += count_get;
  N_gas += count_get_sph;

#ifdef BG_SFR
  N_star += count_get_star;
#endif

  if(NumPart > All.MaxPart)
    {
      printf("Task=%d NumPart=%d All.MaxPart=%d\n", ThisTask, NumPart, All.MaxPart);
      endrun(787878);
    }

  if(N_gas > All.MaxPartSph)
    endrun(787879);

#ifdef BG_SFR
  if(N_star > All.MaxPartStar)
    endrun(787880);
#endif

/*   myfree(keyBuf); */
#ifdef BG_SFR
  myfree(starBuf);
#endif
  myfree(sphBuf);
  myfree(partBuf);

#ifdef BG_SFR
  myfree(offset_recv_star);
#endif
  myfree(offset_recv_sph);
  myfree(offset_recv);
#ifdef BG_SFR
  myfree(count_recv_star);
#endif
  myfree(count_recv_sph);
  myfree(count_recv);

#ifdef BG_SFR
  myfree(offset_star);
#endif
  myfree(offset_sph);
  myfree(offset);
#ifdef BG_SFR
  myfree(count_star);
#endif
  myfree(count_sph);
  myfree(count);

  myfree(requests);

#ifdef BG_SFR
  myfree(list_N_star);
#endif
  myfree(list_N_gas);
  myfree(list_NumPart);
#ifdef BG_SFR
  myfree(local_toGoStar);
#endif
  myfree(local_toGoSph);
  myfree(local_toGo);
#ifdef BG_SFR
  myfree(toGoStar);
#endif
  myfree(toGoSph);
  myfree(toGo);
}



void subfind_countToGo(int mode)
{
  int n, target;

  int ta, i;

  int count_togo, count_toget, count_togo_sph, count_toget_sph;

#ifdef BG_SFR
  int count_togo_star, count_toget_star;
#endif




  for(n = 0; n < NTask; n++)
    {
      local_toGo[n] = 0;
      local_toGoSph[n] = 0;
#ifdef BG_SFR
      local_toGoStar[n] = 0;
#endif
    }

  for(n = 0; n < NumPart; n++)
    {
      if(mode == 0)
	{
	  //      if(P[n].GrNr > Ncollective && P[n].GrNr <= TotNgroups)        /* particle is in small group */
	  if(!(P[n].GrNr <= TotNgroups))	/* particle is in small group */
	    continue;

	  target = (P[n].GrNr - 1) % NTask;
	}
      else if(mode == 1)
	target = P[n].targettask;
      else if(mode == 2)
	target = P[n].origintask;


      if(target != ThisTask)
	{
	  local_toGo[target] += 1;

	  if(P[n].Type == 0)
	    {
	      local_toGoSph[target] += 1;
	    }
#ifdef BG_SFR
	  if(P[n].Type == 4)
	    {
	      local_toGoStar[target] += 1;
	    }
#endif
	  P[n].Type |= 16;	/* flag this particle for export */
	}
    }

  MPI_Allgather(local_toGo, NTask, MPI_INT, toGo, NTask, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(local_toGoSph, NTask, MPI_INT, toGoSph, NTask, MPI_INT, MPI_COMM_WORLD);
#ifdef BG_SFR
  MPI_Allgather(local_toGoStar, NTask, MPI_INT, toGoStar, NTask, MPI_INT, MPI_COMM_WORLD);
#endif


  /* in this case, we are not guaranteed that the temporary state after
     the partial exchange will actually observe the particle limits on all
     processors... we need to test this explicitly and rework the exchange
     such that this is guaranteed. This is actually a rather non-trivial
     constraint. */

  MPI_Allgather(&NumPart, 1, MPI_INT, list_NumPart, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&N_gas, 1, MPI_INT, list_N_gas, 1, MPI_INT, MPI_COMM_WORLD);
#ifdef BG_SFR
  MPI_Allgather(&N_star, 1, MPI_INT, list_N_star, 1, MPI_INT, MPI_COMM_WORLD);
#endif


  for(ta = 0; ta < NTask; ta++)
    {
      count_togo = count_toget = 0;
      count_togo_sph = count_toget_sph = 0;
#ifdef BG_SFR
      count_togo_star = count_toget_star = 0;
#endif

      for(i = 0; i < NTask; i++)
	{
	  count_togo += toGo[ta * NTask + i];
	  count_toget += toGo[i * NTask + ta];
	  count_togo_sph += toGoSph[ta * NTask + i];
	  count_toget_sph += toGoSph[i * NTask + ta];
#ifdef BG_SFR
	  count_togo_star += toGoStar[ta * NTask + i];
	  count_toget_star += toGoStar[i * NTask + ta];
#endif
	}
      for(n = 0; n < NTask; n++)
	{
	  local_toGo[n] = 0;
	  local_toGoSph[n] = 0;
#ifdef BG_SFR
	  local_toGoStar[n] = 0;
#endif
	}

      for(n = 0; n < NumPart; n++)
	{
	  if(P[n].Type & 16)
	    {
	      P[n].Type &= 15;


	      if(mode == 0)
		{
		  //              if(P[n].GrNr > Ncollective && P[n].GrNr <= TotNgroups)        /* particle is in small group */
		  if(!(P[n].GrNr <= TotNgroups))	/* particle is in small group */
		    continue;

		  target = (P[n].GrNr - 1) % NTask;
		}
	      else if(mode == 1)
		target = P[n].targettask;
	      else if(mode == 2)
		target = P[n].origintask;


	      if(P[n].Type == 0)
		{
		  if(local_toGoSph[target] < toGoSph[ThisTask * NTask + target] &&
		     local_toGo[target] < toGo[ThisTask * NTask + target])
		    {
		      local_toGo[target] += 1;
		      local_toGoSph[target] += 1;
		      P[n].Type |= 16;
		    }
		}
#ifdef BG_SFR
	      else if(P[n].Type == 4)
		{
		  if(local_toGoStar[target] < toGoStar[ThisTask * NTask + target] &&
		     local_toGo[target] < toGo[ThisTask * NTask + target])
		    {
		      local_toGo[target] += 1;
		      local_toGoStar[target] += 1;
		      P[n].Type |= 16;
		    }
		}
#endif
	      else
		{
		  if(local_toGo[target] < toGo[ThisTask * NTask + target])
		    {
		      local_toGo[target] += 1;
		      P[n].Type |= 16;
		    }
		}
	    }
	}

      MPI_Allgather(local_toGo, NTask, MPI_INT, toGo, NTask, MPI_INT, MPI_COMM_WORLD);
      MPI_Allgather(local_toGoSph, NTask, MPI_INT, toGoSph, NTask, MPI_INT, MPI_COMM_WORLD);
#ifdef BG_SFR
      MPI_Allgather(local_toGoStar, NTask, MPI_INT, toGoStar, NTask, MPI_INT, MPI_COMM_WORLD);
#endif
    }
}

#endif
