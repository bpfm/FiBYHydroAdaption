#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <sys/stat.h>
#include <sys/types.h>


#include "allvars.h"
#include "proto.h"


#ifdef SUBFIND

#include "fof.h"
#include "subfind.h"


/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct SOdens_in
{
  MyDouble Pos[3];
  MyFloat R200;
  int NodeList[NODELISTLENGTH];
}
 *SOdensIn, *SOdensGet;


static struct SOdens_out
{
  double Mass;
}
 *SOdensResult, *SOdensOut;


static double *R200, *M200;

void subfind_overdensity(void)
{
  long long ntot;

  int i, j, ndone, ndone_flag, npleft, dummy, rep, iter;

  MyFloat *Left, *Right;

  char *Todo;

  int ngrp, sendTask, recvTask, place, nexport, nimport;

  double t0, t1, rguess, overdensity, Deltas[3], rhoback, z, omegaz, x, DeltaMean200, DeltaCrit200,
    DeltaTopHat;


  /* allocate buffers to arrange communication */


  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct SOdens_in) + sizeof(struct SOdens_out) +
					     sizemax(sizeof(struct SOdens_in), sizeof(struct SOdens_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));


  Left = (MyFloat *) mymalloc(sizeof(MyFloat) * Ngroups);
  Right = (MyFloat *) mymalloc(sizeof(MyFloat) * Ngroups);
  R200 = (double *) mymalloc(sizeof(double) * Ngroups);
  M200 = (double *) mymalloc(sizeof(double) * Ngroups);

  Todo = mymalloc(sizeof(char) * Ngroups);

  z = 1 / All.Time - 1;

  rhoback = 3 * All.Omega0 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  omegaz =
    All.Omega0 * pow(1 + z,
		     3) / (All.Omega0 * pow(1 + z, 3) + (1 - All.Omega0 - All.OmegaLambda) * pow(1 + z,
												 2) +
			   All.OmegaLambda);

  DeltaMean200 = 200.0;
  DeltaCrit200 = 200.0 / omegaz;

  x = omegaz - 1;
  DeltaTopHat = 18 * M_PI * M_PI + 82 * x - 39 * x * x;
  DeltaTopHat /= omegaz;

  Deltas[0] = DeltaMean200;	/* standard fixed overdensity with respect to background */
  Deltas[1] = DeltaTopHat;	/* tophat overdensity with respect to background */
  Deltas[2] = DeltaCrit200;	/* overdensity of 200 relative to critical, expressed relative to background density */


  for(rep = 0; rep < 3; rep++)	/* repeat for all three overdensity values */
    {
      for(i = 0; i < Ngroups; i++)
	{
	  if(Group[i].Nsubs > 0)
	    {
	      rguess = pow(All.G * Group[i].Mass / (100 * All.Hubble * All.Hubble), 1.0 / 3);
	      //All.G is G_N * U_m * (U_T)^2 / (U_L)^3
	      //All.Hubble is H * U_T 
	      //rguess is [(G_N / H^2) (U_m)^2 / (U_L)^3]^1/3

	      // Mass/h is normal units. M = int_mass*100/H
	      // 
	      Right[i] = 3 * rguess;
	      Left[i] = 0;

	      Todo[i] = 1;
	    }
	  else
	    {
	      Todo[i] = 0;
	    }
	}

      iter = 0;

      /* we will repeat the whole thing for those groups where we didn't converge to a SO radius yet */
      do
	{
	  t0 = second();

	  i = 0;		/* begin with this index */

	  do
	    {
	      for(j = 0; j < NTask; j++)
		{
		  Send_count[j] = 0;
		  Exportflag[j] = -1;
		}

	      /* do local particles and prepare export list */

	      for(nexport = 0; i < Ngroups; i++)
		{
		  if(Todo[i])
		    {
		      R200[i] = 0.5 * (Left[i] + Right[i]);
		      if(subfind_overdensity_evaluate(i, 0, &nexport, Send_count) < 0)
			break;
		    }
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

	      SOdensGet = (struct SOdens_in *) mymalloc(nimport * sizeof(struct SOdens_in));
	      SOdensIn = (struct SOdens_in *) mymalloc(nexport * sizeof(struct SOdens_in));

	      /* prepare particle data for export */
	      for(j = 0; j < nexport; j++)
		{
		  place = DataIndexTable[j].Index;

		  SOdensIn[j].Pos[0] = Group[place].Pos[0];
		  SOdensIn[j].Pos[1] = Group[place].Pos[1];
		  SOdensIn[j].Pos[2] = Group[place].Pos[2];
		  SOdensIn[j].R200 = R200[place];

		  memcpy(SOdensIn[j].NodeList,
			 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
		}

	      /* exchange data */
	      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
		{
		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
			{
			  /* get the data */
			  MPI_Sendrecv(&SOdensIn[Send_offset[recvTask]],
				       Send_count[recvTask] * sizeof(struct SOdens_in), MPI_BYTE,
				       recvTask, TAG_DENS_A,
				       &SOdensGet[Recv_offset[recvTask]],
				       Recv_count[recvTask] * sizeof(struct SOdens_in), MPI_BYTE,
				       recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		    }
		}

	      myfree(SOdensIn);
	      SOdensResult = (struct SOdens_out *) mymalloc(nimport * sizeof(struct SOdens_out));
	      SOdensOut = (struct SOdens_out *) mymalloc(nexport * sizeof(struct SOdens_out));


	      /* now do the locations that were sent to us */
	      for(j = 0; j < nimport; j++)
		subfind_overdensity_evaluate(j, 1, &dummy, &dummy);

	      if(i >= Ngroups)
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
			  MPI_Sendrecv(&SOdensResult[Recv_offset[recvTask]],
				       Recv_count[recvTask] * sizeof(struct SOdens_out),
				       MPI_BYTE, recvTask, TAG_DENS_B,
				       &SOdensOut[Send_offset[recvTask]],
				       Send_count[recvTask] * sizeof(struct SOdens_out),
				       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		    }
		}

	      /* add the result to the local particles */
	      for(j = 0; j < nexport; j++)
		{
		  place = DataIndexTable[j].Index;

		  M200[place] += SOdensOut[j].Mass;
		}

	      myfree(SOdensOut);
	      myfree(SOdensResult);
	      myfree(SOdensGet);
	    }
	  while(ndone < NTask);


	  /* do final operations on results */
	  for(i = 0, npleft = 0; i < Ngroups; i++)
	    {
	      if(Todo[i])
		{
		  overdensity = M200[i] / (4.0 * M_PI / 3.0 * R200[i] * R200[i] * R200[i]) / rhoback;

		  if((Right[i] - Left[i]) > 1.0e-4 * Left[i])
		    {
		      /* need to redo this group */
		      npleft++;

		      if(overdensity > Deltas[rep])
			Left[i] = R200[i];
		      else
			Right[i] = R200[i];

		      if(iter >= MAXITER - 10)
			{
			  printf
			    ("gr=%d task=%d  R200=%g Left=%g Right=%g Menclosed=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			     i, ThisTask, R200[i], Left[i], Right[i],
			     M200[i], Right[i] - Left[i], Group[i].Pos[0], Group[i].Pos[1], Group[i].Pos[2]);
			  fflush(stdout);
			}
		    }
		  else
		    Todo[i] = 0;
		}
	    }

	  sumup_large_ints(1, &npleft, &ntot);

	  t1 = second();

	  if(ntot > 0)
	    {
	      iter++;

	      if(iter > 0 && ThisTask == 0)
		{
		  printf("SO iteration %d: need to repeat for %d%09d particles. (took %g sec)\n", iter,
			 (int) (ntot / 1000000000), (int) (ntot % 1000000000), timediff(t0, t1));
		  fflush(stdout);
		}

	      if(iter > MAXITER)
		{
		  printf("failed to converge in neighbour iteration in subfind_overdensity()\n");
		  fflush(stdout);
		  endrun(1155);
		}
	    }
	}
      while(ntot > 0);


      for(i = 0; i < Ngroups; i++)
	{
	  if(Group[i].Nsubs > 0)
	    {
	      overdensity = M200[i] / (4.0 * M_PI / 3.0 * R200[i] * R200[i] * R200[i]) / rhoback;

	      if((overdensity - Deltas[rep]) > 0.1 * Deltas[rep])
		{
		  R200[i] = M200[i] = 0;
		}
	      else if(M200[i] < 5 * Group[i].Mass / Group[i].Len)
		{
		  R200[i] = M200[i] = 0;
		}
	    }
	  else
	    R200[i] = M200[i] = 0;

	  switch (rep)
	    {
	    case 0:
	      Group[i].M_Mean200 = M200[i];
	      Group[i].R_Mean200 = R200[i];
	      break;
	    case 1:
	      Group[i].M_TopHat200 = M200[i];
	      Group[i].R_TopHat200 = R200[i];
	      break;
	    case 2:
	      Group[i].M_Crit200 = M200[i];
	      Group[i].R_Crit200 = R200[i];
	      break;
	    }
	}
    }

  myfree(Todo);
  myfree(M200);
  myfree(R200);
  myfree(Right);
  myfree(Left);

  myfree(DataNodeList);
  myfree(DataIndexTable);
}


/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int subfind_overdensity_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int startnode, listindex = 0;

  double h, mass, massret;

  MyDouble *pos;


  if(mode == 0)
    {
      pos = Group[target].Pos;
      h = R200[target];
    }
  else
    {
      pos = SOdensGet[target].Pos;
      h = SOdensGet[target].R200;
    }

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = SOdensGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  mass = 0;

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  massret = subfind_ovderdens_treefind(pos, h, target, &startnode, mode, nexport, nsend_local);

	  if(massret < 0)
	    return -1;

	  mass += massret;	//This is a simple summation of P[i].Mass within a given radius 
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = SOdensGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  if(mode == 0)
    M200[target] = mass;
  else
    SOdensResult[target].Mass = mass;

  return 0;
}


double subfind_ovderdens_treefind(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
				  int mode, int *nexport, int *nsend_local)
{
  int no, p, task, nexport_save;

  struct NODE *current;

  double mass;

  MyDouble dx, dy, dz, dist, r2;

#define FACT2 0.86602540
#ifdef PERIODIC
  MyDouble xtmp;
#endif
  nexport_save = *nexport;

  mass = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

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

#ifdef INCL_DMONLY
	  /* Skip all non-DM particles */
	  if((1 << P[p].Type) & (FOF_SECONDARY_LINK_TYPES))
	    continue;
#endif
	  mass += P[p].Mass;
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(mode == 0)
		{
		  if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
		    {
		      Exportflag[task] = target;
		      Exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(Exportnodecount[task] == NODELISTLENGTH)
		    {
		      if(*nexport >= All.BunchSize)
			{
			  *nexport = nexport_save;
			  if(nexport_save == 0)
			    endrun(13005);	/* in this case, the buffer is too small to process even a single particle */
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
		}

	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
		  return mass;
		}
	    }

	  no = current->u.d.sibling;	/* in case the node can be discarded */
	  dist = hsml + 0.5 * current->len;
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
	  if((r2 = (dx * dx + dy * dy + dz * dz)) > dist * dist)
	    continue;

	  if((current->u.d.bitflags & ((1 << BITFLAG_TOPLEVEL) + (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))) == 0)	/* only use fully local nodes */
	    {
	      /* test whether the node is contained within the sphere */
	      dist = hsml - FACT2 * current->len;
	      if(dist > 0)
		if(r2 < dist * dist)
		  {
		    mass += current->u.d.mass;
		    continue;
		  }
	    }

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return mass;
}


#endif
