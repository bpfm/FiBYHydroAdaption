#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_math.h>
#include <inttypes.h>

#include "allvars.h"
#include "proto.h"
#include "domain.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

/*! \file fof.c
 *  \brief parallel FoF group finder
 */

#ifdef FOF


#include "fof.h"
#ifdef BG_SFR
#include "bg_vars.h"
#include "bg_proto.h"
#endif

int Ngroups, TotNgroups;
long long TotNids;

struct group_properties *Group;

static struct fofdata_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int MinID;
  int MinIDTask;
  int NodeList[NODELISTLENGTH];
}
 *FoFDataIn, *FoFDataGet;

#ifdef FOF_SECONDARY_LINK_TYPES
static struct fofdata_out
{
  MyFloat Distance;
  int MinID;
  int MinIDTask;
}
 *FoFDataResult, *FoFDataOut;
#endif

static struct fof_particle_list
{
  int MinID;
  int MinIDTask;
  int Pindex;
}
 *FOF_PList;

static struct fof_group_list
{
  int MinID;
  int MinIDTask;
  int LocCount;
  int ExtCount;
  int GrNr;
}
 *FOF_GList;

static struct id_list
{
  MyIDType ID;
  int Type;
  int GrNr;
}
 *ID_list;

static double LinkL;
static int NgroupsExt, Nids;
static int *Head, *Len, *Next, *Tail, *MinID, *MinIDTask;
static char *NonlocalFlag;

#ifdef FOF_SECONDARY_LINK_TYPES
static float *fof_nearest_distance;
static float *fof_nearest_hsml;
#endif


void fof_write_hdf5_star_field(double *data, char *tag, int elems, char *description, int blocknr, hid_t grp);


void fof_fof(int num)
{
  int i, ndm, ndmtot, start, lenloc, largestgroup, imax1, imax2;
  double mass, masstot, rhodm, t0, t1;


  if(ThisTask == 0)
    {
      printf("\nBegin to compute FoF group catalogues...  (presently allocated=%g MB)\n",
	     AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  CPU_Step[CPU_MISC] += measure_time();

  All.NumForcesSinceLastDomainDecomp = (long long) (All.TotNumPart * All.TreeDomainUpdateFrequency + 1);

  domain_Decomposition();

  force_treefree();

  for(i = 0, ndm = 0, mass = 0; i < NumPart; i++)
    if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
      {
	ndm++;
	mass += P[i].Mass;
      }

  MPI_Allreduce(&ndm, &ndmtot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  rhodm = (All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  LinkL = LINKLENGTH * pow(masstot / ndmtot / rhodm, 1.0 / 3);

  if(ThisTask == 0)
    {
      printf("\nComoving linking length: %g    ", LinkL);
      printf("(presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  FOF_PList = (struct fof_particle_list *) mymalloc(NumPart * sizeof(struct fof_particle_list));

  if(sizeof(struct fof_particle_list) != 3 * sizeof(int))
    endrun(1232);

  MinID = (int *) FOF_PList;
  MinIDTask = MinID + NumPart;
  Head = MinIDTask + NumPart;
  Len = (int *) mymalloc(NumPart * sizeof(int));
  Next = (int *) mymalloc(NumPart * sizeof(int));
  Tail = (int *) mymalloc(NumPart * sizeof(int));


  if(ThisTask == 0)
    printf("Tree construction.\n");

  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);
  force_treebuild(NumPart, NULL);


  for(i = 0; i < NumPart; i++)
    {
      Head[i] = Tail[i] = i;
      Len[i] = 1;
      Next[i] = -1;
      MinID[i] = P[i].ID;
      MinIDTask[i] = ThisTask;
    }


  t0 = second();

  fof_find_groups();

  t1 = second();
  if(ThisTask == 0)
    printf("group finding took = %g sec\n", timediff(t0, t1));

#ifdef FOF_SECONDARY_LINK_TYPES
  t0 = second();

  fof_find_nearest_dmparticle();

  t1 = second();
  if(ThisTask == 0)
    printf("attaching gas and star particles to nearest dm particles took = %g sec\n", timediff(t0, t1));
#endif

  t0 = second();

  for(i = 0; i < NumPart; i++)
    {
      Next[i] = MinID[Head[i]];
      Tail[i] = MinIDTask[Head[i]];
    }

  for(i = 0; i < NumPart; i++)
    {
      FOF_PList[i].MinID = Next[i];
      FOF_PList[i].MinIDTask = Tail[i];
      FOF_PList[i].Pindex = i;
    }

  force_treefree();

  myfree(Tail);
  myfree(Next);
  myfree(Len);

  FOF_GList = (struct fof_group_list *) mymalloc(sizeof(struct fof_group_list) * NumPart);

  fof_compile_catalogue();

  t1 = second();
  if(ThisTask == 0)
    printf("compiling local group data and catalogue took = %g sec\n", timediff(t0, t1));


  MPI_Allreduce(&Ngroups, &TotNgroups, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  sumup_large_ints(1, &Nids, &TotNids);

  if(TotNgroups > 0)
    {
      int largestloc = 0;

      for(i = 0; i < NgroupsExt; i++)
	if(FOF_GList[i].LocCount + FOF_GList[i].ExtCount > largestloc)
	  largestloc = FOF_GList[i].LocCount + FOF_GList[i].ExtCount;
      MPI_Allreduce(&largestloc, &largestgroup, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
  else
    largestgroup = 0;

  if(ThisTask == 0)
    {
      printf("\nTotal number of groups with at least %d particles: %d\n", GROUP_MIN_LEN, TotNgroups);
      if(TotNgroups > 0)
	{
	  printf("Largest group has %d particles.\n", largestgroup);
	  printf("Total number of particles in groups: %d%09d\n\n",
		 (int) (TotNids / 1000000000), (int) (TotNids % 1000000000));
	}
    }

  t0 = second();

  Group =
    (struct group_properties *) mymalloc(sizeof(struct group_properties) *
					 IMAX(NgroupsExt, TotNgroups / NTask + 1));

  if(ThisTask == 0)
    {
      printf("group properties are now allocated.. (presently allocated=%g MB)\n",
	     AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  for(i = 0, start = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
	{
	  start++;
	  if(start > NumPart)
	    endrun(78);
	}

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
	endrun(123);

      for(lenloc = 0; start + lenloc < NumPart;)
	if(FOF_PList[start + lenloc].MinID == FOF_GList[i].MinID)
	  lenloc++;
	else
	  break;

      Group[i].MinID = FOF_GList[i].MinID;
      Group[i].MinIDTask = FOF_GList[i].MinIDTask;

      fof_compute_group_properties(i, start, lenloc);

      start += lenloc;
    }

  fof_exchange_group_data();

  fof_finish_group_properties();

  t1 = second();
  if(ThisTask == 0)
    printf("computation of group properties took = %g sec\n", timediff(t0, t1));

#ifdef BLACK_HOLES
    fof_make_black_holes();
#endif

#ifdef FOF_OUTPUT
  fof_save_groups(num);
#endif

  MPI_Barrier(MPI_COMM_WORLD);

  myfree(Group);
  myfree(FOF_GList);
  myfree(FOF_PList);

  if(ThisTask == 0)
    {
      printf("Finished computing FoF groups.  (presently allocated=%g MB)\n\n",
	     AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

  if(ThisTask == 0)
    printf("Tree construction.\n");
  force_treebuild(NumPart, NULL);

  TreeReconstructFlag = 0;

  CPU_Step[CPU_FOF] += measure_time();

#ifdef USE_HDF5_FIX
  H5close();
  hdf5_memory_cleanup();
#endif
}


void fof_find_groups(void)
{
  int i, j, ndone_flag, link_count, dummy, nprocessed, ntot;
  int ndone, ngrp, sendTask, recvTask, place, nexport, nimport, link_across, link_across_tot;
  int npart, marked, totnpart, totmarked, *MinIDOld;
  char *FoFDataOut, *FoFDataResult, *MarkedFlag, *ChangedFlag;
  double t0, t1;

  if(ThisTask == 0)
    {
      printf("\nStart linking particles (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }


  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     2 * sizeof(struct fofdata_in)));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));

  NonlocalFlag = (char *) mymalloc(NumPart * sizeof(char));
  MarkedFlag = (char *) mymalloc(NumPart * sizeof(char));
  ChangedFlag = (char *) mymalloc(NumPart * sizeof(char));
  MinIDOld = (int *) mymalloc(NumPart * sizeof(int));
  
  t0 = second();

  /* first, link only among local particles */
  for(i = 0, marked = 0, npart = 0; i < NumPart; i++)
    {
      if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
	{
	  fof_find_dmparticles_evaluate(i, -1, &dummy, &dummy);

	  npart++;

	  if(NonlocalFlag[i])
	    marked++;
	}
    }

  MPI_Allreduce(&marked, &totmarked, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&npart, &totnpart, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  t1 = second();


  if(ThisTask == 0)
    {
      printf("links on local processor done (took %g sec).\nMarked=%d out of the %d primaries which are linked\n", 
             timediff(t0, t1), totmarked,  totnpart);
          
      printf("\nlinking across processors (presently allocated=%g MB) \n",
	     AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  for(i = 0; i < NumPart; i++)
    {
      MinIDOld[i] = MinID[Head[i]];
      MarkedFlag[i] = 1;
    }

  do
    {
      t0 = second();

      for(i = 0; i < NumPart; i++)
        {
          ChangedFlag[i] = MarkedFlag[i];
          MarkedFlag[i] = 0;
        }

      i = 0;			/* begin with this index */
      link_across = 0;
      nprocessed = 0;      

      do
        {
          for(j = 0; j < NTask; j++)
            {
              Send_count[j] = 0;
              Exportflag[j] = -1;
            }
          
          /* do local particles and prepare export list */
          for(nexport = 0; i < NumPart; i++)
            {
              if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
                {
                  if(NonlocalFlag[i] && ChangedFlag[i])
		    {
		      if(fof_find_dmparticles_evaluate(i, 0, &nexport, Send_count) < 0)
			break;

		      nprocessed++;
		    }
                }
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


          FoFDataGet = (struct fofdata_in *) mymalloc(nimport * sizeof(struct fofdata_in));
          FoFDataIn = (struct fofdata_in *) mymalloc(nexport * sizeof(struct fofdata_in));
          

          /* prepare particle data for export */
          for(j = 0; j < nexport; j++)
            {
              place = DataIndexTable[j].Index;
              
              FoFDataIn[j].Pos[0] = P[place].Pos[0];
              FoFDataIn[j].Pos[1] = P[place].Pos[1];
              FoFDataIn[j].Pos[2] = P[place].Pos[2];
              FoFDataIn[j].MinID = MinID[Head[place]];
              FoFDataIn[j].MinIDTask = MinIDTask[Head[place]];
              
              memcpy(FoFDataIn[j].NodeList,
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
                      MPI_Sendrecv(&FoFDataIn[Send_offset[recvTask]],
                                   Send_count[recvTask] * sizeof(struct fofdata_in), MPI_BYTE,
                                   recvTask, TAG_DENS_A,
                                   &FoFDataGet[Recv_offset[recvTask]],
                                   Recv_count[recvTask] * sizeof(struct fofdata_in), MPI_BYTE,
                                   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }

          myfree(FoFDataIn);
          FoFDataResult = (char *) mymalloc(nimport * sizeof(char));
	  FoFDataOut = (char *) mymalloc(nexport * sizeof(char));
              
          /* now do the particles that were sent to us */
          
          for(j = 0; j < nimport; j++)
            {
              link_count = fof_find_dmparticles_evaluate(j, 1, &dummy, &dummy);
              link_across += link_count;
              if(link_count)
                FoFDataResult[j] = 1;
              else
                FoFDataResult[j] = 0;
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
                      /* get the particles */
		      MPI_Sendrecv(&FoFDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(char),
				   MPI_BYTE, recvTask, TAG_DENS_B,
				   &FoFDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(char),
				   MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }

	  /* need to mark the particle if it induced a link */
	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;
              if(FoFDataOut[j])
                MarkedFlag[place] = 1;
            }          
          
          myfree(FoFDataOut);
          myfree(FoFDataResult);
          myfree(FoFDataGet);
          
          if(i >= NumPart)
            ndone_flag = 1;
          else
            ndone_flag = 0;
          
          MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        }
      while(ndone < NTask);
      
      MPI_Allreduce(&link_across, &link_across_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&nprocessed, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      
      t1 = second();

      if(ThisTask == 0)
        {
          printf("have done %d cross links (processed %d, took %g sec)\n", link_across_tot, ntot, timediff(t0, t1));
          fflush(stdout);
        }
      

      /* let's check out which particles have changed their MinID */
      for(i=0; i<NumPart; i++)
        if(NonlocalFlag[i])
          {
            if(MinID[Head[i]] != MinIDOld[i])
              MarkedFlag[i] = 1;

            MinIDOld[i] = MinID[Head[i]];
          }
      
    }
  while(link_across_tot > 0);

  myfree(MinIDOld);
  myfree(ChangedFlag);
  myfree(MarkedFlag);
  myfree(NonlocalFlag);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

  if(ThisTask == 0)
    {
      printf("Local groups found.\n\n");
      fflush(stdout);
    }
}


int fof_find_dmparticles_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n, links, p, s, ss, listindex = 0;
  int startnode, numngb_inbox;
  MyDouble *pos;

  links = 0;

  if(mode == 0 || mode == -1)
    pos = P[target].Pos;
  else
    pos = FoFDataGet[target].Pos;

  if(mode == 0 || mode == -1)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = FoFDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  if(mode == -1)
	    *nexport = 0;

	  numngb_inbox = ngb_treefind_fof_primary(pos, LinkL, target, &startnode, mode, nexport, nsend_local);

	  if(numngb_inbox < 0)
	    return -1;

	  if(mode == -1)
	    {
	      if(*nexport == 0)
		NonlocalFlag[target] = 0;
	      else
		NonlocalFlag[target] = 1;
	    }

	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = Ngblist[n];

	      if(mode == 0 || mode == -1)
		{
		  if(Head[target] != Head[j])	/* only if not yet linked */
		    {
		      if(Len[Head[target]] > Len[Head[j]])	/* p group is longer */
			{
			  p = target;
			  s = j;
			}
		      else
			{
			  p = j;
			  s = target;
			}
		      Next[Tail[Head[p]]] = Head[s];

		      Tail[Head[p]] = Tail[Head[s]];

		      Len[Head[p]] += Len[Head[s]];

		      ss = Head[s];
		      do
			Head[ss] = Head[p];
		      while((ss = Next[ss]) >= 0);

		      if(MinID[Head[s]] < MinID[Head[p]])
			{
			  MinID[Head[p]] = MinID[Head[s]];
			  MinIDTask[Head[p]] = MinIDTask[Head[s]];
			}
		    }
		}
	      else		/* mode is 1 */
		{
		  if(MinID[Head[j]] > FoFDataGet[target].MinID)
		    {
		      MinID[Head[j]] = FoFDataGet[target].MinID;
		      MinIDTask[Head[j]] = FoFDataGet[target].MinIDTask;
		      links++;
		    }
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = FoFDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  return links;
}


void fof_compile_catalogue(void)
{
  int i, j, start, nimport, ngrp, sendTask, recvTask;
  struct fof_group_list *get_FOF_GList;

  /* sort according to MinID */
  qsort(FOF_PList, NumPart, sizeof(struct fof_particle_list), fof_compare_FOF_PList_MinID);

  for(i = 0; i < NumPart; i++)
    {
      FOF_GList[i].MinID = FOF_PList[i].MinID;
      FOF_GList[i].MinIDTask = FOF_PList[i].MinIDTask;
      if(FOF_GList[i].MinIDTask == ThisTask)
	{
	  FOF_GList[i].LocCount = 1;
	  FOF_GList[i].ExtCount = 0;
	}
      else
	{
	  FOF_GList[i].LocCount = 0;
	  FOF_GList[i].ExtCount = 1;
	}
    }

  /* eliminate duplicates in FOF_GList with respect to MinID */

  if(NumPart)
    NgroupsExt = 1;
  else
    NgroupsExt = 0;

  for(i = 1, start = 0; i < NumPart; i++)
    {
      if(FOF_GList[i].MinID == FOF_GList[start].MinID)
	{
	  FOF_GList[start].LocCount += FOF_GList[i].LocCount;
	  FOF_GList[start].ExtCount += FOF_GList[i].ExtCount;
	}
      else
	{
	  start = NgroupsExt;
	  FOF_GList[start] = FOF_GList[i];
	  NgroupsExt++;
	}
    }


  /* sort the remaining ones according to task */
  qsort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinIDTask);

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(i = 0; i < NgroupsExt; i++)
    Send_count[FOF_GList[i].MinIDTask]++;

  MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
      if(j == ThisTask)		/* we will not exchange the ones that are local */
	Recv_count[j] = 0;
      nimport += Recv_count[j];

      if(j > 0)
	{
	  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	}
    }

  get_FOF_GList = (struct fof_group_list *) mymalloc(nimport * sizeof(struct fof_group_list));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the group info */
	      MPI_Sendrecv(&FOF_GList[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct fof_group_list), MPI_BYTE,
			   recvTask, TAG_DENS_A,
			   &get_FOF_GList[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct fof_group_list), MPI_BYTE,
			   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  for(i = 0; i < nimport; i++)
    get_FOF_GList[i].MinIDTask = i;


  /* sort the groups according to MinID */
  qsort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinID);
  qsort(get_FOF_GList, nimport, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinID);

  /* merge the imported ones with the local ones */
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(FOF_GList[start].MinID < get_FOF_GList[i].MinID)
	{
	  start++;
	  if(start >= NgroupsExt)
	    endrun(7973);
	}

      if(get_FOF_GList[i].LocCount != 0)
	endrun(123);

      if(FOF_GList[start].MinIDTask != ThisTask)
	endrun(124);

      FOF_GList[start].ExtCount += get_FOF_GList[i].ExtCount;
    }

  /* copy the size information back into the list, to inform the others */
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(FOF_GList[start].MinID < get_FOF_GList[i].MinID)
	{
	  start++;
	  if(start >= NgroupsExt)
	    endrun(797831);
	}

      get_FOF_GList[i].ExtCount = FOF_GList[start].ExtCount;
      get_FOF_GList[i].LocCount = FOF_GList[start].LocCount;
    }

  /* sort the imported/exported list according to MinIDTask */
  qsort(get_FOF_GList, nimport, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinIDTask);
  qsort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinIDTask);


  for(i = 0; i < nimport; i++)
    get_FOF_GList[i].MinIDTask = ThisTask;

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the group info */
	      MPI_Sendrecv(&get_FOF_GList[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct fof_group_list), MPI_BYTE,
			   recvTask, TAG_DENS_A,
			   &FOF_GList[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct fof_group_list), MPI_BYTE,
			   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  myfree(get_FOF_GList);

  /* eliminate all groups that are too small, and count local groups */
  for(i = 0, Ngroups = 0, Nids = 0; i < NgroupsExt; i++)
    {
      if(FOF_GList[i].LocCount + FOF_GList[i].ExtCount < GROUP_MIN_LEN)
	{
	  FOF_GList[i] = FOF_GList[NgroupsExt - 1];
	  NgroupsExt--;
	  i--;
	}
      else
	{
	  if(FOF_GList[i].MinIDTask == ThisTask)
	    {
	      Ngroups++;
	      Nids += FOF_GList[i].LocCount + FOF_GList[i].ExtCount;
	    }
	}
    }

  /* sort the group list according to MinID */
  qsort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinID);
}


void fof_compute_group_properties(int gr, int start, int len)
{
  int j, k, index;
  double xyz[3];
#ifdef BG_SFR
  int i;
  double initial_mass, birth_aexp;
#endif

  Group[gr].Len = 0;
  Group[gr].Mass = 0;
#ifdef BG_SFR
  Group[gr].Sfr = 0;
#endif

#ifdef BLACK_HOLES
  Group[gr].BH_Mass = 0;
  Group[gr].BH_Mdot = 0;
  Group[gr].index_maxdens = Group[gr].task_maxdens = -1;
  Group[gr].MaxDens = 0;
#endif

  for(k = 0; k < 3; k++)
    {
      Group[gr].CM[k] = 0;
      Group[gr].Vel[k] = 0;
      Group[gr].FirstPos[k] = P[FOF_PList[start].Pindex].Pos[k];
      /* Group[gr].FirstPos[k] = P[start].Pos[k]; */
    }

  for(k = 0; k < 6; k++)
    {
      Group[gr].LenType[k] = 0;
      Group[gr].MassType[k] = 0;
    }

#ifdef BG_SFR /* BG_SFR */
  Group[gr].InitialMassWeightedStellarAge = 0;
  Group[gr].InitialMassWeightedStellarBirthZ = 0;
  Group[gr].StarInitialMass = 0;

  Group[gr].MassWeightedTemperature_sf = 0;
  Group[gr].MassWeightedTemperature_nsf = 0;
  Group[gr].MassWeightedEntropy_sf = 0;
  Group[gr].MassWeightedEntropy_nsf = 0;

  Group[gr].SFMass = 0;
  Group[gr].NSFMass = 0;
  Group[gr].StarMass = 0;

#ifdef BG_STELLAR_EVOLUTION /* BG_STELLAR_EVOLUTION */
  Group[gr].SFMetalMass = 0;
  Group[gr].NSFMetalMass = 0;
  Group[gr].StarMetalMass = 0;

  for(k = 0; k < BG_NELEMENTS; k++)
    {
      Group[gr].SFMetals[k] = 0;
      Group[gr].NSFMetals[k] = 0;
      Group[gr].StarMetals[k] = 0;
    }

#ifdef BG_METALSMOOTHING
  Group[gr].SFMetalMassSmoothed = 0;
  Group[gr].NSFMetalMassSmoothed = 0;
  Group[gr].StarMetalMassSmoothed = 0;

  for(k = 0; k < BG_NELEMENTS; k++)
    {
      Group[gr].SFMetalsSmoothed[k] = 0;
      Group[gr].NSFMetalsSmoothed[k] = 0;
      Group[gr].StarMetalsSmoothed[k] = 0;
    }
#endif

#ifdef BG_SNIA_IRON
  Group[gr].SFIronFromSNIa = 0;
  Group[gr].NSFIronFromSNIa = 0;
  Group[gr].StarIronFromSNIa = 0;

#ifdef BG_METALSMOOTHING
  Group[gr].SFIronFromSNIaSmoothed = 0;
  Group[gr].NSFIronFromSNIaSmoothed = 0;
  Group[gr].StarIronFromSNIaSmoothed = 0;
#endif
#endif

#ifdef BG_Z_WEIGHTED_POTENTIAL
  Group[gr].MetallicityWeightedPotential_sf = 0;
  Group[gr].MetallicityWeightedPotential_nsf = 0;
  Group[gr].MetallicityWeightedPotential_stars = 0;
#endif

#ifdef BG_Z_WEIGHTED_REDSHIFT
  Group[gr].MetallicityWeightedRedshift_sf = 0;
  Group[gr].MetallicityWeightedRedshift_nsf = 0;
  Group[gr].MetallicityWeightedRedshift_stars = 0;
#endif
#endif /* BG_STELLAR_EVOLUTION */

#ifdef BG_EXTRA_ARRAYS
  Group[gr].MassWeightedMaxTemp_sf = 0;
  Group[gr].MassWeightedMaxEntr_sf = 0;
  Group[gr].MassWeightedMaxTempAExp_sf = 0;
  Group[gr].MassWeightedMaxEntrAExp_sf = 0;
  Group[gr].MassWeightedMaxTemp_nsf = 0;
  Group[gr].MassWeightedMaxEntr_nsf = 0;
  Group[gr].MassWeightedMaxTempAExp_nsf = 0;
  Group[gr].MassWeightedMaxEntrAExp_nsf = 0;
  Group[gr].MassWeightedMaxTemp_stars = 0;
  Group[gr].MassWeightedMaxEntr_stars = 0;
  Group[gr].MassWeightedMaxTempAExp_stars = 0;
  Group[gr].MassWeightedMaxEntrAExp_stars = 0;
#endif

#ifdef EVALPOTENTIAL
  Group[gr].MassWeightedPotential_sf = 0;
  Group[gr].MassWeightedPotential_nsf = 0;
  Group[gr].MassWeightedPotential_stars = 0;
#endif
#endif /* BG_SFR */


  for(k = 0; k < len; k++)
    {
      index = FOF_PList[start + k].Pindex;

      Group[gr].Len++;
      Group[gr].Mass += P[index].Mass;
      Group[gr].LenType[P[index].Type]++;
      Group[gr].MassType[P[index].Type] += P[index].Mass;

#ifdef BG_SFR
      if(P[index].Type == 0)
	Group[gr].Sfr += SphP[index].Sfr;
      if(P[index].Type == 4)
	{
	  birth_aexp = StarP[P[index].StarID].StarBirthTime;
	  initial_mass = StarP[P[index].StarID].InitialMass;

	  Group[gr].InitialMassWeightedStellarAge += initial_mass *  bg_get_elapsed_time(birth_aexp, All.Time, 1);
	  Group[gr].InitialMassWeightedStellarBirthZ += initial_mass * (1 / birth_aexp - 1);

	  Group[gr].StarInitialMass += StarP[P[index].StarID].InitialMass;
	}
#endif

#ifdef BLACK_HOLES
      if(P[index].Type == 5)
	{
	  Group[gr].BH_Mdot += P[index].BH_Mdot;
	  Group[gr].BH_Mass += P[index].BH_Mass;
	}
      if(P[index].Type == 0)
	{
	  if(SphP[index].d.Density > Group[gr].MaxDens)
	    {
	      Group[gr].MaxDens = SphP[index].d.Density;
	      Group[gr].index_maxdens = index;
	      Group[gr].task_maxdens = ThisTask;
	    }
	}
#endif

      for(j = 0; j < 3; j++)
	{
	  xyz[j] = P[index].Pos[j];
#ifdef PERIODIC
	  xyz[j] = fof_periodic(xyz[j] - P[FOF_PList[start].Pindex].Pos[j]);
#endif
	  Group[gr].CM[j] += P[index].Mass * xyz[j];
	  Group[gr].Vel[j] += P[index].Mass * P[index].Vel[j];
	}


#ifdef BG_SFR /* BG_SFR */
      if(P[index].Type == 0)
	{
	  if(SphP[index].Sfr > 0)
	    {
	      Group[gr].SFMass += P[index].Mass;

	      Group[gr].MassWeightedTemperature_sf += bg_get_temperature(index) * P[index].Mass;
	      Group[gr].MassWeightedEntropy_sf += SphP[index].Entropy * P[index].Mass;

#ifdef BG_STELLAR_EVOLUTION /* BG_STELLAR_EVOLUTION */
	      Group[gr].SFMetalMass += SphP[index].Metallicity * P[index].Mass;

	      for(i = 0; i < BG_NELEMENTS; i++)
		Group[gr].SFMetals[i] += SphP[index].Metals[i];

#ifdef BG_METALSMOOTHING
	      Group[gr].SFMetalMassSmoothed += SphP[index].MetallicitySmoothed * P[index].Mass;

	      for(i = 0; i < BG_NELEMENTS; i++)
		Group[gr].SFMetalsSmoothed[i] += SphP[index].MetalsSmoothed[i];
#endif

#ifdef BG_SNIA_IRON
	      Group[gr].SFIronFromSNIa += SphP[index].IronFromSNIa;
#ifdef BG_METALSMOOTHING
	      Group[gr].SFIronFromSNIaSmoothed += SphP[index].IronFromSNIaSmoothed;
#endif
#endif

#ifdef BG_Z_WEIGHTED_POTENTIAL
	      Group[gr].MetallicityWeightedPotential_sf +=
		SphP[index].MetallicityWeightedPotential * SphP[index].Metallicity * P[index].Mass;
#endif

#ifdef BG_Z_WEIGHTED_REDSHIFT
	      Group[gr].MetallicityWeightedRedshift_sf +=
		SphP[index].MetallicityWeightedRedshift * SphP[index].Metallicity * P[index].Mass;
#endif
#endif /* BG_STELLAR_EVOLUTION */

#ifdef BG_EXTRA_ARRAYS
	      Group[gr].MassWeightedMaxTemp_sf += SphP[index].MaximumTemperature * P[index].Mass;
	      Group[gr].MassWeightedMaxEntr_sf += SphP[index].MaximumEntropy * P[index].Mass;
	      Group[gr].MassWeightedMaxTempAExp_sf += SphP[index].TimeMaximumTemperature * P[index].Mass;
	      Group[gr].MassWeightedMaxEntrAExp_sf += SphP[index].TimeMaximumEntropy * P[index].Mass;
#endif

#ifdef EVALPOTENTIAL
	      Group[gr].MassWeightedPotential_sf += P[index].p.Potential * P[index].Mass;
#endif
	    }
	  else
	    {
	      Group[gr].NSFMass += P[index].Mass;

	      Group[gr].MassWeightedTemperature_nsf += bg_get_temperature(index) * P[index].Mass;
	      Group[gr].MassWeightedEntropy_nsf += SphP[index].Entropy * P[index].Mass;

#ifdef BG_STELLAR_EVOLUTION /* BG_STELLAR_EVOLUTION */
	      Group[gr].NSFMetalMass += SphP[index].Metallicity * P[index].Mass;

	      for(i = 0; i < BG_NELEMENTS; i++)
		Group[gr].NSFMetals[i] += SphP[index].Metals[i];

#ifdef BG_METALSMOOTHING
	      Group[gr].NSFMetalMassSmoothed +=  SphP[index].MetallicitySmoothed * P[index].Mass;

	      for(i = 0; i < BG_NELEMENTS; i++)
		Group[gr].NSFMetalsSmoothed[i] += SphP[index].MetalsSmoothed[i];
#endif

#ifdef BG_SNIA_IRON
	      Group[gr].NSFIronFromSNIa += SphP[index].IronFromSNIa;
#ifdef BG_METALSMOOTHING
	      Group[gr].NSFIronFromSNIaSmoothed += SphP[index].IronFromSNIaSmoothed;
#endif
#endif

#ifdef BG_Z_WEIGHTED_POTENTIAL
	      Group[gr].MetallicityWeightedPotential_nsf +=
		SphP[index].MetallicityWeightedPotential * SphP[index].Metallicity * P[index].Mass;
#endif

#ifdef BG_Z_WEIGHTED_REDSHIFT
	      Group[gr].MetallicityWeightedRedshift_nsf +=
		SphP[index].MetallicityWeightedRedshift * SphP[index].Metallicity * P[index].Mass;
#endif
#endif /* BG_STELLAR_EVOLUTION */

#ifdef BG_EXTRA_ARRAYS
	      Group[gr].MassWeightedMaxTemp_nsf += SphP[index].MaximumTemperature * P[index].Mass;
	      Group[gr].MassWeightedMaxEntr_nsf += SphP[index].MaximumEntropy * P[index].Mass;
	      Group[gr].MassWeightedMaxTempAExp_nsf += SphP[index].TimeMaximumTemperature * P[index].Mass;
	      Group[gr].MassWeightedMaxEntrAExp_nsf += SphP[index].TimeMaximumEntropy * P[index].Mass;
#endif

#ifdef EVALPOTENTIAL
	      Group[gr].MassWeightedPotential_nsf += P[index].p.Potential * P[index].Mass;
#endif
	    }
	}

      if(P[index].Type == 4)
	{
	  Group[gr].StarMass += P[index].Mass;

#ifdef BG_STELLAR_EVOLUTION /* BG_STELLAR_EVOLUTION */
	  Group[gr].StarMetalMass += StarP[P[index].StarID].Metallicity * P[index].Mass;

	  for(i = 0; i < BG_NELEMENTS; i++)
	    Group[gr].StarMetals[i] += StarP[P[index].StarID].Metals[i];

#ifdef BG_METALSMOOTHING
	  Group[gr].StarMetalMassSmoothed +=  StarP[P[index].StarID].MetallicitySmoothed * P[index].Mass;

	  for(i = 0; i < BG_NELEMENTS; i++)
	    Group[gr].StarMetalsSmoothed[i] += StarP[P[index].StarID].MetalsSmoothed[i];
#endif

#ifdef BG_SNIA_IRON
	  Group[gr].StarIronFromSNIa += StarP[P[index].StarID].IronFromSNIa;
#ifdef BG_METALSMOOTHING
	  Group[gr].StarIronFromSNIaSmoothed += StarP[P[index].StarID].IronFromSNIaSmoothed;
#endif
#endif

#ifdef BG_Z_WEIGHTED_POTENTIAL
	  Group[gr].MetallicityWeightedPotential_stars +=
	    StarP[P[index].StarID].MetallicityWeightedPotential * StarP[P[index].StarID].Metallicity * P[index].Mass;
#endif

#ifdef BG_Z_WEIGHTED_REDSHIFT
	  Group[gr].MetallicityWeightedRedshift_stars +=
	    StarP[P[index].StarID].MetallicityWeightedRedshift * StarP[P[index].StarID].Metallicity * P[index].Mass;
#endif
#endif /* BG_STELLAR_EVOLUTION */

#ifdef BG_EXTRA_ARRAYS
	  Group[gr].MassWeightedMaxTemp_stars += StarP[P[index].StarID].MaximumTemperature * P[index].Mass;
	  Group[gr].MassWeightedMaxEntr_stars += StarP[P[index].StarID].MaximumEntropy * P[index].Mass;
	  Group[gr].MassWeightedMaxTempAExp_stars += StarP[P[index].StarID].TimeMaximumTemperature * P[index].Mass;
	  Group[gr].MassWeightedMaxEntrAExp_stars += StarP[P[index].StarID].TimeMaximumEntropy * P[index].Mass;
#endif

#ifdef EVALPOTENTIAL
	  Group[gr].MassWeightedPotential_stars += P[index].p.Potential * P[index].Mass;
#endif
    	}
#endif /* BG_SFR */
    }
}


void fof_exchange_group_data(void)
{
  struct group_properties *get_Group;
  int i, j, ngrp, sendTask, recvTask, nimport, start;
  double xyz[3];
#ifdef BG_SFR
  int k;
#endif

  /* sort the groups according to task */
  qsort(Group, NgroupsExt, sizeof(struct group_properties), fof_compare_Group_MinIDTask);

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(i = 0; i < NgroupsExt; i++)
    Send_count[FOF_GList[i].MinIDTask]++;

  MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
      if(j == ThisTask)		/* we will not exchange the ones that are local */
	Recv_count[j] = 0;
      nimport += Recv_count[j];

      if(j > 0)
	{
	  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	}
    }

  get_Group = (struct group_properties *) mymalloc(sizeof(struct group_properties) * nimport);

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the group data */
	      MPI_Sendrecv(&Group[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct group_properties), MPI_BYTE,
			   recvTask, TAG_DENS_A,
			   &get_Group[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct group_properties), MPI_BYTE,
			   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  /* sort the groups again according to MinID */
  qsort(Group, NgroupsExt, sizeof(struct group_properties), fof_compare_Group_MinID);
  qsort(get_Group, nimport, sizeof(struct group_properties), fof_compare_Group_MinID);

  /* now add in the partial imported group data to the main ones */
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(Group[start].MinID < get_Group[i].MinID)
	{
	  start++;
	  if(start >= NgroupsExt)
	    endrun(797890);
	}

      Group[start].Len += get_Group[i].Len;
      Group[start].Mass += get_Group[i].Mass;

      for(j = 0; j < 6; j++)
	{
	  Group[start].LenType[j] += get_Group[i].LenType[j];
	  Group[start].MassType[j] += get_Group[i].MassType[j];
	}

#ifdef BG_SFR
      Group[start].Sfr += get_Group[i].Sfr;
#endif
#ifdef BLACK_HOLES
      Group[start].BH_Mdot += get_Group[i].BH_Mdot;
      Group[start].BH_Mass += get_Group[i].BH_Mass;
      if(get_Group[i].MaxDens > Group[start].MaxDens)
	{
	  Group[start].MaxDens = get_Group[i].MaxDens;
	  Group[start].index_maxdens = get_Group[i].index_maxdens;
	  Group[start].task_maxdens = get_Group[i].task_maxdens;
	}
#endif

      for(j = 0; j < 3; j++)
	{
	  xyz[j] = get_Group[i].CM[j] / get_Group[i].Mass + get_Group[i].FirstPos[j];
#ifdef PERIODIC
	  xyz[j] = fof_periodic(xyz[j] - Group[start].FirstPos[j]);
#endif
	  Group[start].CM[j] += get_Group[i].Mass * xyz[j];
	  Group[start].Vel[j] += get_Group[i].Vel[j];
	}

#ifdef BG_SFR
      Group[start].InitialMassWeightedStellarAge += get_Group[i].InitialMassWeightedStellarAge;
      Group[start].InitialMassWeightedStellarBirthZ += get_Group[i].InitialMassWeightedStellarBirthZ;
      Group[start].StarInitialMass += get_Group[i].StarInitialMass;

#ifdef EVALPOTENTIAL
      Group[start].MassWeightedPotential_sf += get_Group[i].MassWeightedPotential_sf;
      Group[start].MassWeightedPotential_nsf += get_Group[i].MassWeightedPotential_nsf;
      Group[start].MassWeightedPotential_stars += get_Group[i].MassWeightedPotential_stars;
#endif

#ifdef BG_Z_WEIGHTED_POTENTIAL
      Group[start].MetallicityWeightedPotential_sf += get_Group[i].MetallicityWeightedPotential_sf;
      Group[start].MetallicityWeightedPotential_nsf += get_Group[i].MetallicityWeightedPotential_nsf;
      Group[start].MetallicityWeightedPotential_stars += get_Group[i].MetallicityWeightedPotential_stars;
#endif

#ifdef BG_Z_WEIGHTED_REDSHIFT
      Group[start].MetallicityWeightedRedshift_sf += get_Group[i].MetallicityWeightedRedshift_sf;
      Group[start].MetallicityWeightedRedshift_nsf += get_Group[i].MetallicityWeightedRedshift_nsf;
      Group[start].MetallicityWeightedRedshift_stars += get_Group[i].MetallicityWeightedRedshift_stars;
#endif
#endif /* BG_SFR */

#ifdef BG_STELLAR_EVOLUTION
      Group[start].SFMass += get_Group[i].SFMass;
      Group[start].NSFMass += get_Group[i].NSFMass;
      Group[start].StarMass += get_Group[i].StarMass;

      Group[start].SFMetalMass += get_Group[i].SFMetalMass;
      Group[start].NSFMetalMass += get_Group[i].NSFMetalMass;
      Group[start].StarMetalMass += get_Group[i].StarMetalMass;

      for(k = 0; k < BG_NELEMENTS; k++)
	{
	  Group[start].SFMetals[k] += get_Group[i].SFMetals[k];
	  Group[start].NSFMetals[k] += get_Group[i].NSFMetals[k];
	  Group[start].StarMetals[k] += get_Group[i].StarMetals[k];
	}

#ifdef BG_SNIA_IRON
      Group[start].SFIronFromSNIa += get_Group[i].SFIronFromSNIa;
      Group[start].NSFIronFromSNIa += get_Group[i].NSFIronFromSNIa;
      Group[start].StarIronFromSNIa += get_Group[i].StarIronFromSNIa;
#ifdef BG_METALSMOOTHING
      Group[start].SFIronFromSNIaSmoothed += get_Group[i].SFIronFromSNIaSmoothed;
      Group[start].NSFIronFromSNIaSmoothed += get_Group[i].NSFIronFromSNIaSmoothed;
      Group[start].StarIronFromSNIaSmoothed += get_Group[i].StarIronFromSNIaSmoothed;
#endif
#endif


#ifdef BG_METALSMOOTHING
      Group[start].SFMetalMassSmoothed += get_Group[i].SFMetalMassSmoothed;
      Group[start].NSFMetalMassSmoothed += get_Group[i].NSFMetalMassSmoothed;
      Group[start].StarMetalMassSmoothed += get_Group[i].StarMetalMassSmoothed;

      for(k = 0; k < BG_NELEMENTS; k++)
	{
	  Group[start].SFMetalsSmoothed[k] += get_Group[i].SFMetalsSmoothed[k];
	  Group[start].NSFMetalsSmoothed[k] += get_Group[i].NSFMetalsSmoothed[k];
	  Group[start].StarMetalsSmoothed[k] += get_Group[i].StarMetalsSmoothed[k];
	}
#endif

      Group[start].MassWeightedTemperature_sf += get_Group[i].MassWeightedTemperature_sf;
      Group[start].MassWeightedTemperature_nsf += get_Group[i].MassWeightedTemperature_nsf;

      Group[start].MassWeightedEntropy_sf += get_Group[i].MassWeightedEntropy_sf;
      Group[start].MassWeightedEntropy_nsf += get_Group[i].MassWeightedEntropy_nsf;

#ifdef BG_EXTRA_ARRAYS
      Group[start].MassWeightedMaxTemp_sf += get_Group[i].MassWeightedMaxTemp_sf;
      Group[start].MassWeightedMaxEntr_sf += get_Group[i].MassWeightedMaxEntr_sf;
      Group[start].MassWeightedMaxTempAExp_sf += get_Group[i].MassWeightedMaxTempAExp_sf;
      Group[start].MassWeightedMaxEntrAExp_sf += get_Group[i].MassWeightedMaxEntrAExp_sf;

      Group[start].MassWeightedMaxTemp_nsf += get_Group[i].MassWeightedMaxTemp_nsf;
      Group[start].MassWeightedMaxEntr_nsf += get_Group[i].MassWeightedMaxEntr_nsf;
      Group[start].MassWeightedMaxTempAExp_nsf += get_Group[i].MassWeightedMaxTempAExp_nsf;
      Group[start].MassWeightedMaxEntrAExp_nsf += get_Group[i].MassWeightedMaxEntrAExp_nsf;

      Group[start].MassWeightedMaxTemp_stars += get_Group[i].MassWeightedMaxTemp_stars;
      Group[start].MassWeightedMaxEntr_stars += get_Group[i].MassWeightedMaxEntr_stars;
      Group[start].MassWeightedMaxTempAExp_stars += get_Group[i].MassWeightedMaxTempAExp_stars;
      Group[start].MassWeightedMaxEntrAExp_stars += get_Group[i].MassWeightedMaxEntrAExp_stars ;
#endif

#endif /* BG_STELLAR_EVOLUTION */
    }


  for(i = 0; i < NgroupsExt; i++)
    {
      if(Group[i].MinIDTask == ThisTask)
	{
#ifdef EVALPOTENTIAL
	  if(Group[i].SFMass > 0)
	    Group[i].MassWeightedPotential_sf /= Group[i].SFMass;
	  if(Group[i].NSFMass > 0)
	    Group[i].MassWeightedPotential_nsf /= Group[i].NSFMass;
	  if(Group[i].StarMass > 0)
	    Group[i].MassWeightedPotential_stars /= Group[i].StarMass;
#endif

#ifdef BG_Z_WEIGHTED_POTENTIAL
	  if(Group[i].SFMetalMass > 0)
	    Group[i].MetallicityWeightedPotential_sf /= Group[i].SFMetalMass;
	  if(Group[i].NSFMetalMass > 0)
	    Group[i].MetallicityWeightedPotential_nsf /= Group[i].NSFMetalMass;
	  if(Group[i].StarMetalMass > 0)
	    Group[i].MetallicityWeightedPotential_stars /= Group[i].StarMetalMass;
#endif

#ifdef BG_Z_WEIGHTED_REDSHIFT
	  if(Group[i].SFMetalMass > 0)
	    Group[i].MetallicityWeightedRedshift_sf /= Group[i].SFMetalMass;
	  if(Group[i].NSFMetalMass > 0)
	    Group[i].MetallicityWeightedRedshift_nsf /= Group[i].NSFMetalMass;
	  if(Group[i].StarMetalMass > 0)
	    Group[i].MetallicityWeightedRedshift_stars /= Group[i].StarMetalMass;
#endif

#ifdef BG_EXTRA_ARRAYS
	  if(Group[i].SFMass > 0)
	    {
	      Group[i].MassWeightedMaxTemp_sf /= Group[i].SFMass;
	      Group[i].MassWeightedMaxEntr_sf /= Group[i].SFMass;
	      Group[i].MassWeightedMaxTempAExp_sf /= Group[i].SFMass;
	      Group[i].MassWeightedMaxEntrAExp_sf /= Group[i].SFMass;
	    }

	  if(Group[i].NSFMass > 0)
	    {
	      Group[i].MassWeightedMaxTemp_nsf /= Group[i].NSFMass;
	      Group[i].MassWeightedMaxEntr_nsf /= Group[i].NSFMass;
	      Group[i].MassWeightedMaxTempAExp_nsf /= Group[i].NSFMass;
	      Group[i].MassWeightedMaxEntrAExp_nsf /= Group[i].NSFMass;
	    }

	  if(Group[i].StarMass > 0)
	    {
	      Group[i].MassWeightedMaxTemp_stars /= Group[i].StarMass;
	      Group[i].MassWeightedMaxEntr_stars /= Group[i].StarMass;
	      Group[i].MassWeightedMaxTempAExp_stars /= Group[i].StarMass;
	      Group[i].MassWeightedMaxEntrAExp_stars /= Group[i].StarMass;
	    }
#endif

#ifdef BG_SFR
	  if(Group[i].StarInitialMass > 0)
	    {
	      Group[i].InitialMassWeightedStellarAge /= Group[i].StarInitialMass;
	      Group[i].InitialMassWeightedStellarBirthZ /= Group[i].StarInitialMass;
	    }
#endif

#ifdef BG_STELLAR_EVOLUTION
	  if(Group[i].SFMass > 0)
	    Group[i].SFMetalMass /= Group[i].SFMass;
	  if(Group[i].NSFMass > 0)
	    Group[i].NSFMetalMass /= Group[i].NSFMass;
	  if(Group[i].StarMass > 0)
	    Group[i].StarMetalMass /= Group[i].StarMass;

	  for(k = 0; k < BG_NELEMENTS; k++)
	    {
	      if(Group[i].SFMass > 0)
		Group[i].SFMetals[k] /= Group[i].SFMass;
	      if(Group[i].NSFMass > 0)
		Group[i].NSFMetals[k] /= Group[i].NSFMass;
	      if(Group[i].StarMass > 0)
		Group[i].StarMetals[k] /= Group[i].StarMass;
	    }

	  if(Group[i].SFMass > 0)
	    Group[i].MassWeightedTemperature_sf /= Group[i].SFMass;
	  if(Group[i].NSFMass > 0)
	    Group[i].MassWeightedTemperature_nsf /= Group[i].NSFMass;

	  if(Group[i].SFMass > 0)
	    Group[i].MassWeightedEntropy_sf /= Group[i].SFMass;
	  if(Group[i].NSFMass > 0)
	    Group[i].MassWeightedEntropy_nsf /= Group[i].NSFMass;
#endif
#ifdef BG_SNIA_IRON
	  if(Group[i].SFMass > 0)
	    Group[i].SFIronFromSNIa /= Group[i].SFMass;
	  if(Group[i].NSFMass > 0)
	    Group[i].NSFIronFromSNIa /= Group[i].NSFMass;
	  if(Group[i].StarMass > 0)
	    Group[i].StarIronFromSNIa /= Group[i].StarMass;

#ifdef BG_METALSMOOTHING
	  if(Group[i].SFMass > 0)
	    Group[i].SFIronFromSNIaSmoothed /= Group[i].SFMass;
	  if(Group[i].NSFMass > 0)
	    Group[i].NSFIronFromSNIaSmoothed /= Group[i].NSFMass;
	  if(Group[i].StarMass > 0)
	    Group[i].StarIronFromSNIaSmoothed /= Group[i].StarMass;
#endif
#endif

#ifdef BG_METALSMOOTHING
	  if(Group[i].SFMass > 0)
	    Group[i].SFMetalMassSmoothed /= Group[i].SFMass;
	  if(Group[i].NSFMass > 0)
	    Group[i].NSFMetalMassSmoothed /= Group[i].NSFMass;
	  if(Group[i].StarMass > 0)
	    Group[i].StarMetalMassSmoothed /= Group[i].StarMass;

	  for(k = 0; k < BG_NELEMENTS; k++)
	    {
	      if(Group[i].SFMass > 0)
		Group[i].SFMetalsSmoothed[k] /= Group[i].SFMass;
	      if(Group[i].NSFMass > 0)
		Group[i].NSFMetalsSmoothed[k] /= Group[i].NSFMass;
	      if(Group[i].StarMass > 0)
		Group[i].StarMetalsSmoothed[k] /= Group[i].StarMass;
	    }
#endif
	}
    }

  myfree(get_Group);
}


void fof_finish_group_properties(void)
{
  double cm[3];
  int i, j, ngr;

  for(i = 0; i < NgroupsExt; i++)
    {
      if(Group[i].MinIDTask == ThisTask)
	{
	  for(j = 0; j < 3; j++)
	    {
	      Group[i].Vel[j] /= Group[i].Mass;

	      cm[j] = Group[i].CM[j] / Group[i].Mass;
#ifdef PERIODIC
	      cm[j] = fof_periodic_wrap(cm[j] + Group[i].FirstPos[j]);
#endif
	      Group[i].CM[j] = cm[j];
	    }
	}
    }

  /* eliminate the non-local groups */
  for(i = 0, ngr = NgroupsExt; i < ngr; i++)
    {
      if(Group[i].MinIDTask != ThisTask)
	{
	  Group[i] = Group[ngr - 1];
	  i--;
	  ngr--;
	}
    }

  if(ngr != Ngroups)
    endrun(876889);

  qsort(Group, Ngroups, sizeof(struct group_properties), fof_compare_Group_MinID);
}


void fof_save_groups(int num)
{
  int i, j, k, start, lenloc, nprocgroup, masterTask, groupTask, ngr, totlen;
  long long totNids;
  char buf[500];
  double t0, t1;

  if(ThisTask == 0)
    {
      printf("start global sorting of group catalogues\n");
      fflush(stdout);
    }

  t0 = second();

  /* assign group numbers (at this point, both Group and FOF_GList are sorted by MinID) */
  for(i = 0; i < NgroupsExt; i++)
    {
      FOF_GList[i].LocCount += FOF_GList[i].ExtCount;	/* total length */
      FOF_GList[i].ExtCount = ThisTask;	/* original task */
    }

  parallel_sort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list),
		fof_compare_FOF_GList_LocCountTaskDiffMinID);

  for(i = 0, ngr = 0; i < NgroupsExt; i++)
    {
      if(FOF_GList[i].ExtCount == FOF_GList[i].MinIDTask)
	ngr++;

      FOF_GList[i].GrNr = ngr;
    }

  MPI_Allgather(&ngr, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  for(j = 1, Send_offset[0] = 0; j < NTask; j++)
    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];

  for(i = 0; i < NgroupsExt; i++)
    FOF_GList[i].GrNr += Send_offset[ThisTask];


  MPI_Allreduce(&ngr, &i, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(i != TotNgroups)
    {
      printf("i=%d\n", i);
      endrun(123123);
    }

  /* bring the group list back into the original order */
  parallel_sort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_ExtCountMinID);

  /* Assign the group numbers to the group properties array */
  for(i = 0, start = 0; i < Ngroups; i++)
    {
      while(FOF_GList[start].MinID < Group[i].MinID)
	{
	  start++;
	  if(start >= NgroupsExt)
	    endrun(7297890);
	}
      Group[i].GrNr = FOF_GList[start].GrNr;
    }

  /* sort the groups according to group-number */
  parallel_sort(Group, Ngroups, sizeof(struct group_properties), fof_compare_Group_GrNr);

  /* fill in the offset-values */
  for(i = 0, totlen = 0; i < Ngroups; i++)
    {
      if(i > 0)
	Group[i].Offset = Group[i - 1].Offset + Group[i - 1].Len;
      else
	Group[i].Offset = 0;
      totlen += Group[i].Len;
    }

  MPI_Allgather(&totlen, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  for(j = 1, Send_offset[0] = 0; j < NTask; j++)
    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];

  for(i = 0; i < Ngroups; i++)
    Group[i].Offset += Send_offset[ThisTask];


  /* fill in the offset-values per type */
  for(k = 0; k < 6; k++)
    {
      for(i = 0, totlen = 0; i < Ngroups; i++)
	{
	  if(i > 0)
	    Group[i].OffsetType[k] = Group[i - 1].OffsetType[k] + Group[i - 1].LenType[k];
	  else
	    Group[i].OffsetType[k] = 0;
	  totlen += Group[i].LenType[k];
	}

      MPI_Barrier(MPI_COMM_WORLD);

      MPI_Allgather(&totlen, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
      for(j = 1, Send_offset[0] = 0; j < NTask; j++)
	Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];

      for(i = 0; i < Ngroups; i++)
	Group[i].OffsetType[k] += Send_offset[ThisTask];
    }


  /* prepare list of ids with assigned group numbers */
  ID_list = mymalloc(sizeof(struct id_list) * NumPart);

  for(i = 0, start = 0, Nids = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
	{
	  start++;
	  if(start > NumPart)
	    endrun(78);
	}

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
	endrun(1313);

      for(lenloc = 0; start + lenloc < NumPart;)
	if(FOF_PList[start + lenloc].MinID == FOF_GList[i].MinID)
	  {
	    ID_list[Nids].GrNr = FOF_GList[i].GrNr;
	    ID_list[Nids].Type = P[FOF_PList[start + lenloc].Pindex].Type;
	    ID_list[Nids].ID = P[FOF_PList[start + lenloc].Pindex].ID;
	    Nids++;
	    lenloc++;
	  }
	else
	  break;

      start += lenloc;
    }

  sumup_large_ints(1, &Nids, &totNids);

  MPI_Allgather(&Nids, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  for(j = 1, Send_offset[0] = 0; j < NTask; j++)
    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];


  if(totNids != TotNids)
    {
      printf("Task=%d Nids=%d totNids=%d TotNids=%d\n", ThisTask, Nids, (int) totNids, (int) TotNids);
      endrun(12);
    }


  if(ThisTask == 0)
    {
	printf("before ID sort  FreeBytes=%g MB  AllocatedBytes=%g MB.\n",
	       FreeBytes / (1024.0 * 1024.0),  AllocatedBytes / (1024.0 * 1024.0));
	fflush(stdout);
    }

  /* parallel sort IDs */
  parallel_sort(ID_list, Nids, sizeof(struct id_list), fof_compare_ID_list_GrNrTypeID);


  t1 = second();
  if(ThisTask == 0)
    {
      printf("Group catalogues globally sorted. took = %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }


  if(ThisTask == 0)
    {
      printf("starting saving of group catalogue\n");
      fflush(stdout);
    }

  t0 = second();

  if(ThisTask == 0)
    {
       sprintf(buf, "%s/groups_%03d", All.OutputDir, num);
       mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);


  if(NTask < All.NumFilesWrittenInParallel)
    {
      printf
	("Fatal error.\nNumber of processors must be smaller or equal than `NumFilesWrittenInParallel'.\n");
      endrun(241931);
    }


  /* Fill global header */

  int ntot_type[6];
  long long ntot_type_all[6];

  /* determine global and local particle numbers */
  for(i = 0; i < 6; i++)
    ntot_type[i] = 0;

  for(i = 0; i < NumPart; i++)
    ntot_type[P[i].Type]++;

  /* fill file header */
  for(i = 0; i < 6; i++)
    {
      header.npart[i] = ntot_type[i];
      header.npartTotal[i] = (unsigned int) ntot_type_all[i];
      header.npartTotalHighWord[i] = (unsigned int) (ntot_type_all[i] >> 32);
    }

  for(i = 0; i < 6; i++)
    header.mass[i] = All.MassTable[i];

  header.time = All.Time;

  if(All.ComovingIntegrationOn)
    header.redshift = 1.0 / All.Time - 1;
  else
    header.redshift = 0;

  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;
  header.flag_cooling = 1;
#ifdef BG_SFR
  header.flag_sfr = 1;
  header.flag_feedback = 1;
  header.flag_stellarage = 1;
  header.flag_metals = 1;
#endif
  header.num_files = All.NumFilesPerSnapshot;
  header.BoxSize = All.BoxSize;
  header.Omega0 = All.Omega0;
  header.OmegaLambda = All.OmegaLambda;
  header.HubbleParam = All.HubbleParam;

  nprocgroup = NTask / All.NumFilesWrittenInParallel;
  if((NTask % All.NumFilesWrittenInParallel))
    nprocgroup++;
  masterTask = (ThisTask / nprocgroup) * nprocgroup;
  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))	/* ok, it's this processor's turn */
	{
	  fof_save_local_catalogue(num);
	}

      MPI_Barrier(MPI_COMM_WORLD);	/* wait inside the group */
    }

  myfree(ID_list);

  t1 = second();

  if(ThisTask == 0)
    {
      printf("Group catalogues saved. took = %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }
}


void fof_save_local_catalogue(int num)
{
  FILE *fd;
  float *mass, *cm, *vel;
  char fname[500];
  int i, j, *len;
  MyIDType *ids;
#ifdef HAVE_HDF5
  char buf[500];

  hid_t file = 0;
  hid_t headergrp = 0, paramgrp = 0;
  hid_t constgrp = 0, unitsgrp = 0;
  hid_t groupgrp = 0;
  hid_t dataspace = 0, attribute = 0;
#ifdef BG_SFR
  hid_t elementgrp = 0, typegrp = 0;
#endif
#endif


#ifndef  HAVE_HDF5
  if(All.SnapFormat == 3)
    {
      if(ThisTask == 0)
	printf("Code wasn't compiled with HDF5 support enabled!\n");
      endrun(0);
    }
#endif


  if(All.SnapFormat != 3)
    {
#ifdef INCL_DMONLY
      sprintf(fname, "%s/groups_%03d_DMONLY/%s_%03d.%d", All.OutputDir, num, "group_tab", num, ThisTask);
#else
      sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_tab", num, ThisTask);
#endif
      if(!(fd = fopen(fname, "w")))
	{
	  printf("can't open file `%s`\n", fname);
	  endrun(1183);
	}

      my_fwrite(&Ngroups, sizeof(int), 1, fd);
      my_fwrite(&TotNgroups, sizeof(int), 1, fd);
      my_fwrite(&Nids, sizeof(int), 1, fd);
      my_fwrite(&TotNids, sizeof(long long), 1, fd);
      my_fwrite(&NTask, sizeof(int), 1, fd);

      /* group len */
      len = mymalloc(Ngroups * sizeof(int));
      for(i = 0; i < Ngroups; i++)
	len[i] = Group[i].Len;
      my_fwrite(len, Ngroups, sizeof(int), fd);
      myfree(len);

      /* offset into id-list */
      len = mymalloc(Ngroups * sizeof(int));
      for(i = 0; i < Ngroups; i++)
	len[i] = Group[i].Offset;
      my_fwrite(len, Ngroups, sizeof(int), fd);
      myfree(len);

      /* mass */
      mass = mymalloc(Ngroups * sizeof(float));
      for(i = 0; i < Ngroups; i++)
	mass[i] = Group[i].Mass;
      my_fwrite(len, Ngroups, sizeof(float), fd);
      myfree(mass);

      /* CM */
      cm = mymalloc(Ngroups * 3 * sizeof(float));
      for(i = 0; i < Ngroups; i++)
	for(j = 0; j < 3; j++)
	  cm[i * 3 + j] = Group[i].CM[j];
      my_fwrite(cm, Ngroups, 3 * sizeof(float), fd);
      myfree(cm);

      /* vel */
      vel = mymalloc(Ngroups * 3 * sizeof(float));
      for(i = 0; i < Ngroups; i++)
	for(j = 0; j < 3; j++)
	  vel[i * 3 + j] = Group[i].Vel[j];
      my_fwrite(vel, Ngroups, 3 * sizeof(float), fd);
      myfree(vel);

      /* group len for each type */
      len = mymalloc(Ngroups * 6 * sizeof(int));
      for(i = 0; i < Ngroups; i++)
	for(j = 0; j < 6; j++)
	  len[i * 6 + j] = Group[i].LenType[j];
      my_fwrite(len, Ngroups, 6 * sizeof(int), fd);
      myfree(len);

      /* group mass for each type */
      mass = mymalloc(Ngroups * 6 * sizeof(int));
      for(i = 0; i < Ngroups; i++)
	for(j = 0; j < 6; j++)
	  mass[i * 6 + j] = Group[i].MassType[j];
      my_fwrite(mass, Ngroups, 6 * sizeof(float), fd);
      myfree(mass);

#ifdef BLACK_HOLES
      /* BH_Mass */
      mass = mymalloc(Ngroups * sizeof(float));
      for(i = 0; i < Ngroups; i++)
	mass[i] = Group[i].BH_Mass;
      my_fwrite(mass, Ngroups, sizeof(float), fd);
      myfree(mass);

      /* BH_Mdot */
      mass = mymalloc(Ngroups * sizeof(float));
      for(i = 0; i < Ngroups; i++)
	mass[i] = Group[i].BH_Mdot;
      my_fwrite(mass, Ngroups, sizeof(float), fd);
      myfree(mass);
#endif

#ifdef BG_SFR
      /* sfr */
      mass = mymalloc(Ngroups * sizeof(float));
      for(i = 0; i < Ngroups; i++)
	mass[i] = Group[i].Sfr;
      my_fwrite(mass, Ngroups, sizeof(float), fd);
      myfree(mass);
#endif


      fclose(fd);


      ids = (MyIDType *) mymalloc(Nids * sizeof(MyIDType));
      for(i = 0; i < Nids; i++)
	ids[i] = ID_list[i].ID;

      sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_ids", num, ThisTask);

      if(!(fd = fopen(fname, "w")))
	{
	  printf("can't open file `%s`\n", fname);
	  endrun(1184);
	}

      my_fwrite(&Ngroups, sizeof(int), 1, fd);
      my_fwrite(&TotNgroups, sizeof(int), 1, fd);
      my_fwrite(&Nids, sizeof(int), 1, fd);
      my_fwrite(&TotNids, sizeof(long long), 1, fd);
      my_fwrite(&NTask, sizeof(int), 1, fd);
      my_fwrite(&Send_offset[ThisTask], sizeof(int), 1, fd);	/* this is the number of IDs in previous files */
      my_fwrite(ids, sizeof(int), Nids, fd);
      fclose(fd);
    }
  else
    {
      sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group", num, ThisTask);
      sprintf(buf, "%s.hdf5", fname);

      file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      /* write usual headers */
      headergrp = H5Gcreate(file, "/Header", 0);
      paramgrp = H5Gcreate(file, "/Parameters", 0);
      unitsgrp = H5Gcreate(file, "/Units", 0);
      constgrp = H5Gcreate(file, "/Constants", 0);

      if(ThisTask == 0)
	printf("Writing Header\n");

      write_header_attributes_in_hdf5(headergrp);

      if(ThisTask == 0)
	printf("Writing Parameters\n");

      write_parameters_attributes_in_hdf5(paramgrp);

      if(ThisTask == 0)
	printf("Writing Units\n");

      write_units_attributes_in_hdf5(unitsgrp);

      if(ThisTask == 0)
	printf("Writing Constants\n");

      write_constants_attributes_in_hdf5(constgrp);

      H5Gclose(headergrp);
      H5Gclose(paramgrp);
      H5Gclose(unitsgrp);
      H5Gclose(constgrp);

      if(ThisTask == 0)
	printf("Writing FOF catalogue\n");

      /* specific for group files */
      groupgrp = H5Gcreate(file, "/FOF", 0);
      dataspace = H5Screate(H5S_SCALAR);


#ifdef FOF_OUTPUT_PID_ARRAY
      if(ThisTask == 0)
	printf(" --> Number_of_IDs\n");

      /* Nids */
      attribute = H5Acreate(groupgrp, "Number_of_IDs", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
      H5Awrite(attribute, H5T_NATIVE_INT, &Nids);
      H5Aclose(attribute);
#endif

      if(ThisTask == 0)
	printf(" --> Number_of_groups\n");

      /* Ngroups */
      attribute = H5Acreate(groupgrp, "Number_of_groups", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
      H5Awrite(attribute, H5T_NATIVE_INT, &Ngroups);
      H5Aclose(attribute);

      if(ThisTask == 0)
	printf(" --> Total_Number_of_groups\n");

      /* TotNgroups */
      attribute = H5Acreate(groupgrp, "Total_Number_of_groups", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
      H5Awrite(attribute, H5T_NATIVE_INT, &TotNgroups);
      H5Aclose(attribute);

      if(ThisTask == 0)
	printf(" --> NTask\n");

      /* NTask */
      attribute = H5Acreate(groupgrp, "NTask", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
      H5Awrite(attribute, H5T_NATIVE_INT, &NTask);
      H5Aclose(attribute);

      H5Sclose(dataspace);

#ifdef FOF_OUTPUT_PID_ARRAY
      if(ThisTask == 0)
	printf("Writing ParticleIDs\n");

      if(Nids)
	{
	  ids = (MyIDType *) mymalloc(Nids * sizeof(MyIDType));
	  for(i = 0; i < Nids; i++)
	    ids[i] = ID_list[i].ID;

	  /* write in the group ID's here */
	  fof_write_hdf5_myIDType_field(&ids[0], "ParticleIDs", Nids, "Unique particle identifier", groupgrp);

	  myfree(ids);
	}
#endif

      /* only if actually any groups found */
      if(Ngroups > 0)
	{  
	  if(ThisTask == 0)
	    printf("Writing FOF header\n");
	  fof_write_hdf5_int_field(&Group[0].Len, "Length", 1, "Number of particles in each group", groupgrp);

	  fof_write_hdf5_int_field(&Group[0].Offset,
				   "Offset", 1, "Particle number offset of each group", groupgrp);

	  fof_write_hdf5_int_field(&Group[0].LenType[0],
				   "LengthType", 6, "Number of particles of different type in each group",
				   groupgrp);
	  fof_write_hdf5_double_field(&Group[0].MassType[0], "MassType", 6, "", IO_MASS, groupgrp);

	  fof_write_hdf5_double_field(&Group[0].Mass, "Mass", 1, "", IO_MASS, groupgrp);

	  fof_write_hdf5_double_field(&Group[0].CM[0], "CenterOfMass", 3, "", IO_POS, groupgrp);

	  fof_write_hdf5_double_field(&Group[0].Vel[0], "CenterOfMassVelocity", 3, "", IO_VEL, groupgrp);

#if defined (BG_SFR) && defined (FOF_OUTPUT_OWLS_ARRAYS)
	  fof_write_hdf5_double_field(&Group[0].Sfr, "StarFormationRate", 1, "", IO_SFR, groupgrp);

	  if(ThisTask == 0)
	    printf("Writing /FOF/SF\n");

	  typegrp = H5Gcreate(file, "/FOF/SF", 0);

	  fof_write_hdf5_double_field(&Group[0].MassWeightedTemperature_sf,
				     "Temperature", 1, "", IO_BG_TEMP, typegrp);

	  fof_write_hdf5_double_field(&Group[0].MassWeightedEntropy_sf,
				     "Entropy", 1, "", IO_BG_MAX_ENTROPY, typegrp);

	  /* GroupSFMass */
	  fof_write_hdf5_double_field(&Group[0].SFMass, "Mass", 1, "", IO_MASS, typegrp);

#ifdef BG_STELLAR_EVOLUTION /* BG_STELLAR_EVOLUTION */
	  if(ThisTask == 0)
	    printf("Writing ElementAbundance\n");

	  /* GroupSFMetals */
	  elementgrp = H5Gcreate(typegrp, "ElementAbundance", 0);

	  for(i = 0; i < BG_NELEMENTS; i++)
	    fof_write_hdf5_double_field(&Group[0].SFMetals[i],
					ElementNames[i], 1, "", IO_BG_METALS, elementgrp);

	  H5Gclose(elementgrp);

	  /* GroupSFMetallicity */
	  fof_write_hdf5_double_field(&Group[0].SFMetalMass, "Metallicity", 1, "", IO_BG_METALLICITY, typegrp);

#ifdef BG_METALSMOOTHING
	  if(ThisTask == 0)
	    printf("Writing SmoothedElementAbundance\n");

	  /* GroupSFMetalsSmoothed */
	  elementgrp = H5Gcreate(typegrp, "SmoothedElementAbundance", 0);
	  for(i = 0; i < BG_NELEMENTS; i++)
	    fof_write_hdf5_double_field(&Group[0].SFMetalsSmoothed[i],
				       ElementNames[i], 1, "", IO_BG_METALS, elementgrp);
	  H5Gclose(elementgrp);

	  /* GroupSFMetallicitySmoothed */
	  fof_write_hdf5_double_field(&Group[0].SFMetalMassSmoothed,
				     "SmoothedMetallicity", 1, "", IO_BG_METALLICITY_SMOOTHED, typegrp);
#endif

#ifdef BG_SNIA_IRON
#ifdef BG_METALSMOOTHING
	  /* GroupSFIronFromSNIaSmoothed */
	  fof_write_hdf5_double_field(&Group[0].SFIronFromSNIaSmoothed,
				     "SmoothedIronFromSNIa", 1, "", IO_BG_IRON_FROM_SNIA_SMOOTHED, typegrp);
#endif
	  /* GroupSFIronFromSNIa */
	  fof_write_hdf5_double_field(&Group[0].SFIronFromSNIa,
				     "IronFromSNIa", 1, "", IO_BG_IRON_FROM_SNIA, typegrp);
#endif

#ifdef BG_Z_WEIGHTED_REDSHIFT
	  fof_write_hdf5_double_field(&Group[0].MetallicityWeightedRedshift_sf,
				     "MetallicityWeightedRedshift", 1, "",
				     IO_BG_METALLICITY_WEIGHTED_REDSHIFT, typegrp);
#endif

#ifdef BG_Z_WEIGHTED_POTENTIAL
	  fof_write_hdf5_double_field(&Group[0].MetallicityWeightedPotential_sf,
				     "MetallicityWeightedPotential", 1, "",
				     IO_BG_METALLICITY_WEIGHTED_POTENTIAL, typegrp);
#endif
#endif /* BG_STELLAR_EVOLUTION */

#ifdef BG_EXTRA_ARRAYS
	  fof_write_hdf5_double_field(&Group[0].MassWeightedMaxTemp_sf,
				     "MaximumTemperature", 1, "", IO_BG_MAX_TEMPERATURE, typegrp);

	  fof_write_hdf5_double_field(&Group[0].MassWeightedMaxEntr_sf,
				     "MaximumEntropy", 1, "", IO_BG_MAX_ENTROPY, typegrp);

	  fof_write_hdf5_double_field(&Group[0].MassWeightedMaxTempAExp_sf,
				     "AExpMaximumTemperature", 1, "", IO_BG_TIME_MAX_TEMPERATURE, typegrp);

	  fof_write_hdf5_double_field(&Group[0].MassWeightedMaxEntrAExp_sf,
				     "AExpMaximumEntropy", 1, "", IO_BG_TIME_MAX_ENTROPY, typegrp);
#endif

#ifdef EVALPOTENTIAL
	  fof_write_hdf5_double_field(&Group[0].MassWeightedPotential_sf,
				      "Potential", 1, "", IO_POT, typegrp);
#endif
	  H5Gclose(typegrp);

	  if(ThisTask == 0)
	    printf("Writing /FOF/NSF\n");

	  typegrp = H5Gcreate(file, "/FOF/NSF", 0);

	  fof_write_hdf5_double_field(&Group[0].MassWeightedTemperature_nsf,
				     "Temperature", 1, "", IO_BG_TEMP, typegrp);

	  fof_write_hdf5_double_field(&Group[0].MassWeightedEntropy_nsf,
				     "Entropy", 1, "", IO_BG_MAX_ENTROPY, typegrp);

	  /* GroupNsfMass */
	  fof_write_hdf5_double_field(&Group[0].NSFMass, "Mass", 1, "", IO_MASS, typegrp);

#ifdef BG_STELLAR_EVOLUTION /* BG_STELLAR_EVOLUTION */
	  if(ThisTask == 0)
	    printf("Writing ElementAbundance\n");

	  /* GroupNsfMetals */
	  elementgrp = H5Gcreate(typegrp, "ElementAbundance", 0);
	  for(i = 0; i < BG_NELEMENTS; i++)
	    fof_write_hdf5_double_field(&Group[0].NSFMetals[i],
				       ElementNames[i], 1, "", IO_BG_METALS, elementgrp);
	  H5Gclose(elementgrp);

	  /* GroupNSFMetallicity */
	  fof_write_hdf5_double_field(&Group[0].NSFMetalMass,
				     "Metallicity", 1, "", IO_BG_METALLICITY, typegrp);

#ifdef BG_METALSMOOTHING
	  if(ThisTask == 0)
	    printf("Writing SmoothedElementAbundance\n");

	  /* GroupNSFMetalsSmoothed */
	  elementgrp = H5Gcreate(typegrp, "SmoothedElementAbundance", 0);
	  for(i = 0; i < BG_NELEMENTS; i++)
	    fof_write_hdf5_double_field(&Group[0].NSFMetalsSmoothed[i],
				       ElementNames[i], 1, "", IO_BG_METALS, elementgrp);
	  H5Gclose(elementgrp);


	  /* GroupNSFMetallicitySmoothed */
	  fof_write_hdf5_double_field(&Group[0].NSFMetalMassSmoothed,
				     "SmoothedMetallicity", 1, "", IO_BG_METALLICITY_SMOOTHED, typegrp);
#endif

#ifdef BG_SNIA_IRON
#ifdef BG_METALSMOOTHING
	  /* GroupNSFIronFromSNIaSmoothed */
	  fof_write_hdf5_double_field(&Group[0].SFIronFromSNIaSmoothed,
				     "SmoothedIronFromSNIa", 1, "", IO_BG_IRON_FROM_SNIA_SMOOTHED, typegrp);
#endif
	  /* GroupNSFIronFromSNIa */
	  fof_write_hdf5_double_field(&Group[0].NSFIronFromSNIa,
				     "IronFromSNIa", 1, "", IO_BG_IRON_FROM_SNIA, typegrp);
#endif

#ifdef BG_Z_WEIGHTED_REDSHIFT
	  fof_write_hdf5_double_field(&Group[0].MetallicityWeightedRedshift_nsf,
				     "MetallicityWeightedRedshift", 1, "",
				     IO_BG_METALLICITY_WEIGHTED_REDSHIFT, typegrp);
#endif

#ifdef BG_Z_WEIGHTED_POTENTIAL
	  fof_write_hdf5_double_field(&Group[0].MetallicityWeightedPotential_nsf,
				     "MetallicityWeightedPotential", 1, "",
				     IO_BG_METALLICITY_WEIGHTED_POTENTIAL, typegrp);
#endif
#endif /* BG_STELLAR_EVOLUTION */

#ifdef BG_EXTRA_ARRAYS
	  fof_write_hdf5_double_field(&Group[0].MassWeightedMaxTemp_nsf,
				     "MaximumTemperature", 1, "", IO_BG_MAX_TEMPERATURE, typegrp);

	  fof_write_hdf5_double_field(&Group[0].MassWeightedMaxEntr_nsf,
				     "MaximumEntropy", 1, "", IO_BG_MAX_ENTROPY, typegrp);

	  fof_write_hdf5_double_field(&Group[0].MassWeightedMaxTempAExp_nsf,
				     "AExpMaximumTemperature", 1, "", IO_BG_TIME_MAX_TEMPERATURE, typegrp);

	  fof_write_hdf5_double_field(&Group[0].MassWeightedMaxEntrAExp_nsf,
				     "AExpMaximumEntropy", 1, "", IO_BG_TIME_MAX_ENTROPY, typegrp);
#endif

#ifdef EVALPOTENTIAL
	  fof_write_hdf5_double_field(&Group[0].MassWeightedPotential_nsf,
				     "Potential", 1, "", IO_POT, typegrp);
#endif
	  H5Gclose(typegrp);

	  if(ThisTask == 0)
	    printf("Writing /FOF/Stars\n");

	  typegrp = H5Gcreate(file, "/FOF/Stars", 0);

	  /* InitialMassWeightedStellarAge */
	  fof_write_hdf5_double_field(&Group[0].InitialMassWeightedStellarAge,
				     "InitialMassWeightedStellarAge", 1,
				     "", IO_BG_STELLAR_AGE, typegrp);

	  /* InitialMassWeightedStellarBirthZ */
	  fof_write_hdf5_double_field(&Group[0].InitialMassWeightedStellarBirthZ,
				     "InitialMassWeightedStellarBirthZ", 1,
				     "Mass weighted formation redshift of stars in the group", -1, typegrp);

	  fof_write_hdf5_double_field(&Group[0].StarMass, "Mass", 1, "", IO_MASS, typegrp);

	  fof_write_hdf5_double_field(&Group[0].StarInitialMass, "InitialMass", 1, "", IO_MASS, typegrp);

#ifdef BG_STELLAR_EVOLUTION /* BG_STELLAR_EVOLUTION */
	  if(ThisTask == 0)
	    printf("Writing ElementAbundance\n");

	  /* GroupStarMetals */
	  elementgrp = H5Gcreate(typegrp, "ElementAbundance", 0);
	  for(i = 0; i < BG_NELEMENTS; i++)
	    fof_write_hdf5_double_field(&Group[0].StarMetals[i],
				       ElementNames[i], 1, "", IO_BG_METALS, elementgrp);
	  H5Gclose(elementgrp);


	  fof_write_hdf5_double_field(&Group[0].StarMetalMass,
				     "Metallicity", 1, "", IO_BG_METALLICITY, typegrp);

#ifdef BG_METALSMOOTHING
	  if(ThisTask == 0)
	    printf("Writing SmoothedElementAbundance\n");

	  /* GroupStarMetalsSmoothed */
	  elementgrp = H5Gcreate(typegrp, "SmoothedElementAbundance", 0);
	  for(i = 0; i < BG_NELEMENTS; i++)
	    fof_write_hdf5_double_field(&Group[0].StarMetalsSmoothed[i],
				       ElementNames[i], 1, "", IO_BG_METALS, elementgrp);
	  H5Gclose(elementgrp);

	  /* GroupStarMetallicitySmoothed */
	  fof_write_hdf5_double_field(&Group[0].StarMetalMassSmoothed,
				     "SmoothedMetallicity", 1, "", IO_BG_METALLICITY_SMOOTHED, typegrp);
#endif

#ifdef BG_SNIA_IRON
#ifdef BG_METALSMOOTHING
	  /* GroupStarIronFromSNIaSmoothed */
	  fof_write_hdf5_double_field(&Group[0].StarIronFromSNIaSmoothed,
				     "SmoothedIronFromSNIa", 1, "", IO_BG_IRON_FROM_SNIA_SMOOTHED, typegrp);
#endif
	  /* GroupStarIronFromSNIa */
	  fof_write_hdf5_double_field(&Group[0].StarIronFromSNIa,
				     "IronFromSNIa", 1, "", IO_BG_IRON_FROM_SNIA, typegrp);
#endif

#ifdef BG_Z_WEIGHTED_REDSHIFT
	  fof_write_hdf5_double_field(&Group[0].MetallicityWeightedRedshift_stars,
				     "MetallicityWeightedRedshift", 1, "",
				     IO_BG_METALLICITY_WEIGHTED_REDSHIFT, typegrp);
#endif

#ifdef BG_Z_WEIGHTED_POTENTIAL
	  fof_write_hdf5_double_field(&Group[0].MetallicityWeightedPotential_stars,
				     "MetallicityWeightedPotential", 1, "",
				     IO_BG_METALLICITY_WEIGHTED_POTENTIAL, typegrp);
#endif
#endif /* BG_STELLAR_EVOLUTION */

#ifdef BG_EXTRA_ARRAYS
	  fof_write_hdf5_double_field(&Group[0].MassWeightedMaxTemp_stars,
				     "MaximumTemperature", 1, "", IO_BG_MAX_TEMPERATURE, typegrp);

	  fof_write_hdf5_double_field(&Group[0].MassWeightedMaxEntr_stars,
				     "MaximumEntropy", 1, "", IO_BG_MAX_ENTROPY, typegrp);

	  fof_write_hdf5_double_field(&Group[0].MassWeightedMaxTempAExp_stars,
				     "AExpMaximumTemperature", 1, "", IO_BG_TIME_MAX_TEMPERATURE, typegrp);

	  fof_write_hdf5_double_field(&Group[0].MassWeightedMaxEntrAExp_stars,
				     "AExpMaximumEntropy", 1, "", IO_BG_TIME_MAX_ENTROPY, typegrp);
#endif

#ifdef EVALPOTENTIAL
	  fof_write_hdf5_double_field(&Group[0].MassWeightedPotential_stars,
				      "Potential", 1, "", IO_POT, typegrp);
#endif
	  H5Gclose(typegrp);
#endif /* BG_SFR */


#ifdef BLACK_HOLES
          fof_write_hdf5_double_field(&Group[0].BH_Mass, "BH_Mass", 1, 
"", -1, groupgrp);

          fof_write_hdf5_double_field(&Group[0].BH_Mdot, "BH_Mdot", 1, 
"", -1, groupgrp);
#endif
	}			/* end of if(Ngroups > 0){} */

      /* Close group grp */
      H5Gclose(groupgrp);
      
      /* Close file */
      H5Fclose(file);
    }
}


#ifdef FOF_SECONDARY_LINK_TYPES
void fof_find_nearest_dmparticle(void)
{
  int i, j, n, ntot, dummy;
  int ndone, ndone_flag, ngrp, sendTask, recvTask, place, nexport, nimport, npleft, iter;

  if(ThisTask == 0)
    {
      printf("Start finding nearest dm-particle (presently allocated=%g MB)\n",
	     AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  fof_nearest_distance = (float *) mymalloc(sizeof(float) * NumPart);
  fof_nearest_hsml = (float *) mymalloc(sizeof(float) * NumPart);

  for(n = 0; n < NumPart; n++)
    {
      if(((1 << P[n].Type) & (FOF_SECONDARY_LINK_TYPES)))
	{
	  fof_nearest_distance[n] = 1.0e30;
	  if(P[n].Type == 0)
	    fof_nearest_hsml[n] = PPP[n].Hsml;
	  else
	    fof_nearest_hsml[n] = 0.1 * LinkL;
	}
    }

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct fofdata_in) + sizeof(struct fofdata_out) +
					     sizemax(sizeof(struct fofdata_in), sizeof(struct fofdata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));


  iter = 0;
  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      i = 0;			/* beginn with this index */

      do
	{
	  for(j = 0; j < NTask; j++)
	    {
	      Send_count[j] = 0;
	      Exportflag[j] = -1;
	    }

	  /* do local particles and prepare export list */
	  for(nexport = 0; i < NumPart; i++)
	    if(((1 << P[i].Type) & (FOF_SECONDARY_LINK_TYPES)))
	      {
		if(fof_nearest_distance[i] > 1.0e29)
		  {
		    if(fof_find_nearest_dmparticle_evaluate(i, 0, &nexport, Send_count) < 0)
		      break;
		  }
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

	  FoFDataGet = (struct fofdata_in *) mymalloc(nimport * sizeof(struct fofdata_in));
	  FoFDataIn = (struct fofdata_in *) mymalloc(nexport * sizeof(struct fofdata_in));

	  if(ThisTask == 0)
	    {
	      printf("still finding nearest... (presently allocated=%g MB)\n",
		     AllocatedBytes / (1024.0 * 1024.0));
	      fflush(stdout);
	    }

	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      FoFDataIn[j].Pos[0] = P[place].Pos[0];
	      FoFDataIn[j].Pos[1] = P[place].Pos[1];
	      FoFDataIn[j].Pos[2] = P[place].Pos[2];
	      FoFDataIn[j].Hsml = fof_nearest_hsml[place];

	      memcpy(FoFDataIn[j].NodeList,
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
		      MPI_Sendrecv(&FoFDataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct fofdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A,
				   &FoFDataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct fofdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }

	  myfree(FoFDataIn);
	  FoFDataResult = (struct fofdata_out *) mymalloc(nimport * sizeof(struct fofdata_out));
	  FoFDataOut = (struct fofdata_out *) mymalloc(nexport * sizeof(struct fofdata_out));

	  for(j = 0; j < nimport; j++)
	    {
	      fof_find_nearest_dmparticle_evaluate(j, 1, &dummy, &dummy);
	    }

	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&FoFDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct fofdata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B,
				   &FoFDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct fofdata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}

	    }

	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      if(FoFDataOut[j].Distance < fof_nearest_distance[place])
		{
		  fof_nearest_distance[place] = FoFDataOut[j].Distance;
		  MinID[place] = FoFDataOut[j].MinID;
		  MinIDTask[place] = FoFDataOut[j].MinIDTask;
		}
	    }

	  if(i >= NumPart)
	    ndone_flag = 1;
	  else
	    ndone_flag = 0;

	  MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	  myfree(FoFDataOut);
	  myfree(FoFDataResult);
	  myfree(FoFDataGet);
	}
      while(ndone < NTask);

      /* do final operations on results */
      for(i = 0, npleft = 0; i < NumPart; i++)
	{
	  if(((1 << P[i].Type) & (FOF_SECONDARY_LINK_TYPES)))
	    {
	      if(fof_nearest_distance[i] > 1.0e29)
		{
		  /* need to redo this particle */
		  npleft++;
		  fof_nearest_hsml[i] *= 2.0;
		  if(iter >= MAXITER - 10)
		    {
		      printf("i=%d task=%d ID=%d Hsml=%g  pos=(%g|%g|%g)\n",
			     i, ThisTask, (int) P[i].ID, fof_nearest_hsml[i],
			     P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
		      fflush(stdout);
		    }
		}
	    }
	}

      MPI_Allreduce(&npleft, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if(ntot > 0)
	{
	  iter++;
	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("fof-nearest iteration %d: need to repeat for %d particles.\n", iter, ntot);
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in fof-nearest\n");
	      fflush(stdout);
	      endrun(1159);
	    }
	}
    }
  while(ntot > 0);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

  myfree(fof_nearest_hsml);
  myfree(fof_nearest_distance);

  if(ThisTask == 0)
    {
      printf("done finding nearest dm-particle\n");
      fflush(stdout);
    }
}
#endif


#ifdef FOF_SECONDARY_LINK_TYPES
int fof_find_nearest_dmparticle_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n, index, listindex = 0;
  int startnode, numngb_inbox;
  double h, r2max;
  double dx, dy, dz, r2;
  MyDouble *pos;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = fof_nearest_hsml[target];
    }
  else
    {
      pos = FoFDataGet[target].Pos;
      h = FoFDataGet[target].Hsml;
    }

  index = -1;
  r2max = 1.0e30;

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = FoFDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb_inbox = ngb_treefind_fof_primary(pos, h, target, &startnode, mode, nexport, nsend_local);

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
	      if(r2 < r2max && r2 < h * h)
		{
		  index = j;
		  r2max = r2;
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = FoFDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }


  if(mode == 0)
    {
      if(index >= 0)
	{
	  fof_nearest_distance[target] = sqrt(r2max);
	  MinID[target] = MinID[Head[index]];
	  MinIDTask[target] = MinIDTask[Head[index]];
	}
    }
  else
    {
      if(index >= 0)
	{
	  FoFDataResult[target].Distance = sqrt(r2max);
	  FoFDataResult[target].MinID = MinID[Head[index]];
	  FoFDataResult[target].MinIDTask = MinIDTask[Head[index]];
	}
      else
	FoFDataResult[target].Distance = 2.0e30;
    }
  return 0;
}
#endif


#ifdef BLACK_HOLES

void fof_make_black_holes(void)
{ 
  int i, j, n, ntot;
  int nexport, nimport, sendTask, recvTask, level;
  int *import_indices, *export_indices;
  int self_count;

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

  for(i = 0; i < Ngroups; i++)
    {
      if(Group[i].LenType[1] >= All.MinFoFSizeForNewSeed)
   if(Group[i].LenType[5] == 0)
     {
       if(Group[i].index_maxdens >= 0)
         Send_count[Group[i].task_maxdens]++;
     }
    }

  MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];

      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
   {
     Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
     Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
   }
    }
  import_indices = mymalloc(nimport * sizeof(int));
  export_indices = mymalloc(nexport * sizeof(int));

  for(i=0;i<nimport;i++)
    import_indices[i] = -1;

  for(i=0;i<nexport;i++)
    export_indices[i] = -1;

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

  self_count = 0;
  for(i = 0; i < Ngroups; i++)
    {
   if(Group[i].LenType[1]  >= All.MinFoFSizeForNewSeed && Group[i].LenType[5] == 0)
     {
       if(Group[i].index_maxdens >= 0) {

         if (Group[i].task_maxdens == ThisTask) {
      import_indices[Recv_offset[ThisTask]+self_count++] = Group[i].index_maxdens;
         }
         export_indices[Send_offset[Group[i].task_maxdens] +
              Send_count[Group[i].task_maxdens]++] = Group[i].index_maxdens;
       }
     }
    }

  memcpy(&import_indices[Recv_offset[ThisTask]], &export_indices[Send_offset[ThisTask]], Send_count[ThisTask] * sizeof(int));

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  fflush(stdout);

  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      fflush(stdout);
      if(recvTask < NTask) {
          MPI_Sendrecv(&export_indices[Send_offset[recvTask]],
             Send_count[recvTask] * sizeof(int),
             MPI_BYTE, recvTask, TAG_FOF_E,
             &import_indices[Recv_offset[recvTask]],
             Recv_count[recvTask] * sizeof(int),
             MPI_BYTE, recvTask, TAG_FOF_E, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          }
    }

  MPI_Allreduce(&nimport, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("\nMaking %d new black hole particles\n\n", ntot);
      fprintf(FdBlackHolesSeeds,"%f %d\n",All.Time,ntot);
      fflush(stdout);
      fflush(FdBlackHolesSeeds);
    }


  All.TotBHs += ntot;

  for(n = 0; n < nimport; n++)
    {

      if(P[import_indices[n]].Type != 0)
   endrun(7772);

      P[import_indices[n]].Type    = 5;   /* make it a black hole particle */
      P[import_indices[n]].BH_BirthTime = All.Time;
#ifdef BH_THERMALFEEDBACK
      P[import_indices[n]].BH_Energy = 0;
#endif
      P[import_indices[n]].b1.BH_Density = 1e-10;
      P[import_indices[n]].b2.BH_Entropy = 0.0;
      P[import_indices[n]].b3.BH_SurroundingGasVel[0] = 0.0;
      P[import_indices[n]].b3.BH_SurroundingGasVel[1] = 0.0;
      P[import_indices[n]].b3.BH_SurroundingGasVel[2] = 0.0;
      P[import_indices[n]].BH_Mass = All.SeedBHMassOverGasMass * All.OrigGasMass;
      P[import_indices[n]].BH_Mdot = 0;
      Stars_converted++;
    }

  All.TotN_gas -= ntot;
  myfree(export_indices);
  myfree(import_indices);

  fflush(stdout);
}

#endif




double fof_periodic(double x)
{
  if(x >= 0.5 * All.BoxSize)
    x -= All.BoxSize;
  if(x < -0.5 * All.BoxSize)
    x += All.BoxSize;
  return x;
}


double fof_periodic_wrap(double x)
{
  while(x >= All.BoxSize)
    x -= All.BoxSize;
  while(x < 0)
    x += All.BoxSize;
  return x;
}


int fof_compare_FOF_PList_MinID(const void *a, const void *b)
{
  if(((struct fof_particle_list *) a)->MinID < ((struct fof_particle_list *) b)->MinID)
    return -1;

  if(((struct fof_particle_list *) a)->MinID > ((struct fof_particle_list *) b)->MinID)
    return +1;

  return 0;
}


int fof_compare_FOF_GList_MinID(const void *a, const void *b)
{
  if(((struct fof_group_list *) a)->MinID < ((struct fof_group_list *) b)->MinID)
    return -1;

  if(((struct fof_group_list *) a)->MinID > ((struct fof_group_list *) b)->MinID)
    return +1;

  return 0;
}


int fof_compare_FOF_GList_MinIDTask(const void *a, const void *b)
{
  if(((struct fof_group_list *) a)->MinIDTask < ((struct fof_group_list *) b)->MinIDTask)
    return -1;

  if(((struct fof_group_list *) a)->MinIDTask > ((struct fof_group_list *) b)->MinIDTask)
    return +1;

  return 0;
}


int fof_compare_FOF_GList_LocCountTaskDiffMinID(const void *a, const void *b)
{
  if(((struct fof_group_list *) a)->LocCount > ((struct fof_group_list *) b)->LocCount)
    return -1;

  if(((struct fof_group_list *) a)->LocCount < ((struct fof_group_list *) b)->LocCount)
    return +1;

  if(((struct fof_group_list *) a)->MinID < ((struct fof_group_list *) b)->MinID)
    return -1;

  if(((struct fof_group_list *) a)->MinID > ((struct fof_group_list *) b)->MinID)
    return +1;

  if(abs(((struct fof_group_list *) a)->ExtCount - ((struct fof_group_list *) a)->MinIDTask) <
     abs(((struct fof_group_list *) b)->ExtCount - ((struct fof_group_list *) b)->MinIDTask))
    return -1;

  if(abs(((struct fof_group_list *) a)->ExtCount - ((struct fof_group_list *) a)->MinIDTask) >
     abs(((struct fof_group_list *) b)->ExtCount - ((struct fof_group_list *) b)->MinIDTask))
    return +1;

  return 0;
}


int fof_compare_FOF_GList_ExtCountMinID(const void *a, const void *b)
{
  if(((struct fof_group_list *) a)->ExtCount < ((struct fof_group_list *) b)->ExtCount)
    return -1;

  if(((struct fof_group_list *) a)->ExtCount > ((struct fof_group_list *) b)->ExtCount)
    return +1;

  if(((struct fof_group_list *) a)->MinID < ((struct fof_group_list *) b)->MinID)
    return -1;

  if(((struct fof_group_list *) a)->MinID > ((struct fof_group_list *) b)->MinID)
    return +1;

  return 0;
}


int fof_compare_Group_MinID(const void *a, const void *b)
{
  if(((struct group_properties *) a)->MinID < ((struct group_properties *) b)->MinID)
    return -1;

  if(((struct group_properties *) a)->MinID > ((struct group_properties *) b)->MinID)
    return +1;

  return 0;
}

int fof_compare_Group_GrNr(const void *a, const void *b)
{
  if(((struct group_properties *) a)->GrNr < ((struct group_properties *) b)->GrNr)
    return -1;

  if(((struct group_properties *) a)->GrNr > ((struct group_properties *) b)->GrNr)
    return +1;

  return 0;
}

int fof_compare_Group_MinIDTask(const void *a, const void *b)
{
  if(((struct group_properties *) a)->MinIDTask < ((struct group_properties *) b)->MinIDTask)
    return -1;

  if(((struct group_properties *) a)->MinIDTask > ((struct group_properties *) b)->MinIDTask)
    return +1;

  return 0;
}

int fof_compare_Group_MinIDTask_MinID(const void *a, const void *b)
{
  if(((struct group_properties *) a)->MinIDTask < ((struct group_properties *) b)->MinIDTask)
    return -1;

  if(((struct group_properties *) a)->MinIDTask > ((struct group_properties *) b)->MinIDTask)
    return +1;

  if(((struct group_properties *) a)->MinID < ((struct group_properties *) b)->MinID)
    return -1;

  if(((struct group_properties *) a)->MinID > ((struct group_properties *) b)->MinID)
    return +1;

  return 0;
}


int fof_compare_Group_Len(const void *a, const void *b)
{
  if(((struct group_properties *) a)->Len > ((struct group_properties *) b)->Len)
    return -1;

  if(((struct group_properties *) a)->Len < ((struct group_properties *) b)->Len)
    return +1;

  return 0;
}


int fof_compare_ID_list_GrNrID(const void *a, const void *b)
{
  if(((struct id_list *) a)->GrNr < ((struct id_list *) b)->GrNr)
    return -1;

  if(((struct id_list *) a)->GrNr > ((struct id_list *) b)->GrNr)
    return +1;

  if(((struct id_list *) a)->ID < ((struct id_list *) b)->ID)
    return -1;

  if(((struct id_list *) a)->ID > ((struct id_list *) b)->ID)
    return +1;

  return 0;
}


int fof_compare_ID_list_GrNrTypeID(const void *a, const void *b)
{
  if(((struct id_list *) a)->GrNr < ((struct id_list *) b)->GrNr)
    return -1;

  if(((struct id_list *) a)->GrNr > ((struct id_list *) b)->GrNr)
    return +1;

  if(((struct id_list *) a)->Type < ((struct id_list *) b)->Type)
    return -1;

  if(((struct id_list *) a)->Type > ((struct id_list *) b)->Type)
    return +1;

  if(((struct id_list *) a)->ID < ((struct id_list *) b)->ID)
    return -1;

  if(((struct id_list *) a)->ID > ((struct id_list *) b)->ID)
    return +1;

  return 0;
}


/*
#ifdef BG_SFR
int fof_compare_Star_list_GrNrID(const void *a, const void *b)
{
  if(((struct star_list *) a)->GrNr < ((struct star_list *) b)->GrNr)
    return -1;

  if(((struct star_list *) a)->GrNr > ((struct star_list *) b)->GrNr)
    return +1;

  if(((struct star_list *) a)->ID < ((struct star_list *) b)->ID)
    return -1;

  if(((struct star_list *) a)->ID > ((struct star_list *) b)->ID)
    return +1;

  return 0;
}
#endif
*/

#ifdef HAVE_HDF5
void fof_write_hdf5_double_field(double *data, char *tag, int elems, char *description, int blocknr, hid_t grp)
{
  int i, j;
  double *tmp, *d;
  hid_t dataspace, dataset;
  hsize_t adim[1];

  if(ThisTask == 0)
    printf(" --> %s\n", tag);

  tmp = (double *) mymalloc(sizeof(double) * elems * Ngroups);
  for(i = 0, d = tmp; i < Ngroups; i++)
    {
      for(j = 0; j < elems; j++)
	*d++ = data[j];
      data += (sizeof(struct group_properties) / sizeof(double));
    }

  adim[0] = elems * Ngroups;
  dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(dataspace, 1, adim, NULL);
  dataset = H5Dcreate(grp, tag, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
  if(blocknr >= 0)
    write_attributes_in_hdf5((enum iofields) blocknr, dataset);
  if(*description)
    write_dummy_attributes_in_hdf5(description, dataset);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  myfree(tmp);
}

void fof_write_hdf5_int_field(int *data, char *tag, int elems, char *description, hid_t grp)
{
  int i, j;
  int *tmp, *d;
  hid_t dataspace, dataset;
  hsize_t adim[1];

  if(ThisTask == 0)
    printf(" --> %s\n", tag);

  tmp = (int *) mymalloc(sizeof(int) * elems * Ngroups);
  for(i = 0, d = tmp; i < Ngroups; i++)
    {
      for(j = 0; j < elems; j++)
	*d++ = data[j];
      data += (sizeof(struct group_properties) / sizeof(int));
    }

  adim[0] = elems * Ngroups;
  dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(dataspace, 1, adim, NULL);
  dataset = H5Dcreate(grp, tag, H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
  if(*description)
    write_dummy_attributes_in_hdf5(description, dataset);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  myfree(tmp);
}


void fof_write_hdf5_myIDType_field(MyIDType *data, char *tag, int elems, char *description, hid_t grp)
{
  hid_t dataspace, dataset;
  hsize_t adim[1];

  if(ThisTask == 0)
    printf(" --> %s\n", tag);

  //if elems is Nids
  adim[0] = elems;
  dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(dataspace, 1, adim, NULL);
  dataset = H5Dcreate(grp, tag, H5T_NATIVE_UINT, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  if(*description)
    write_dummy_attributes_in_hdf5(description, dataset);
  H5Dclose(dataset);
  H5Sclose(dataspace);
}


void fof_write_hdf5_star_field(double *data, char *tag, int elems, char *description, int blocknr, hid_t grp)
{
  hid_t dataspace, dataset;
  hsize_t adim[1];

  if(ThisTask == 0)
    printf(" --> %s\n", tag);

  adim[0] = elems;
  dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(dataspace, 1, adim, NULL);

  dataset = H5Dcreate(grp, tag, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  if(blocknr >= 0)
    write_attributes_in_hdf5((enum iofields) blocknr, dataset);
  if(*description)
    write_dummy_attributes_in_hdf5(description, dataset);

  H5Dclose(dataset);
  H5Sclose(dataspace);
}

#endif /* HDF5 routines over */

#endif /* of FOF */
