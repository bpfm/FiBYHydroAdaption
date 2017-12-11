#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#ifdef SUBFIND

#include "fof.h"

#include "allvars.h"
#include "proto.h"
#include "domain.h"
#include "subfind.h"

static struct id_list
{
  long Index;
  MyIDType ID;
  int GrNr;
  int SubNr;
  int Type;
  float BindingEgy;
}
 *ID_list;

static int Nids;

void test_particle_sequence(char *message);

void subfind(int num)
{
  double t0, t1, tstart, tend;

  int i, gr, nlocid, offset, limit, ncount;

  if(ThisTask == 0)
    printf("\nWe now execute a parallel version of SUBFIND.\n");

  tstart = second();

  if(!All.ComovingIntegrationOn)
    {
      if(ThisTask == 0)
	printf("works only for comoving integration.\n");
      endrun(0);
    }

  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

  t0 = second();
  if(ThisTask == 0)
    printf("Tree construction.\n");

  force_treebuild(NumPart, NULL);

  t1 = second();
  if(ThisTask == 0)
    printf("tree build took %g sec\n", timediff(t0, t1));

  /* let's determine the local dark matter densities */
  t0 = second();
  subfind_setup_smoothinglengths();	//Modified to include the baryons for the calculation of hsml - Alan
  subfind_density();		//Modified to include baryons in the evaluation of the subhaloes as well as movement of particles across processors - Alan

  t1 = second();
#ifndef INCL_BARYONS
  if(ThisTask == 0)
    printf("dark matter density() took %g sec\n", timediff(t0, t1));
#else
  if(ThisTask == 0)
    printf("DM + Baryon density() took %g sec\n", timediff(t0, t1));
#endif
  force_treefree();

  /* let's save the densities to a file (for making images) */
  t0 = second();
  subfind_save_densities(num);
  t1 = second();
  if(ThisTask == 0)
    printf("saving densities took %g sec\n", timediff(t0, t1));
  //The density here has been modified to run over all particle types, outputting all particle densities into a new HDF5 file if the correct makefile options are selected - Alan

  /* count how many groups we have that should be done collectively */
  limit = 0.1 * All.TotNumPart;

  for(i = 0, ncount = 0; i < Ngroups; i++)
    if(Group[i].Len >= limit)
      ncount++;
  MPI_Allreduce(&ncount, &Ncollective, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("Number of FOF halos treated with collective SubFind code = %d\n", Ncollective);
      printf("the other %d FOF halos are treated in parallel with serial code\n", TotNgroups - Ncollective);
    }

  if(Ncollective > 0)
    {
      printf("Ncollective = %d!!!! This is wrong, I stop now!!!\n", Ncollective);
      endrun(6492);
    }

  /*  to decide on which task a group should be:
   *  if   GrNr <= Ncollective:  collective groupfinding.
   *  the task where the group info is put is TaskNr = (GrNr - 1) % NTask
   */

  /* now we distribute the particles such that small groups are assigned in
   *  total to certain CPUs, and big groups are left where they are 
   */

  /*  test_particle_sequence("[subfind] before subfind_distribute_particles"); */

  t0 = second();
  subfind_distribute_particles(0);
  t1 = second();
  if(ThisTask == 0)
    printf("distribute_particles()() took %g sec\n", timediff(t0, t1));

  /*  test_particle_sequence("[subfind] after subfind_distribute_particles"); */

  subfind_distribute_groups();	//Until here there is no use of the rearrange_particles subroutine - Alan
  // The distribution of particles and groups is insensitive to particle type - Alan
  qsort(Group, Ngroups, sizeof(struct group_properties), fof_compare_Group_GrNr);

  for(i = 0; i < NumPart; i++)
    if(P[i].GrNr > Ncollective && P[i].GrNr <= TotNgroups)
      if(((P[i].GrNr - 1) % NTask) != ThisTask)
	{
	  printf("i=%d %d task=%d\n", i, P[i].GrNr, ThisTask);
	  endrun(87);
	}

/* #ifdef IGNOREFUZZ */
/*   sizelimit = NumPart * All.GroupSizeLimit; */

/*   printf("sizelimit is %f while NumPart is %d and All.GroupSizeLimit is %f \n",sizelimit,NumPart,All.GroupSizeLimit); */

/*   skipped = 0; */
/* #endif */
  /* lets estimate the maximum number of substructures we need to store on the local CPU */
  for(i = 0, nlocid = 0; i < Ngroups; i++)
    {
      nlocid += Group[i].Len;
/* #ifdef IGNOREFUZZ */
/*       if(Group[i].Len >= sizelimit) */
/* 	skipped++; */
/* #endif */
    }
/* #ifdef IGNOREFUZZ */
/*   printf("There are %d groups which exceed sizelimit \n",skipped); */
/* #endif   */

#ifndef LOW_MEM_USE
  MaxNsubgroups = nlocid / All.DesLinkNgb;	/* this is a quite conservative upper limit */
#else
  MaxNsubgroups = nlocid / (5.0 * All.DesLinkNgb);	// Use less memory by allocating less space for subfind objects - Alan
#endif
  Nsubgroups = 0;
  SubGroup = (struct subgroup_properties *) mymalloc(MaxNsubgroups * sizeof(struct subgroup_properties));
  printf("Size of subgroup_properties %d bytes\n", (int) sizeof(struct subgroup_properties));
  printf("Total memory allocation for subgroup_properties %g MB \n",
	 (MaxNsubgroups * sizeof(struct subgroup_properties)) / (1024.0 * 1024.0));

  initialize_SubGroup(MaxNsubgroups);	//Just sets everything to zero - Alan

  for(i = 0; i < NumPart; i++)
    P[i].SubNr = (1 << 30);	/* default */

  t0 = second();
  for(GrNr = 1; GrNr <= Ncollective; GrNr++)
    subfind_process_group_collectively();
  t1 = second();
  if(ThisTask == 0)
    printf("processing of collective halos took %g sec\n", timediff(t0, t1));

  Ind = (struct my_serial_index *) mymalloc(NumPart * sizeof(struct my_serial_index));

  for(i = 0; i < NumPart; i++)
    {
      Ind[i].Index = i;
      Ind[i].Density = P[i].u.DM_Density;
      Ind[i].GrNr = P[i].GrNr;
    }

  t0 = second();
  qsort(Ind, NumPart, sizeof(struct my_serial_index), subfind_compare_Ind_GrNr_DM_Density);
  t1 = second();
  if(ThisTask == 0)
    printf("sort of local particles()() took %g sec\n", timediff(t0, t1));

  Rank = (int *) mymalloc(NumPart * sizeof(int));

  for(i = 0; i < NumPart; i++)
    Rank[Ind[i].Index] = i;

  /* let's count how many local particles we have in small groups */
  for(i = 0, nlocid = 0; i < NumPart; i++)
    if(P[Ind[i].Index].GrNr > Ncollective && P[i].GrNr <= Ngroups)	/* particle is in small group */
      nlocid++;

  if(ThisTask == 0)
    printf("contructing tree for serial subfind of local groups\n");

  subfind_loctree_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

  if(ThisTask == 0)
    printf("Start to do local groups with serial subfind algorithm\n");

  t0 = second();

  int iter, skip;

  double t2, t3;

  /* we now apply a serial version of subfind to the local groups */
  for(gr = 0, offset = 0; gr < Ngroups; gr++)
    {
      if(Group[gr].GrNr > Ncollective)
	{
	  t2 = second();
	  if(((Group[gr].GrNr - 1) % NTask) == ThisTask)
	    offset = subfind_process_group_serial(gr, offset, &iter, &skip);
	  t3 = second();
	  printf("%d\t%d\t%f\t%d\t%d\n", ThisTask, Group[gr].GrNr, timediff(t2, t3), iter, Group[gr].Len);

	  fflush(stdout);
	}
    }

/* #ifdef IGNOREFUZZ */
/*   if(skip != skipped) */
/*     { */
/*       printf("Predicted that %d haloes would be skipped, instead %d were missed \n",skipped,skip); */
/*       endrun(816); */
/*     } */
/* #endif */

  t1 = second();
  if(ThisTask == 0)
    printf("processing of local groups took %g sec\n", timediff(t0, t1));

  t0 = second();

  subfind_loctree_treefree();

  GrNr = -1;			/* to ensure that domain decomposition acts normally again */

  /* now determine the remaining spherical overdensity values for the non-local groups */
  domain_free_trick();

  domain_Decomposition();

  force_treebuild(NumPart, NULL);

  t1 = second();
  if(ThisTask == 0)
    printf("freeing/building trees and domains took %g sec\n", timediff(t0, t1));

  /* compute spherical overdensities for FOF groups */
  t0 = second();

  subfind_overdensity();

  t1 = second();
  if(ThisTask == 0)
    printf("determining spherical overdensity masses took %g sec\n", timediff(t0, t1));

#ifdef CONTAMINATION
  /* determine which halos are contaminated by boundary particles */
  t0 = second();

  subfind_contamination();

  t1 = second();
  if(ThisTask == 0)
    printf("determining contamination of halos took %g sec\n", timediff(t0, t1));
#endif

  force_treefree();
  domain_free();

  domain_allocate_trick();

  /* now assemble final output */
  subfind_save_final(num);

  tend = second();

  if(ThisTask == 0)
    printf("\nFinished with SUBFIND.  (total time=%g sec)\n\n", timediff(tstart, tend));

/* #ifdef IGNOREFUZZ */
/*   if(skipped > 0) */
/*     myfree(SkipGrUnbind); */
/* #endif */
  myfree(Rank);
  myfree(Ind);

  myfree(SubGroup);

}

void subfind_save_final(int num)
{
  int i, j, totsubs, masterTask, groupTask, nprocgroup;

  char buf[1000];

  double t0, t1;

  /* prepare list of ids with assigned group numbers */

  parallel_sort(Group, Ngroups, sizeof(struct group_properties), fof_compare_Group_GrNr);
  parallel_sort(SubGroup, Nsubgroups, sizeof(struct subgroup_properties),
		subfind_compare_SubGroup_GrNr_SubNr);

  ID_list = mymalloc(sizeof(struct id_list) * NumPart);

  MPI_Allreduce(&Nsubgroups, &TotNsubgroups, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(i = 0, Nids = 0; i < NumPart; i++)
    {
      if(P[i].GrNr <= TotNgroups)
	{
	  ID_list[Nids].Index = i;
	  ID_list[Nids].GrNr = P[i].GrNr;
	  ID_list[Nids].SubNr = P[i].SubNr;
	  ID_list[Nids].Type = P[i].Type;
	  ID_list[Nids].BindingEgy = P[i].v.DM_BindingEnergy;
	  ID_list[Nids].ID = P[i].ID;

	  Nids++;
	}
    }

  //Nids is calculated here, as well as the ID_list struct - Alan
  // Now sorts GrNr, SubNr, Type, BE
  parallel_sort(ID_list, Nids, sizeof(struct id_list), subfind_compare_ID_list);
  /* fill in the FirstSub-values */
  for(i = 0, totsubs = 0; i < Ngroups; i++)
    {
      if(i > 0)
	Group[i].FirstSub = Group[i - 1].FirstSub + Group[i - 1].Nsubs;
      else
	Group[i].FirstSub = 0;
      totsubs += Group[i].Nsubs;
    }

  MPI_Allgather(&totsubs, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  for(j = 1, Send_offset[0] = 0; j < NTask; j++)
    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];

  for(i = 0; i < Ngroups; i++)
    Group[i].FirstSub += Send_offset[ThisTask];

  MPI_Allgather(&Nids, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  for(j = 1, Send_offset[0] = 0; j < NTask; j++)
    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];

  if(ThisTask == 0)
    {
#ifdef INCL_DMONLY
      sprintf(buf, "%s/groups_%03d_DMONLY", All.OutputDir, num);
#else
      sprintf(buf, "%s/groups_%03d", All.OutputDir, num);
#endif
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);

  if(NTask < All.NumFilesWrittenInParallel)
    {
      printf
	("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");
      endrun(241931);
    }

  t0 = second();

  nprocgroup = NTask / All.NumFilesWrittenInParallel;
  if((NTask % All.NumFilesWrittenInParallel))
    nprocgroup++;
  masterTask = (ThisTask / nprocgroup) * nprocgroup;
  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))	/* ok, it's this processor's turn */
	{
	  if(All.SnapFormat != 3)
	    subfind_save_local_catalogue(num);
	  else
	    subfind_save_local_catalogue_hdf5(num);
	}
      MPI_Barrier(MPI_COMM_WORLD);	/* wait inside the group */
    }

  t1 = second();

  if(ThisTask == 0)
    {
      printf("Subgroup catalogues saved. took = %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }

  myfree(ID_list);
}


void subfind_save_local_catalogue(int num)
{
  FILE *fd;

  char buf[500], fname[500];

  MyDouble *mass, *pos, *vel;

  float *spin;

  int i, j, *len;

  MyIDType *ids;

#ifdef INCL_DMONLY
  sprintf(fname, "%s/groups_%03d_DMONLY/%s_%03d.%d", All.OutputDir, num, "subhalo_tab", num, ThisTask);
#else
  sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_tab", num, ThisTask);
#endif
  strcpy(buf, fname);
  if(!(fd = fopen(buf, "w")))
    {
      printf("can't open file `%s`\n", buf);
      endrun(1183);
    }

  my_fwrite(&Ngroups, sizeof(int), 1, fd);
  my_fwrite(&TotNgroups, sizeof(int), 1, fd);
  my_fwrite(&Nids, sizeof(int), 1, fd);
  my_fwrite(&TotNids, sizeof(long long), 1, fd);
  my_fwrite(&NTask, sizeof(int), 1, fd);
  my_fwrite(&Nsubgroups, sizeof(int), 1, fd);
  my_fwrite(&TotNsubgroups, sizeof(int), 1, fd);

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

  /* location (potential minimum) */
  pos = mymalloc(Ngroups * 3 * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    for(j = 0; j < 3; j++)
      pos[i * 3 + j] = Group[i].Pos[j];
  my_fwrite(pos, Ngroups, 3 * sizeof(float), fd);
  myfree(pos);

  /* M_Mean200 */
  mass = mymalloc(Ngroups * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    mass[i] = Group[i].M_Mean200;
  my_fwrite(mass, Ngroups, sizeof(float), fd);
  myfree(mass);

  /* R_Mean200 */
  mass = mymalloc(Ngroups * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    mass[i] = Group[i].R_Mean200;
  my_fwrite(mass, Ngroups, sizeof(float), fd);
  myfree(mass);

  /* M_Crit200 */
  mass = mymalloc(Ngroups * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    mass[i] = Group[i].M_Crit200;
  my_fwrite(mass, Ngroups, sizeof(float), fd);
  myfree(mass);

  /* R_Crit200 */
  mass = mymalloc(Ngroups * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    mass[i] = Group[i].R_Crit200;
  my_fwrite(mass, Ngroups, sizeof(float), fd);
  myfree(mass);

  /* M_TopHat200 */
  mass = mymalloc(Ngroups * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    mass[i] = Group[i].M_TopHat200;
  my_fwrite(mass, Ngroups, sizeof(float), fd);
  myfree(mass);

  /* R_TopHat200 */
  mass = mymalloc(Ngroups * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    mass[i] = Group[i].R_TopHat200;
  my_fwrite(mass, Ngroups, sizeof(float), fd);
  myfree(mass);

#ifdef CONTAMINATION
  /* contamination particle count */
  len = mymalloc(Ngroups * sizeof(int));
  for(i = 0; i < Ngroups; i++)
    len[i] = Group[i].ContaminationLen;
  my_fwrite(len, Ngroups, sizeof(int), fd);
  myfree(len);

  /* contamination mass */
  mass = mymalloc(Ngroups * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    mass[i] = Group[i].ContaminationMass;
  my_fwrite(mass, Ngroups, sizeof(float), fd);
  myfree(mass);
#endif

  /* number of substructures in FOF group  */
  len = mymalloc(Ngroups * sizeof(int));
  for(i = 0; i < Ngroups; i++)
    len[i] = Group[i].Nsubs;
  my_fwrite(len, Ngroups, sizeof(int), fd);
  myfree(len);

  /* first substructure in FOF group  */
  len = mymalloc(Ngroups * sizeof(int));
  for(i = 0; i < Ngroups; i++)
    len[i] = Group[i].FirstSub;
  my_fwrite(len, Ngroups, sizeof(int), fd);
  myfree(len);

  /* ------------------------------ */

  /* Len of substructure  */
  len = mymalloc(Nsubgroups * sizeof(int));
  for(i = 0; i < Nsubgroups; i++)
    len[i] = SubGroup[i].Len;
  my_fwrite(len, Nsubgroups, sizeof(int), fd);
  myfree(len);

  /* Len per type of substructure  */
  len = mymalloc(Nsubgroups * 6 * sizeof(int));
  for(i = 0; i < Nsubgroups; i++)
    for(j = 0; j < 6; j++)
      len[i * 6 + j] = SubGroup[i].LenType[j];
  my_fwrite(len, Nsubgroups, 6 * sizeof(int), fd);
  myfree(len);

  /* offset of substructure  */
  len = mymalloc(Nsubgroups * sizeof(int));
  for(i = 0; i < Nsubgroups; i++)
    len[i] = SubGroup[i].Offset;
  my_fwrite(len, Nsubgroups, sizeof(int), fd);
  myfree(len);

  /* offset per type of substructure  */
  len = mymalloc(Nsubgroups * 6 * sizeof(int));
  for(i = 0; i < Nsubgroups; i++)
    for(j = 0; j < 6; j++)
      len[i * 6 + j] = SubGroup[i].OffsetType[j];
  my_fwrite(len, Nsubgroups, 6 * sizeof(int), fd);
  myfree(len);

  /* parent of substructure  */
  len = mymalloc(Nsubgroups * sizeof(int));
  for(i = 0; i < Nsubgroups; i++)
    len[i] = SubGroup[i].SubParent;
  my_fwrite(len, Nsubgroups, sizeof(int), fd);
  myfree(len);

  /* Mass of substructure  */
  mass = mymalloc(Nsubgroups * sizeof(float));
  for(i = 0; i < Nsubgroups; i++)
    mass[i] = SubGroup[i].Mass;
  my_fwrite(mass, Nsubgroups, sizeof(float), fd);
  myfree(mass);

  /* Mass of substructure  */
  mass = mymalloc(Nsubgroups * 6 * sizeof(float));
  for(i = 0; i < Nsubgroups; i++)
    for(j = 0; j < 6; j++)
      mass[i * 6 + j] = SubGroup[i].MassType[j];
  my_fwrite(mass, Nsubgroups, 6 * sizeof(float), fd);
  myfree(mass);

  /* Pos of substructure  */
  pos = mymalloc(Nsubgroups * 3 * sizeof(float));
  for(i = 0; i < Nsubgroups; i++)
    for(j = 0; j < 3; j++)
      pos[i * 3 + j] = SubGroup[i].Pos[j];
  my_fwrite(pos, Nsubgroups, 3 * sizeof(float), fd);
  myfree(pos);

  /* Vel of substructure  */
  vel = mymalloc(Nsubgroups * 3 * sizeof(float));
  for(i = 0; i < Nsubgroups; i++)
    for(j = 0; j < 3; j++)
      vel[i * 3 + j] = SubGroup[i].Vel[j];
  my_fwrite(vel, Nsubgroups, 3 * sizeof(float), fd);
  myfree(vel);

  /* Center of mass of substructure  */
  pos = mymalloc(Nsubgroups * 3 * sizeof(float));
  for(i = 0; i < Nsubgroups; i++)
    for(j = 0; j < 3; j++)
      pos[i * 3 + j] = SubGroup[i].CM[j];
  my_fwrite(pos, Nsubgroups, 3 * sizeof(float), fd);
  myfree(pos);

  /* Spin of substructure  */
  spin = mymalloc(Nsubgroups * 3 * sizeof(float));
  for(i = 0; i < Nsubgroups; i++)
    for(j = 0; j < 3; j++)
      spin[i * 3 + j] = SubGroup[i].Spin[j];
  my_fwrite(spin, Nsubgroups, 3 * sizeof(float), fd);
  myfree(spin);

  /* velocity dispersion  */
  mass = mymalloc(Nsubgroups * sizeof(float));
  for(i = 0; i < Nsubgroups; i++)
    mass[i] = SubGroup[i].SubVelDisp;
  my_fwrite(mass, Nsubgroups, sizeof(float), fd);
  myfree(mass);

  /* stellar velocity dispersion  */
  mass = mymalloc(Nsubgroups * sizeof(float));
  for(i = 0; i < Nsubgroups; i++)
    mass[i] = SubGroup[i].SubStellarVelDisp;
  my_fwrite(mass, Nsubgroups, sizeof(float), fd);
  myfree(mass);

  /* stellar velocity dispersion within projected half mass  */
  mass = mymalloc(Nsubgroups * sizeof(float));
  for(i = 0; i < Nsubgroups; i++)
    mass[i] = SubGroup[i].SubStellarVelDispHalfProj;
  my_fwrite(mass, Nsubgroups, sizeof(float), fd);
  myfree(mass);

  /* maximum circular velocity  */
  mass = mymalloc(Nsubgroups * sizeof(float));
  for(i = 0; i < Nsubgroups; i++)
    mass[i] = SubGroup[i].SubVmax;
  my_fwrite(mass, Nsubgroups, sizeof(float), fd);
  myfree(mass);

  /* radius of maximum circular velocity  */
  mass = mymalloc(Nsubgroups * sizeof(float));
  for(i = 0; i < Nsubgroups; i++)
    mass[i] = SubGroup[i].SubVmaxRad;
  my_fwrite(mass, Nsubgroups, sizeof(float), fd);
  myfree(mass);

  /* radius of half the mass  */
  mass = mymalloc(Nsubgroups * sizeof(float));
  for(i = 0; i < Nsubgroups; i++)
    for(j = 0; j < 6; j++)
      mass[i * 6 + j] = SubGroup[i].SubHalfMass[j];
  my_fwrite(mass, Nsubgroups, 6 * sizeof(float), fd);
  myfree(mass);

  /* projected radius of half the mass  */
  mass = mymalloc(Nsubgroups * sizeof(float));
  for(i = 0; i < Nsubgroups; i++)
    for(j = 0; j < 6; j++)
      mass[i * 6 + j] = SubGroup[i].SubHalfMassProj[j];
  my_fwrite(mass, Nsubgroups, 6 * sizeof(float), fd);
  myfree(mass);

  /* ID of most bound particle  */
  len = mymalloc(Nsubgroups * sizeof(MyIDType));
  for(i = 0; i < Nsubgroups; i++)
    len[i] = SubGroup[i].SubMostBoundID;
  my_fwrite(len, Nsubgroups, sizeof(MyIDType), fd);
  myfree(len);

  /* GrNr of substructure  */
  len = mymalloc(Nsubgroups * sizeof(int));
  for(i = 0; i < Nsubgroups; i++)
    len[i] = SubGroup[i].GrNr;
  my_fwrite(len, Nsubgroups, sizeof(int), fd);
  myfree(len);

  fclose(fd);

  // Output particle IDs
  ids = (MyIDType *) ID_list;

  for(i = 0; i < Nids; i++)
    ids[i] = ID_list[i].ID;
#ifdef INCL_DMONLY
  sprintf(buf, "%s/groups_%03d_DMONLY/%s_%03d.%d", All.OutputDir, num, "subhalo_ids", num, ThisTask);
#else
  sprintf(buf, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_ids", num, ThisTask);
#endif
  if(!(fd = fopen(buf, "w")))
    {
      printf("can't open file `%s`\n", buf);
      endrun(1184);
    }

  my_fwrite(&Ngroups, sizeof(int), 1, fd);
  my_fwrite(&TotNgroups, sizeof(int), 1, fd);
  my_fwrite(&Nids, sizeof(int), 1, fd);
  my_fwrite(&TotNids, sizeof(long long), 1, fd);
  my_fwrite(&NTask, sizeof(int), 1, fd);
  my_fwrite(&Send_offset[ThisTask], sizeof(int), 1, fd);
  my_fwrite(ids, sizeof(int), Nids, fd);

  // Output particle masses
  mass = (MyDouble *) mymalloc(Nids * sizeof(MyDouble));
  for(i = 0; i < Nids; i++)
    mass[i] = P[ID_list[i].Index].Mass;
#ifdef INCL_DMONLY
  sprintf(buf, "%s/groups_%03d_DMONLY/%s_%03d.%d", All.OutputDir, num, "subhalo_mass", num, ThisTask);
#else
  sprintf(buf, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_mass", num, ThisTask);
#endif
  if(!(fd = fopen(buf, "w")))
    {
      printf("can't open file `%s`\n", buf);
      endrun(1184);
    }
  my_fwrite(&Ngroups, sizeof(int), 1, fd);
  my_fwrite(&TotNgroups, sizeof(int), 1, fd);
  my_fwrite(&Nids, sizeof(int), 1, fd);
  my_fwrite(&TotNids, sizeof(long long), 1, fd);
  my_fwrite(&NTask, sizeof(int), 1, fd);
  my_fwrite(&Send_offset[ThisTask], sizeof(int), 1, fd);
  my_fwrite(mass, sizeof(MyDouble), Nids, fd);
  myfree(mass);

  // Output particle positions
  pos = (MyDouble *) mymalloc(Nids * sizeof(MyDouble));
  for(i = 0; i < Nids; i++)
    for(j = 0; j < 3; j++)
      pos[3 * i + j] = P[ID_list[i].Index].Pos[j];
#ifdef INCL_DMONLY
  sprintf(buf, "%s/groups_%03d_DMONLY/%s_%03d.%d", All.OutputDir, num, "subhalo_pos", num, ThisTask);
#else
  sprintf(buf, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_pos", num, ThisTask);
#endif
  if(!(fd = fopen(buf, "w")))
    {
      printf("can't open file `%s`\n", buf);
      endrun(1184);
    }
  my_fwrite(&Ngroups, sizeof(int), 1, fd);
  my_fwrite(&TotNgroups, sizeof(int), 1, fd);
  my_fwrite(&Nids, sizeof(int), 1, fd);
  my_fwrite(&TotNids, sizeof(long long), 1, fd);
  my_fwrite(&NTask, sizeof(int), 1, fd);
  my_fwrite(&Send_offset[ThisTask], sizeof(int), 1, fd);
  my_fwrite(pos, sizeof(MyDouble), Nids, fd);
  myfree(pos);

  // Output particle velocities
  vel = (MyDouble *) mymalloc(Nids * sizeof(MyDouble));
  for(i = 0; i < Nids; i++)
    for(j = 0; j < 3; j++)
      vel[3 * i + j] = P[ID_list[i].Index].Vel[j];
#ifdef INCL_DMONLY
  sprintf(buf, "%s/groups_%03d_DMONLY/%s_%03d.%d", All.OutputDir, num, "subhalo_vel", num, ThisTask);
#else
  sprintf(buf, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_vel", num, ThisTask);
#endif
  if(!(fd = fopen(buf, "w")))
    {
      printf("can't open file `%s`\n", buf);
      endrun(1184);
    }
  my_fwrite(&Ngroups, sizeof(int), 1, fd);
  my_fwrite(&TotNgroups, sizeof(int), 1, fd);
  my_fwrite(&Nids, sizeof(int), 1, fd);
  my_fwrite(&TotNids, sizeof(long long), 1, fd);
  my_fwrite(&NTask, sizeof(int), 1, fd);
  my_fwrite(&Send_offset[ThisTask], sizeof(int), 1, fd);
  my_fwrite(vel, sizeof(MyDouble), Nids, fd);
  myfree(vel);

  fclose(fd);
}

#ifdef HAVE_HDF5
void subfind_save_local_catalogue_hdf5(int num)
{
  char buf[500], fname[500];

  int i, j;

  MyDouble *pos, *vel, *mass;

  MyIDType *ids;

  hid_t file_id = 0;

  hid_t headergrp = 0, paramgrp = 0;

  hid_t constgrp = 0, unitsgrp = 0;

  hid_t groupgrp = 0;

  hid_t dataspace = 0, attribute = 0;

#ifdef INCL_DMONLY
  sprintf(fname, "%s/groups_%03d_DMONLY/%s_%03d.%d", All.OutputDir, num, "subhalo", num, ThisTask);
#else
  sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo", num, ThisTask);
#endif
  sprintf(buf, "%s.hdf5", fname);

  file_id = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  if(file_id < 0)
    {
      printf("can't open file `%s`\n", buf);
      exit(0);
    }

  fflush(stdout);

  /* write usual headers */
  headergrp = H5Gcreate(file_id, "/Header", 0);
  paramgrp = H5Gcreate(file_id, "/Parameters", 0);
  unitsgrp = H5Gcreate(file_id, "/Units", 0);
  constgrp = H5Gcreate(file_id, "/Constants", 0);

  /* determine global and local particle numbers */
  int ntot_type[6];		//, *temp;

  long long ntot_type_all[6];

  for(i = 0; i < 6; i++)
    ntot_type[i] = 0;

  for(i = 0; i < NumPart; i++)
    ntot_type[P[i].Type]++;

  for(i = 0; i < 6; i++)
    {
      ntot_type_all[i] = ntot_type[i];
    }

  /* Note: this function will not be called by all processors, hence a collective
   * communication will lead to deadlock
   */

  /*  
     temp = (int *) mymalloc(NTask * 6 * sizeof(int));
     MPI_Allgather(ntot_type, 6, MPI_INT, temp, 6, MPI_INT, MPI_COMM_WORLD);
     for(i = 0; i < 6; i++)
     {
     ntot_type_all[i] = 0;
     for(j = 0; j < NTask; j++)
     ntot_type_all[i] += temp[j * 6 + i];
     }
     myfree(temp);
   */


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
  header.num_files = All.NumFilesPerSnapshot;
  header.BoxSize = All.BoxSize;
  header.Omega0 = All.Omega0;
  header.OmegaLambda = All.OmegaLambda;
  header.HubbleParam = All.HubbleParam;

  write_header_attributes_in_hdf5(headergrp);
  write_parameters_attributes_in_hdf5(paramgrp);
  write_units_attributes_in_hdf5(unitsgrp);
  write_constants_attributes_in_hdf5(constgrp);

  H5Gclose(headergrp);
  H5Gclose(paramgrp);
  H5Gclose(unitsgrp);
  H5Gclose(constgrp);

  /* specific for group files */
  groupgrp = H5Gcreate(file_id, "/FOF", 0);
  dataspace = H5Screate(H5S_SCALAR);

  /* Ngroups */
  attribute = H5Acreate(groupgrp, "Number_of_FOF_Groups", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_INT, &Ngroups);
  H5Aclose(attribute);

  /* TotNgroups */
  attribute = H5Acreate(groupgrp, "Total_Number_of_FOF_Groups", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_INT, &TotNgroups);
  H5Aclose(attribute);

  /* Nids */
  attribute = H5Acreate(groupgrp, "Number_of_IDs", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_INT, &Nids);
  H5Aclose(attribute);

  /* Nids */
  attribute = H5Acreate(groupgrp, "Total_Number_of_IDs", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_INT, &TotNids);
  H5Aclose(attribute);

  /* NTask */
  attribute = H5Acreate(groupgrp, "NTask", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_INT, &NTask);
  H5Aclose(attribute);

  /* Nsubgroups */
  attribute = H5Acreate(groupgrp, "Number_of_groups", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_INT, &Nsubgroups);
  H5Aclose(attribute);

  /* TotNsubgroups */
  attribute = H5Acreate(groupgrp, "Total_Number_of_groups", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_INT, &TotNsubgroups);
  H5Aclose(attribute);

  H5Sclose(dataspace);

  /* only if any groups have actually been found */
  if(Ngroups > 0)
    {
      fof_write_hdf5_int_field(&Group[0].Len, "FOF_Length", 1, "Number of particles in each group", groupgrp);

      fof_write_hdf5_int_field(&Group[0].Offset, "FOF_Offset", 1, "Particle number offset of each group",
			       groupgrp);

      fof_write_hdf5_double_field(&Group[0].Mass, "FOF_Mass", 1, "", IO_MASS, groupgrp);

      fof_write_hdf5_double_field(&Group[0].Pos[0], "FOF_Position", 3, "", IO_POS, groupgrp);

      fof_write_hdf5_double_field(&Group[0].M_Mean200, "Halo_M_Mean200", 1, "", IO_M_MEAN, groupgrp);

      fof_write_hdf5_double_field(&Group[0].R_Mean200, "Halo_R_Mean200", 1, "", IO_R_MEAN, groupgrp);

      fof_write_hdf5_double_field(&Group[0].M_Crit200, "Halo_M_Crit200", 1, "", IO_M_CRIT, groupgrp);

      fof_write_hdf5_double_field(&Group[0].R_Crit200, "Halo_R_Crit200", 1, "", IO_R_CRIT, groupgrp);

      fof_write_hdf5_double_field(&Group[0].M_TopHat200, "Halo_M_TopHat200", 1, "", IO_M_TOPHAT, groupgrp);

      fof_write_hdf5_double_field(&Group[0].R_TopHat200, "Halo_R_TopHat200", 1, "", IO_R_TOPHAT, groupgrp);

#ifdef CONTAMINATION
      fof_write_hdf5_int_field(&Group[0].ContaminationLen, "ContaminationLen", 1,
			       "Number of contamination particles within each FOF group", groupgrp);

      fof_write_hdf5_double_field(&Group[0].ContaminationMass, "ContaminationMass", 1, "", IO_MASS, groupgrp);
#endif

      fof_write_hdf5_int_field(&Group[0].Nsubs, "NsubPerHalo", 1, "Number of subhaloes within each FOF group",
			       groupgrp);

      fof_write_hdf5_int_field(&Group[0].FirstSub, "FirstSubOfHalo", 1,
			       "Subhalo list number of the first substructure within each FOF group",
			       groupgrp);
    }

  if(Nsubgroups > 0)
    {
      subfind_write_hdf5_int_field(&SubGroup[0].Len, "Length", 1, "Number of particles in each subhalo",
				   groupgrp);

      subfind_write_hdf5_int_field(&SubGroup[0].LenType[0], "LengthType", 6,
				   "Number of particles per type in each subhalo", groupgrp);

      subfind_write_hdf5_int_field(&SubGroup[0].Offset, "Offset", 1, "Particle number offset of each subhalo",
				   groupgrp);

      subfind_write_hdf5_int_field(&SubGroup[0].OffsetType[0], "OffsetType", 6,
				   "Particle number offset per type of each subhalo", groupgrp);

      subfind_write_hdf5_int_field(&SubGroup[0].SubParent, "SubParentHalo", 1,
				   "Subhalo parent for subsequent substructure, zero if no further levels",
				   groupgrp);

      subfind_write_hdf5_int_field(&SubGroup[0].GrNr, "GrNr", 1, "Group Number of each subhalo", groupgrp);

      subfind_write_hdf5_MyIDType_field(&SubGroup[0].SubMostBoundID, "SubMostBoundID", 1,
					"Most bound particle of each subhalo", groupgrp);

      subfind_write_hdf5_double_field(&SubGroup[0].Mass, "Mass", 1, "", IO_SUBFIND_MASS, groupgrp);

      subfind_write_hdf5_double_field(&SubGroup[0].MassType[0], "MassType", 6, "", IO_SUBFIND_MASSTYPE,
				      groupgrp);

      subfind_write_hdf5_double_field(&SubGroup[0].Pos[0], "Position", 3, "", IO_SUBFIND_POS, groupgrp);

      subfind_write_hdf5_double_field(&SubGroup[0].Vel[0], "CenterOfMassVelocity", 3, "", IO_SUBFIND_CMVEL,
				      groupgrp);

      subfind_write_hdf5_double_field(&SubGroup[0].CM[0], "CenterOfMass", 3, "", IO_SUBFIND_CMPOS, groupgrp);

      subfind_write_hdf5_double_field(&SubGroup[0].SubVelDisp, "SubVelDisp", 1, "", IO_SUBFIND_VELDISP,
				      groupgrp);

      subfind_write_hdf5_double_field(&SubGroup[0].SubStellarVelDisp, "SubStellarVelDisp", 1, "",
				      IO_SUBFIND_STELLARVELDISP, groupgrp);

      subfind_write_hdf5_double_field(&SubGroup[0].SubStellarVelDispHalfProj, "SubStellarVelDispHalfProj", 1,
				      "", IO_SUBFIND_STELLARVELDISPHALFPROJ, groupgrp);

      subfind_write_hdf5_double_field(&SubGroup[0].SubHalfMass[0], "SubHalfMass", 6, "", IO_SUBFIND_HALFMASS,
				      groupgrp);

      subfind_write_hdf5_double_field(&SubGroup[0].SubHalfMassProj[0], "SubHalfMassProj", 6, "",
				      IO_SUBFIND_HALFMASSPROJ, groupgrp);

      subfind_write_hdf5_double_field(&SubGroup[0].SubVmax, "SubVmax", 1, "", IO_SUBFIND_VMAX, groupgrp);

      subfind_write_hdf5_double_field(&SubGroup[0].SubVmaxRad, "SubVmaxRad", 1, "", IO_SUBFIND_VMAXRAD,
				      groupgrp);

      subfind_write_hdf5_double_field(&SubGroup[0].Spin[0], "SubSpin", 3, "", IO_SUBFIND_SPIN, groupgrp);
    }

  if(Nids)
    {
      // Write the IDs
      ids = (MyIDType *) mymalloc(Nids * sizeof(MyIDType));
      for(i = 0; i < Nids; i++)
	ids[i] = ID_list[i].ID;
      subfind_write_hdf5_MyIDType_particle(&ids[0], "ParticleIDs", 1, "Unique particle identifier", groupgrp);
      myfree(ids);

      // Write the masses of the particles
      mass = (MyDouble *) mymalloc(Nids * sizeof(MyDouble));
      for(i = 0; i < Nids; i++)
	mass[i] = P[ID_list[i].Index].Mass;
      subfind_write_hdf5_MyDouble_particle(&mass[0], "ParticleMasses", 1, "", IO_SUBFIND_PARTMASS, groupgrp);
      myfree(mass);

      // Write the positions of the particles
      pos = (MyDouble *) mymalloc(Nids * 3 * sizeof(MyDouble));
      for(i = 0; i < Nids; i++)
	for(j = 0; j < 3; j++)
	  pos[i * 3 + j] = P[ID_list[i].Index].Pos[j];
      subfind_write_hdf5_MyDouble_particle(&pos[0], "ParticlePositions", 3, "", IO_SUBFIND_PARTPOS, groupgrp);
      myfree(pos);

      // Write the velocities of the particles
      vel = (MyDouble *) mymalloc(Nids * 3 * sizeof(MyDouble));
      for(i = 0; i < Nids; i++)
	for(j = 0; j < 3; j++)
	  vel[i * 3 + j] = P[ID_list[i].Index].Vel[j];
      subfind_write_hdf5_MyDouble_particle(&vel[0], "ParticleVelocities", 3, "", IO_SUBFIND_PARTVEL,
					   groupgrp);
      myfree(vel);
    }

  H5Gclose(groupgrp);
  H5Fclose(file_id);
}

void subfind_write_hdf5_double_field(double *data, char *tag, int elems, char *description, int blocknr,
				     hid_t grp)
{
  int i, j;

  double *tmp, *d;

  hid_t dataspace, dataset;

  hsize_t adim[1];

  // Previously this was written for MaxNsubgroups but this is wasteful in terms of outputted memory
  //  tmp = (double *) mymalloc(sizeof(double) * elems * MaxNsubgroups);
  // Much better to write in terms of Nsubgroups
  tmp = (double *) mymalloc(sizeof(double) * elems * Nsubgroups);
  for(i = 0, d = tmp; i < Nsubgroups; i++)
    {
      for(j = 0; j < elems; j++)
	*d++ = data[j];
      data += (sizeof(struct subgroup_properties) / sizeof(double));
    }

  adim[0] = elems * Nsubgroups;
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

void subfind_write_hdf5_int_field(int *data, char *tag, int elems, char *description, hid_t grp)
{
  int i, j;

  int *tmp, *d;

  hid_t dataspace, dataset;

  hsize_t adim[1];

  // Previously this was written for MaxNsubgroups but this is wasteful in terms of outputted memory
  //  tmp = (double *) mymalloc(sizeof(double) * elems * MaxNsubgroups);
  // Much better to write in terms of Nsubgroups
  tmp = (int *) mymalloc(sizeof(int) * elems * Nsubgroups);
  for(i = 0, d = tmp; i < Nsubgroups; i++)
    {
      for(j = 0; j < elems; j++)
	*d++ = data[j];
      data += (sizeof(struct subgroup_properties) / sizeof(int));
    }

  adim[0] = elems * Nsubgroups;
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

void subfind_write_hdf5_MyIDType_field(MyIDType * data, char *tag, int elems, char *description, hid_t grp)
{
  int i, j;

  MyIDType *tmp, *d;

  hid_t dataspace, dataset;

  hsize_t adim[1];

  tmp = (MyIDType *) mymalloc(sizeof(MyIDType) * elems * Nsubgroups);
  for(i = 0, d = tmp; i < Nsubgroups; i++)
    {
      for(j = 0; j < elems; j++)
	*d++ = data[j];
      data += (sizeof(struct subgroup_properties) / sizeof(MyIDType));
    }

  adim[0] = elems * Nsubgroups;
  dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(dataspace, 1, adim, NULL);

#ifndef LONGIDS
  dataset = H5Dcreate(grp, tag, H5T_NATIVE_UINT, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
#else
  dataset = H5Dcreate(grp, tag, H5T_NATIVE_ULLONG, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
#endif
  if(*description)
    write_dummy_attributes_in_hdf5(description, dataset);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  myfree(tmp);
}

void subfind_write_hdf5_MyIDType_particle(MyIDType * data, char *tag, int elems, char *description, hid_t grp)
{
  hid_t dataspace, dataset;

  hsize_t adim[1];

  adim[0] = elems * Nids;
  dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(dataspace, 1, adim, NULL);

#ifndef LONGIDS
  dataset = H5Dcreate(grp, tag, H5T_NATIVE_UINT, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
#else
  dataset = H5Dcreate(grp, tag, H5T_NATIVE_ULLONG, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
#endif
  if(*description)
    write_dummy_attributes_in_hdf5(description, dataset);
  H5Dclose(dataset);
  H5Sclose(dataspace);
}

void subfind_write_hdf5_MyDouble_particle(MyDouble * data, char *tag, int elems, char *description,
					  int blocknr, hid_t grp)
{
  hid_t dataspace, dataset;

  hsize_t adim[1];

  adim[0] = elems * Nids;
  dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(dataspace, 1, adim, NULL);

#ifndef DOUBLEPRECISION		/* default is single-precision */
  dataset = H5Dcreate(grp, tag, H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
#else
#if (DOUBLEPRECISION == 1)	/* everything double-precision */
  dataset = H5Dcreate(grp, tag, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
#else
#if (DOUBLEPRECISION == 2)	/* mixed precision */
  dataset = H5Dcreate(grp, tag, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
#endif
#endif
#endif
  if(blocknr >= 0)
    write_attributes_in_hdf5((enum iofields) blocknr, dataset);
  if(*description)
    write_dummy_attributes_in_hdf5(description, dataset);
  H5Dclose(dataset);
  H5Sclose(dataspace);
}

void subfind_write_hdf5_int_particle(int *data, char *tag, int elems, char *description, int blocknr,
				     hid_t grp)
{
  hid_t dataspace, dataset;

  hsize_t adim[1];

  adim[0] = elems * Nids;
  dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(dataspace, 1, adim, NULL);
  dataset = H5Dcreate(grp, tag, H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if(blocknr >= 0)
    write_attributes_in_hdf5((enum iofields) blocknr, dataset);
  if(*description)
    write_dummy_attributes_in_hdf5(description, dataset);
  H5Dclose(dataset);
  H5Sclose(dataspace);
}

void initialize_SubGroup(int size)
{
  int i, j;

  for(i = 0; i < size; i++)
    {
      SubGroup[i].Len = 0;
      for(j = 0; j < 6; j++)
	SubGroup[i].LenType[j] = 0;

      SubGroup[i].GrNr = 0;
      SubGroup[i].SubNr = 0;
      SubGroup[i].SubParent = 0;

      SubGroup[i].Offset = 0;
      for(j = 0; j < 6; j++)
	SubGroup[i].OffsetType[j] = 0;

      SubGroup[i].SubMostBoundID = 0;

      SubGroup[i].Mass = 0.0;
      for(j = 0; j < 6; j++)
	SubGroup[i].MassType[j] = 0.0;

      SubGroup[i].SubVelDisp = 0.0;
      SubGroup[i].SubStellarVelDisp = 0.0;
      SubGroup[i].SubStellarVelDispHalfProj = 0.0;
      SubGroup[i].SubVmax = 0.0;
      SubGroup[i].SubVmaxRad = 0.0;

      for(j = 0; j < 6; j++)
	{
	  SubGroup[i].SubHalfMass[j] = 0.0;
	  SubGroup[i].SubHalfMassProj[j] = 0.0;
	}

      for(j = 0; j < 3; j++)
	{
	  SubGroup[i].Pos[j] = 0.0;
	  SubGroup[i].CM[j] = 0.0;
	  SubGroup[i].Vel[j] = 0.0;
	  SubGroup[i].Spin[j] = 0.0;
	}
    }
}
#endif

int subfind_compare_ID_list(const void *a, const void *b)
{
  if(((struct id_list *) a)->GrNr < ((struct id_list *) b)->GrNr)
    return -1;

  if(((struct id_list *) a)->GrNr > ((struct id_list *) b)->GrNr)
    return +1;

  if(((struct id_list *) a)->SubNr < ((struct id_list *) b)->SubNr)
    return -1;

  if(((struct id_list *) a)->SubNr > ((struct id_list *) b)->SubNr)
    return +1;

  // Add another ordering over type
  if(((struct id_list *) a)->Type < ((struct id_list *) b)->Type)
    return -1;

  if(((struct id_list *) a)->Type > ((struct id_list *) b)->Type)
    return +1;

  if(((struct id_list *) a)->BindingEgy < ((struct id_list *) b)->BindingEgy)
    return -1;

  if(((struct id_list *) a)->BindingEgy > ((struct id_list *) b)->BindingEgy)
    return +1;

  return 0;
}

int subfind_compare_SubGroup_GrNr_SubNr(const void *a, const void *b)
{
  if(((struct subgroup_properties *) a)->GrNr < ((struct subgroup_properties *) b)->GrNr)
    return -1;

  if(((struct subgroup_properties *) a)->GrNr > ((struct subgroup_properties *) b)->GrNr)
    return +1;

  if(((struct subgroup_properties *) a)->SubNr < ((struct subgroup_properties *) b)->SubNr)
    return -1;

  if(((struct subgroup_properties *) a)->SubNr > ((struct subgroup_properties *) b)->SubNr)
    return +1;

  return 0;
}


int subfind_compare_P_GrNr_DM_Density(const void *a, const void *b)
{
  if(((struct particle_data *) a)->GrNr < (((struct particle_data *) b)->GrNr))
    return -1;

  if(((struct particle_data *) a)->GrNr > (((struct particle_data *) b)->GrNr))
    return +1;

  if(((struct particle_data *) a)->u.DM_Density > (((struct particle_data *) b)->u.DM_Density))
    return -1;

  if(((struct particle_data *) a)->u.DM_Density < (((struct particle_data *) b)->u.DM_Density))
    return +1;

  return 0;
}


int subfind_compare_Ind_GrNr_DM_Density(const void *a, const void *b)
{
  if(((struct my_serial_index *) a)->GrNr < (((struct my_serial_index *) b)->GrNr))
    return -1;

  if(((struct my_serial_index *) a)->GrNr > (((struct my_serial_index *) b)->GrNr))
    return +1;

  if(((struct my_serial_index *) a)->Density > (((struct my_serial_index *) b)->Density))
    return -1;

  if(((struct my_serial_index *) a)->Density < (((struct my_serial_index *) b)->Density))
    return +1;

  return 0;
}


#ifdef BG_SFR
void test_particle_sequence(char *message)
{
  int i;

  if(ThisTask == 0)
    printf("%s\n", message);

  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 4)
      {
	if((i != (int) StarP[P[i].StarID].PID) || ((int) P[i].StarID >= N_star))
	  {
	    printf
	      ("************\n  [rearrange_particle_sequence] Warning in rearrange on Task %i at particle %i\n************\n"
	       "P.StarID = %u  StarP.ID = %u N_star = %i\n", ThisTask, i, P[i].StarID, StarP[P[i].StarID].PID,
	       N_star);
	    fflush(stdout);
	    endrun(1212);
	  }
      }
}
#endif

#endif
