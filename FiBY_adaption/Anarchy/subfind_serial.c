#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#ifdef SUBFIND
#include "subfind.h"
#include "fof.h"

/* this file processes the local groups in serial mode */

#define NEAREST(x) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))


static int *Head, *Next, *Tail, *Len;

static struct cand_dat
{
  int head;
  int len;
  int lentype[6];
  int nsub;
  int rank, subnr, parent;
  int bound_length;
}
 *candidates;


int subfind_process_group_serial(int gr, int Offs, int *iter, int *skip)
{
  int i, j, k, p, len, subnr, totlen, ss, ngbs, ndiff, N, head = 0, head_attach, count_cand;

  int lentype[6], totlentype;

  int skipunbind, stmp;

  int extraDesLinkNgb;

  int listofdifferent[2], count, prev;

  int ngb_index, part_index, nsubs, rank;

  double SubMass, SubMassType[6], SubPos[3], SubVel[3], SubCM[3], SubVelDisp, SubStellarVelDisp,
    SubStellarVelDispHalfProj, SubVmax, SubVmaxRad, SubSpin[3], SubHalfMass[6], SubHalfMassProj[6];
  MyIDType SubMostBoundID;

  static struct unbind_data *ud;

  while(P[Ind[Offs].Index].GrNr != Group[gr].GrNr)
    {
      Offs++;
      if(Offs >= NumPart)
	{
	  printf("don't find a particle for groupnr=%d\n", Group[gr].GrNr);
	  endrun(312);
	}
    }

  N = Group[gr].Len;
  GrNr = Group[gr].GrNr;

  for(i = 0; i < N; i++)
    {
      if(P[Ind[Offs + i].Index].GrNr != GrNr)
	{
	  printf
	    ("task=%d, gr=%d: don't have the number of particles for GrNr=%d group-len=%d found=%d before=%d\n",
	     ThisTask, gr, Group[gr].GrNr, N, P[Ind[Offs + i].Index].GrNr, P[Ind[Offs - 1].Index].GrNr);
	  endrun(312);
	}

      if(P[Ind[Offs + i].Index].u.DM_Density < 0)
	{
	  printf("task=%d i=%d/%d %g\n", ThisTask, Ind[Offs + i].Index, NumPart,
		 P[Ind[Offs + i].Index].u.DM_Density);
	  endrun(7120);
	}
    }


  candidates = mymalloc(N * sizeof(struct cand_dat));

  Head = mymalloc(N * sizeof(int));
  Next = mymalloc(N * sizeof(int));
  Tail = mymalloc(N * sizeof(int));
  Len = mymalloc(N * sizeof(int));
  ud = mymalloc(N * sizeof(struct unbind_data));

  Head -= Offs;
  Next -= Offs;
  Tail -= Offs;
  Len -= Offs;

  for(i = 0; i < N; i++)
    ud[i].index = Ind[Offs + i].Index;

  subfind_loctree_findExtent(N, ud);

  subfind_loctree_treebuild(N, ud);	/* build tree for all particles of this group */

  for(i = Offs; i < Offs + N; i++)
    Head[i] = Next[i] = Tail[i] = -1;


  /* note: particles are already ordered in the order of decreasing density */
  for(i = 0, count_cand = 0; i < N; i++)
    {
      part_index = Offs + i;

      subfind_locngb_treefind(P[Ind[part_index].Index].Pos, All.DesLinkNgb, P[Ind[part_index].Index].DM_Hsml);

      // If N < DesLinkNgb then the search for more particles will enter an infinite loop, so use the minimum FOF size instead (which is obviously small enough to ensure that enough particles will be found) - this isn't included in parameter file but in allvars.h as GROUP_MIN_LEN - Alan
      /* note: returned neighbours are already sorted by distance */
      for(k = 0, ndiff = 0, ngbs = 0; k < All.DesLinkNgb && ngbs < 2; k++)
	{
	  ngb_index = Rank[R2list[k].index];

	  if(ngb_index != part_index)	/* to exclude the particle itself */
	    {
	      /* we only look at neighbours that are denser */
	      if(P[Ind[ngb_index].Index].u.DM_Density > P[Ind[part_index].Index].u.DM_Density)
		{
		  ngbs++;

		  if(Head[ngb_index] >= 0)	/* neighbor is attached to a group */
		    {
		      if(ndiff == 1)
			if(listofdifferent[0] == Head[ngb_index])
			  continue;

		      /* a new group has been found */
		      listofdifferent[ndiff++] = Head[ngb_index];
		    }
		  else
		    {
		      printf("this may not occur.\n");
		      printf
			("ThisTask=%d gr=%d k=%d i=%d part_index=%d ngb_index = %d  head[ngb_index]=%d P[part_index].DM_Density=%g %g GrNrs= %d %d \n",
			 ThisTask, gr, k, i, Ind[part_index].Index, ngb_index, Head[ngb_index],
			 P[Ind[part_index].Index].u.DM_Density, P[Ind[ngb_index].Index].u.DM_Density,
			 P[Ind[part_index].Index].GrNr, P[Ind[ngb_index].Index].GrNr);
		      endrun(2);
		    }
		}
	    }
	}

      switch (ndiff)		/* treat the different possible cases */
	{
	case 0:		/* this appears to be a lonely maximum -> new group */
	  head = part_index;
	  Head[part_index] = Tail[part_index] = part_index;
	  Len[part_index] = 1;
	  Next[part_index] = -1;
	  break;

	case 1:		/* the particle is attached to exactly one group */
	  head = listofdifferent[0];
	  Head[part_index] = head;
	  Next[Tail[head]] = part_index;
	  Tail[head] = part_index;
	  Len[head]++;
	  Next[part_index] = -1;
	  break;

	case 2:		/* the particle merges two groups together */
	  head = listofdifferent[0];
	  head_attach = listofdifferent[1];

	  if(Len[head_attach] > Len[head])	/* other group is longer, swap them */
	    {
	      head = listofdifferent[1];
	      head_attach = listofdifferent[0];
	    }

	  /* only in case the attached group is long enough we bother to register is 
	     as a subhalo candidate */

	  // Have the routine count different particle types here - Alan
	  if(Len[head_attach] >= All.DesLinkNgb)
	    {
	      candidates[count_cand].len = Len[head_attach];
	      candidates[count_cand].head = Head[head_attach];
	      count_cand++;
	    }

	  /* now join the two groups */
	  Next[Tail[head]] = head_attach;
	  Tail[head] = Tail[head_attach];
	  Len[head] += Len[head_attach];

	  ss = head_attach;
	  do
	    {
	      Head[ss] = head;
	    }
	  while((ss = Next[ss]) >= 0);

	  /* finally, attach the particle */
	  Head[part_index] = head;
	  Next[Tail[head]] = part_index;
	  Tail[head] = part_index;
	  Len[head]++;
	  Next[part_index] = -1;	//If nothing else is attached then this group has a definite end (marked by -1) - Alan
	  break;

	default:
	  printf("can't be! (a)\n");
	  endrun(1);
	  break;
	}
    }


  /* add the full thing as a subhalo candidate */
  for(i = 0, prev = -1; i < N; i++)
    {
      if(Head[Offs + i] == Offs + i)
	if(Next[Tail[Offs + i]] == -1)
	  {
	    if(prev < 0)
	      head = Offs + i;
	    if(prev >= 0)
	      Next[prev] = Offs + i;

	    prev = Tail[Offs + i];
	  }
    }

  candidates[count_cand].len = N;
  candidates[count_cand].head = head;
  count_cand++;

  /* go through them once and assign the rank */
  for(i = 0, p = head, rank = 0; i < N; i++)
    {
      Len[p] = rank++;
      p = Next[p];
    }

  /* for each candidate, we now pull out the rank of its head */
  for(k = 0; k < count_cand; k++)
    candidates[k].rank = Len[candidates[k].head];

  for(i = Offs; i < Offs + N; i++)
    Tail[i] = -1;

  stmp = *skip;
  for(k = 0, nsubs = 0; k < count_cand; k++)
    {
      for(j = 0; j < 6; j++)
	lentype[j] = 0;		// Internal counter for particle number per type, will then be subtracted from len (ignored thus) depending on make option
      for(i = 0, p = candidates[k].head, len = 0; i < candidates[k].len; i++, p = Next[p])
	if(Tail[p] < 0)
	  ud[len++].index = Ind[p].Index;

      if(len >= All.DesLinkNgb)
	{
	  len = subfind_unbind(ud, len, &lentype[0], iter, &skipunbind);
	  if(skipunbind == 1)
	    SkipGrUnbind[stmp++] = gr;	// Group number for which unbind has been skipped, i.e. the halo retains its fuzz.    
	}

      extraDesLinkNgb = 0;
#ifdef DONOT_COUNT_GAS
      extraDesLinkNgb += lentype[0];
#endif
#ifdef DONOT_COUNT_DM
      extraDesLinkNgb += lentype[1];
#endif
#ifdef DONOT_COUNT_STARS
      extraDesLinkNgb += lentype[4];
#endif

      if(len >= (All.DesLinkNgb + extraDesLinkNgb))
	{
	  /* ok, we found a substructure */
	  for(i = 0; i < 6; i++)
	    candidates[k].lentype[i] = 0;	//lentype is set to zero here - Alan

	  for(i = 0; i < len; i++)
	    {
	      Tail[Rank[ud[i].index]] = nsubs;	/* we use this to flag the substructures */
	      candidates[k].lentype[P[ud[i].index].Type]++;	//lentype is created here - Alan
	    }
	  // Just to check that lentype and candidates.lentype agree
	  for(i = 0; i < 6; i++)
	    if(candidates[k].lentype[i] != lentype[i])
	      {
		printf("candidates[%d].lentype[%d] = %d, lentype[%d] = %d; they should be equal! \n", k, i,
		       candidates[k].lentype[i], i, lentype[i]);
		for(j = 0; j < 6; j++)
		  printf("lentype[%d] = %d \n", j, lentype[j]);
		for(j = 0; j < len; j++)
		  printf("ud[%d].index = %d, P[%d].Type is %d \n", j, ud[j].index, ud[j].index,
			 P[ud[j].index].Type);
		endrun(811);
	      }
	  // End of check
	  candidates[k].nsub = nsubs;
	  candidates[k].bound_length = len;
	  nsubs++;
	}
      else
	{
	  candidates[k].nsub = -1;
	  candidates[k].bound_length = 0;
	}
    }

#ifdef VERBOSE
  printf("\nGroupLen=%d  (gr=%d)\n", N, gr);
  printf("Number of substructures: %d\n", nsubs);
#endif

  Group[gr].Nsubs = nsubs;

  qsort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_boundlength);

  /* now we determine the parent subhalo for each candidate */

  for(k = 0; k < count_cand; k++)
    {
      candidates[k].subnr = k;
      candidates[k].parent = 0;
    }

  qsort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_rank);


  for(k = 0; k < count_cand; k++)
    {
      for(j = k + 1; j < count_cand; j++)
	{
	  if(candidates[j].rank > candidates[k].rank + candidates[k].len)
	    break;

	  if(candidates[k].rank + candidates[k].len >= candidates[j].rank + candidates[j].len)
	    {
	      if(candidates[k].bound_length >= All.DesLinkNgb)
		candidates[j].parent = candidates[k].subnr;
	    }
	  else
	    {
	      printf("k=%d|%d has rank=%d and len=%d.  j=%d has rank=%d and len=%d bound=%d\n",
		     k, count_cand, (int) candidates[k].rank, candidates[k].len,
		     (int) candidates[k].bound_length, candidates[j].rank,
		     (int) candidates[j].len, candidates[j].bound_length);
	      endrun(121235513);
	    }

	}
    }
  *skip = stmp;

  qsort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_subnr);

  /* now determine the properties */

  for(k = 0, subnr = 0, totlen = 0, totlentype = 0; k < nsubs; k++)
    {
      len = candidates[k].bound_length;

      for(j = 0; j < 6; j++)
	{
	  lentype[j] = candidates[k].lentype[j];

	  totlentype += lentype[j];
	}
      //Ouput the total the particle numbers per type per subhalo into lentype[] - Alan

#ifdef VERBOSE
      printf("subnr=%d  SubLen=%d\n", subnr, len);
#endif
      totlen += len;

      if(totlentype != totlen)
	printf("Error! totlen = %d while totlentype = %d \n", totlen, totlentype);
      //Quick check to ensure the individual types are not going awry - Alan

      for(i = 0, p = candidates[k].head, count = 0; i < candidates[k].len; i++)
	{
	  if(Tail[p] == candidates[k].nsub)
	    ud[count++].index = Ind[p].Index;

	  p = Next[p];
	}

      if(count != len)
	endrun(12);

      subfind_determine_sub_halo_properties(ud, len, &SubMass, &SubMassType[0],
					    &SubPos[0], &SubVel[0], &SubCM[0],
					    &SubVelDisp, &SubStellarVelDisp,
					    &SubStellarVelDispHalfProj,
					    &SubVmax, &SubVmaxRad, &SubSpin[0],
					    &SubMostBoundID, &SubHalfMass[0], &SubHalfMassProj[0]);

      if(Nsubgroups >= MaxNsubgroups)
	endrun(899);

      if(subnr == 0)
	{
	  for(j = 0; j < 3; j++)
	    Group[gr].Pos[j] = SubPos[j];
	}

      SubGroup[Nsubgroups].Len = len;
      for(j = 0; j < 6; j++)
	SubGroup[Nsubgroups].LenType[j] = lentype[j];	//Final array for LenType - Alan


      //      if(Nsubgroups == 0)
      //        SubGroup[Nsubgroups].Offset = 0;
      if(subnr == 0)
	SubGroup[Nsubgroups].Offset = Group[gr].Offset;	// - Previously was set like this
      else
	SubGroup[Nsubgroups].Offset = SubGroup[Nsubgroups - 1].Offset + SubGroup[Nsubgroups - 1].Len;

      // Set the subgroup offset to zero - Alan
      for(j = 0; j < 6; j++)
	SubGroup[Nsubgroups].OffsetType[j] = 0;

      //      if(Nsubgroups == 0)
      if(subnr == 0)
	for(j = 0; j < 6; j++)
	  {
	    if(j == 0)
	      SubGroup[Nsubgroups].OffsetType[j] = Group[gr].Offset;	// This means that one must read in the ParticleID list from the FOF group files 
	    //      if(j == 0)
	    //        SubGroup[Nsubgroups].OffsetType[j] = 0; // This means that one must read in the ParticleID list from the FOF group files 
	    else		// or I could output only the remaining IDs like in L-SubFind and modify this line 
	      SubGroup[Nsubgroups].OffsetType[j] =
		SubGroup[Nsubgroups].OffsetType[j - 1] + SubGroup[Nsubgroups].LenType[j - 1];
	  }
      else
	for(j = 0; j < 6; j++)
	  {
	    if(j == 0)
	      SubGroup[Nsubgroups].OffsetType[j] =
		SubGroup[Nsubgroups - 1].OffsetType[5] + SubGroup[Nsubgroups - 1].LenType[5];
	    else
	      SubGroup[Nsubgroups].OffsetType[j] =
		SubGroup[Nsubgroups].OffsetType[j - 1] + SubGroup[Nsubgroups].LenType[j - 1];
	  }
      //Create the Offset per type for each subhalo in the particle ID list - Alan
      //      for(j = 0; j < 6; j++)
      //        printf("SubGroup[%d].Offset = %d while SubGroup[%d].OffsetType[%d] = %d \n",Nsubgroups,SubGroup[Nsubgroups].Offset,Nsubgroups,j,SubGroup[Nsubgroups].OffsetType[j]);

      if(SubGroup[Nsubgroups].Offset != SubGroup[Nsubgroups].OffsetType[0])
	printf("Problem! SubGroup[%d].Offset = %d while (SubGroup[%d].OffsetType[0] = %d) \n", Nsubgroups,
	       SubGroup[Nsubgroups].Offset, (Nsubgroups), SubGroup[Nsubgroups].OffsetType[0]);
      //Quick check to ensure the individual offset types make sense - Alan

      SubGroup[Nsubgroups].GrNr = GrNr - 1;
      SubGroup[Nsubgroups].SubNr = subnr;
      SubGroup[Nsubgroups].SubParent = candidates[k].parent;
      SubGroup[Nsubgroups].Mass = SubMass;

      for(j = 0; j < 6; j++)
	{
	  SubGroup[Nsubgroups].MassType[j] = SubMassType[j];
	  SubGroup[Nsubgroups].SubHalfMass[j] = SubHalfMass[j];
	  SubGroup[Nsubgroups].SubHalfMassProj[j] = SubHalfMassProj[j];
	}

      SubGroup[Nsubgroups].SubMostBoundID = SubMostBoundID;
      SubGroup[Nsubgroups].SubVelDisp = SubVelDisp;
      SubGroup[Nsubgroups].SubStellarVelDisp = SubStellarVelDisp;
      SubGroup[Nsubgroups].SubStellarVelDispHalfProj = SubStellarVelDispHalfProj;
      SubGroup[Nsubgroups].SubVmax = SubVmax;
      SubGroup[Nsubgroups].SubVmaxRad = SubVmaxRad;

      for(j = 0; j < 3; j++)
	{
	  SubGroup[Nsubgroups].Pos[j] = SubPos[j];
	  SubGroup[Nsubgroups].CM[j] = SubCM[j];
	  SubGroup[Nsubgroups].Vel[j] = SubVel[j];
	  SubGroup[Nsubgroups].Spin[j] = SubSpin[j];
	}

      Nsubgroups++;

      /* Let's now assign the subgroup number */

      for(i = 0; i < len; i++)
	P[ud[i].index].SubNr = subnr;

      subnr++;
    }


#ifdef VERBOSE
  printf("Fuzz=%d\n", N - totlen);
#endif

  myfree(ud);
  myfree(Len + Offs);
  myfree(Tail + Offs);
  myfree(Next + Offs);
  myfree(Head + Offs);

  myfree(candidates);

  return Offs;
}




int subfind_unbind(struct unbind_data *ud, int len, int *lentype, int *iter, int *skipunbind)
{
  double *bnd_energy, energy_limit, weakly_bound_limit = 0;

  int i, j, p, minindex, unbound, phaseflag;

  double ddxx, s[3], dx[3], v[3], dv[3], pos[3];

  double vel_to_phys, H_of_a, pot, minpot = 0;

  double boxsize, boxhalf;

  double TotMass;

  double u, a3inv;

  int itmp = 0;

  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;

  vel_to_phys = 1.0 / All.Time;

  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  H_of_a =
    All.Omega0 / (All.Time * All.Time * All.Time) + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time)
#ifdef DARKENERGY
    + DarkEnergy_a(All.Time);
#else
    + All.OmegaLambda;
#endif
  H_of_a = All.Hubble * sqrt(H_of_a);

  bnd_energy = (double *) mymalloc(len * sizeof(double));

  phaseflag = 0;		/* this means we will recompute the potential for all particles */

  *iter = 0;
  for(i = 0; i < 6; i++)
    lentype[i] = 0;

  do
    {
      itmp++;

      subfind_loctree_treebuild(len, ud);

      /* let's compute the potential  */

      if(phaseflag == 0)	/* redo it for all the particles */
	{
	  for(i = 0, minindex = -1, minpot = 1.0e30; i < len; i++)
	    {
	      p = ud[i].index;

	      pot = subfind_loctree_treeevaluate_potential(p);
	      /* note: add self-energy */
	      P[p].u.DM_Potential = pot + P[p].Mass / All.SofteningTable[P[p].Type];
	      P[p].u.DM_Potential *= All.G / All.Time;
	      //This includes baryons already - provided there exists a ForceSoftening lengthtable for the particle types - Alan

	      /* if(P[p].Type > 0) *///Don't want gas particle being minimum potential point
	      if(P[p].u.DM_Potential < minpot || minindex == -1)
		{
		  minpot = P[p].u.DM_Potential;
		  minindex = p;
		}
	    }

	  for(j = 0; j < 3; j++)
	    pos[j] = P[minindex].Pos[j];	/* position of minimum potential */
	}
      else
	{
	  /* we only repeat for those close to the unbinding threshold */
	  for(i = 0; i < len; i++)
	    {
	      p = ud[i].index;

	      if(P[p].v.DM_BindingEnergy >= weakly_bound_limit)
		{
		  //              pot = subfind_loctree_treeevaluate_potential(p);
		  pot = subfind_loctree_treeevaluate_accurate_potential(p);	// Use a more accurate version of the potential calculation
		  /* note: add self-energy */
		  P[p].u.DM_Potential = pot + P[p].Mass / All.SofteningTable[P[p].Type];
		  P[p].u.DM_Potential *= All.G / All.Time;
		}
	    }
	}

      /* let's get bulk velocity and the center-of-mass */

      v[0] = v[1] = v[2] = 0;
      s[0] = s[1] = s[2] = 0;

      for(i = 0, TotMass = 0; i < len; i++)
	{
	  p = ud[i].index;

	  for(j = 0; j < 3; j++)
	    {
#ifdef PERIODIC
	      ddxx = NEAREST(P[p].Pos[j] - pos[j]);
#else
	      ddxx = P[p].Pos[j] - pos[j];
#endif
	      s[j] += P[p].Mass * ddxx;
	      v[j] += P[p].Mass * P[p].Vel[j];
	    }
	  TotMass += P[p].Mass;
	}

      for(j = 0; j < 3; j++)
	{
	  v[j] /= TotMass;
	  s[j] /= TotMass;	/* center-of-mass */

	  s[j] += pos[j];

#ifdef PERIODIC
	  while(s[j] < 0)
	    s[j] += boxsize;
	  while(s[j] >= boxsize)
	    s[j] -= boxsize;
#endif
	}

      for(i = 0; i < 6; i++)
	lentype[i] = 0;

      for(i = 0; i < len; i++)
	{
	  p = ud[i].index;
	  lentype[P[p].Type]++;

	  for(j = 0; j < 3; j++)
	    {
	      dv[j] = vel_to_phys * (P[p].Vel[j] - v[j]);
#ifdef PERIODIC
	      dx[j] = All.Time * NEAREST(P[p].Pos[j] - s[j]);
#else
	      dx[j] = All.Time * (P[p].Pos[j] - s[j]);
#endif
	      dv[j] += H_of_a * dx[j];
	    }

	  P[p].v.DM_BindingEnergy =
	    P[p].u.DM_Potential + 0.5 * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
	  //specific potential plus kinetic energy terms - Alan

	  //Added this internal energy term - Alan
	  if(P[p].Type == 0)
	    {
	      u = SphP[p].Entropy / GAMMA_MINUS1 * pow(SphP[p].d.Density * a3inv, GAMMA_MINUS1);
	      P[p].v.DM_BindingEnergy += u;
	    }

	  bnd_energy[i] = P[p].v.DM_BindingEnergy;
	}

      qsort(bnd_energy, len, sizeof(double), subfind_compare_binding_energy);	/* largest comes first! */

      *skipunbind = 0;
/* #ifdef IGNOREFUZZ */
/*       // Check that the object isn't too large; if so then ignore unbinding of Fuzz */
/*       if(len <= sizelimit) */
/* 	{ */
/* #endif */
      energy_limit = bnd_energy[(int) (0.25 * len)];

      for(i = 0, unbound = 0; i < len - 1; i++)
	{
	  if(bnd_energy[i] > 0)
	    unbound++;
	  else
	    unbound--;

	  if(unbound <= 0)
	    break;
	}
      weakly_bound_limit = bnd_energy[i];

      /* now omit unbound particles, but at most 1/4 of the original size */

      for(i = 0, unbound = 0; i < len; i++)
	{
	  p = ud[i].index;
	  if(P[p].v.DM_BindingEnergy > 0 && P[p].v.DM_BindingEnergy > energy_limit)
	    {
	      unbound++;
	      ud[i] = ud[len - 1];
	      i--;
	      len--;
	      lentype[P[p].Type]--;
	    }
	}

      if(len < All.DesLinkNgb)
	break;

      if(phaseflag == 0)
	{
	  if(unbound > 0)
	    phaseflag = 1;
	}
      else
	{
	  if(unbound == 0)
	    {
	      phaseflag = 0;	/* this will make us repeat everything once more for all particles */
	      unbound = 1;
	    }
	}
/* #ifdef IGNOREFUZZ */
/*       	} // If group has len < sizelimit then it will run through these calculations, and unbind fuzz, else it calculates potential once and retains unbound particles */
/*       else */
/*       	{ */
/*       	  unbound = 0; */
/*       	  *skipunbind = 1; */
/*       	  printf("\n"); */
/*       	} */
/* #endif */
    }
  while(unbound > 0);

  myfree(bnd_energy);

  *iter = itmp;

  return (len);
}



int subfind_compare_grp_particles(const void *a, const void *b)
{
  if(((struct particle_data *) a)->GrNr < ((struct particle_data *) b)->GrNr)
    return -1;

  if(((struct particle_data *) a)->GrNr > ((struct particle_data *) b)->GrNr)
    return +1;

  if(((struct particle_data *) a)->SubNr < ((struct particle_data *) b)->SubNr)
    return -1;

  if(((struct particle_data *) a)->SubNr > ((struct particle_data *) b)->SubNr)
    return +1;

  if(((struct particle_data *) a)->v.DM_BindingEnergy < ((struct particle_data *) b)->v.DM_BindingEnergy)
    return -1;

  if(((struct particle_data *) a)->v.DM_BindingEnergy > ((struct particle_data *) b)->v.DM_BindingEnergy)
    return +1;

  return 0;
}



void subfind_determine_sub_halo_properties(struct unbind_data *d, int num, double *totmass,
					   double *totmasstype, double *pos, double *vel, double *cm,
					   double *veldisp, double *stellarveldisp,
					   double *halfproj_stellarveldisp, double *vmax, double *vmaxrad,
					   double *spin, MyIDType * mostboundid, double *halfmassrad,
					   double *proj_halfmassrad)
{
  int i, j, p;

  double s[3], v[3], max, vel_to_phys, H_of_a, minpot;

  double lx, ly, lz, dv[3], dx[3], disp[6], masstype[6];

  double xdisp[6], ydisp[6], zdisp[6], xhalfmassrad[6], yhalfmassrad[6], zhalfmassrad[6];

  double xy, xz, yz, sum[6];

  double boxhalf, boxsize, ddxx;

  sort_r2list *rr_list = 0;

  int minindex;

  double mass, maxrad;


  boxsize = All.BoxSize;
  boxhalf = 0.5 * boxsize;

  vel_to_phys = 1.0 / All.Time;

  H_of_a =
    All.Hubble * sqrt(All.Omega0 / (All.Time * All.Time * All.Time) +
		      (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) + All.OmegaLambda);


  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
      p = d[i].index;
      /* if(P[p].Type > 0) *///Don't want gas particle being minimum potential point
      if(P[p].u.DM_Potential < minpot || minindex == -1)
	{
	  minpot = P[p].u.DM_Potential;
	  minindex = p;
	}
    }

  for(j = 0; j < 3; j++)
    pos[j] = P[minindex].Pos[j];	//Comoving coordinates - Alan


  /* pos[] now holds the position of minimum potential */
  /* we take it that as the center */


  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
      p = d[i].index;
      /* if(P[p].Type > 0) *///Don't want gas particle being minimum potential point
      if(P[p].v.DM_BindingEnergy < minpot || minindex == -1)
	{
	  minpot = P[p].v.DM_BindingEnergy;
	  minindex = p;
	}
    }

  *mostboundid = P[minindex].ID;
  //  printf("mostboundid is %d \n",*mostboundid);

  /* let's get bulk velocity and the center-of-mass */

  for(j = 0; j < 3; j++)
    s[j] = v[j] = 0;

  // Have total mass per type
  for(j = 0; j < 6; j++)
    masstype[j] = 0.0;

  for(i = 0, mass = 0; i < num; i++)
    {
      p = d[i].index;
      for(j = 0; j < 3; j++)
	{
#ifdef PERIODIC
	  ddxx = NEAREST(P[p].Pos[j] - pos[j]);
#else
	  ddxx = P[p].Pos[j] - pos[j];
#endif
	  s[j] += P[p].Mass * ddxx;
	  v[j] += P[p].Mass * P[p].Vel[j];
	}
      mass += P[p].Mass;
      masstype[P[p].Type] += P[p].Mass;
    }

  *totmass = mass;

  // analagous to totmass variable, but per particle type
  for(j = 0; j < 6; j++)
    totmasstype[j] = masstype[j];

  for(j = 0; j < 3; j++)
    {
      s[j] /= mass;		/* center of mass */
      v[j] /= mass;

      vel[j] = vel_to_phys * v[j];
    }

  for(j = 0; j < 3; j++)
    {
      s[j] += pos[j];

#ifdef PERIODIC
      while(s[j] < 0)
	s[j] += boxsize;
      while(s[j] >= boxsize)
	s[j] -= boxsize;
#endif

      cm[j] = s[j];		//Comoving coordinates
    }

  // Have dispersion for each particle type
  for(j = 0; j < 6; j++)
    disp[j] = 0.0;

  lx = ly = lz = 0;

  rr_list = mymalloc(sizeof(sort_r2list) * num);

  for(i = 0; i < num; i++)
    {
      p = d[i].index;
      rr_list[i].index = p;	// Original particle index
      rr_list[i].r = 0.0;
      // 3 cartesian planes for projected radius
      rr_list[i].xproj = 0.0;
      rr_list[i].yproj = 0.0;
      rr_list[i].zproj = 0.0;
      // x,y,z above refer to the LoS axis
      rr_list[i].mass = P[p].Mass;
      rr_list[i].type = P[p].Type;
      xy = xz = yz = 0.;

      for(j = 0; j < 3; j++)
	{
#ifdef PERIODIC
	  ddxx = NEAREST(P[p].Pos[j] - s[j]);
#else
	  ddxx = P[p].Pos[j] - s[j];
#endif
	  dx[j] = All.Time * ddxx;
	  dv[j] = vel_to_phys * (P[p].Vel[j] - v[j]);
	  dv[j] += H_of_a * dx[j];	/* BUG HERE */

	  disp[P[p].Type] += P[p].Mass * dv[j] * dv[j];
	  /* for rotation curve computation, take minimum of potential as center */
#ifdef PERIODIC
	  ddxx = NEAREST(P[p].Pos[j] - pos[j]);
#else
	  ddxx = P[p].Pos[j] - pos[j];
#endif
	  ddxx = All.Time * ddxx;

	  // Construct the xy, xz and yz plane distances for the projected distribution
	  if(j == 0)
	    {
	      xy += ddxx * ddxx;
	      xz += ddxx * ddxx;
	    }
	  if(j == 1)
	    {
	      xy += ddxx * ddxx;
	      yz += ddxx * ddxx;
	    }
	  if(j == 2)
	    {
	      xz += ddxx * ddxx;
	      yz += ddxx * ddxx;
	    }
	  rr_list[i].r += ddxx * ddxx;
	  rr_list[i].xproj = yz;	// Proper distance of projected radius with LoS along x-axis   
	  rr_list[i].yproj = xz;	// Proper distance of projected radius with LoS along y-axis
	  rr_list[i].zproj = xy;	// Proper distance of projected radius with LoS along z-axis
	}
      //dx is a*x, proper, dv is v_int/a and, since v_int is a*a, this is peculiar velocity so lx is proper as well.
      lx += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
      ly += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
      lz += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

      rr_list[i].r = sqrt(rr_list[i].r);	//Proper distance
    }

  if(masstype[1] > 0.0)
    *veldisp = sqrt(disp[1] / (3. * masstype[1]));	/* convert to 1d velocity dispersion -DM only */
  else
    *veldisp = 0.0;		/* there were no DM particles in this halo so set to zero */

  if(masstype[4] > 0.0)
    *stellarveldisp = sqrt(disp[4] / (3. * masstype[4]));	/* convert to 1d stellar velocity dispersion */
  else
    *stellarveldisp = 0.0;	/* there were no star particles, so set dispersion to zero */

  spin[0] = lx / mass;		// from lx there is only dx and dv terms, both are proper although the vel and pos both carry h-1
  spin[1] = ly / mass;
  spin[2] = lz / mass;

  // Summation of the masses to determine projected half-mass radii
  // yz plane; x LoS 
  qsort(rr_list, num, sizeof(sort_r2list), subfind_compare_xproj_rotcurve);
  for(j = 0; j < 6; j++)
    xhalfmassrad[j] = sum[j] = xdisp[j] = 0.0;
  for(j = 0; j < 6; j++)
    {
      if(masstype[j] == 0.)
	continue;		// No reason to continue loop over particle types if this is zero
      for(i = 0; i < num; i++)
	{
	  if(j != rr_list[i].type)
	    continue;		// Only loop when particles are of the correct type
	  p = rr_list[i].index;
#ifdef PERIODIC
	  ddxx = NEAREST(P[p].Pos[0] - s[0]);
#else
	  ddxx = P[p].Pos[0] - s[0];
#endif
	  dx[0] = All.Time * ddxx;
	  dv[0] = vel_to_phys * (P[p].Vel[0] - v[0]);
	  dv[0] += H_of_a * dx[0];
	  xdisp[j] += P[p].Mass * dv[0] * dv[0];
	  sum[j] += rr_list[i].mass;	// Begin summating the masses in halo of type 'j' from centre out 
	  if(sum[j] > (masstype[j] / 2.))
	    {
	      if((sum[j] - (masstype[j] / 2.)) + ((sum[j] - rr_list[i].mass) - masstype[j] / 2.) > 0)	// The new particle takes the average mass further from mean than not adding it
		{
		  sum[j] -= rr_list[i].mass;	// Remove this particle from mass total; halfmassrad is still at previous particle value
		  break;
		}
	      else
		{
		  xhalfmassrad[j] = rr_list[i].xproj;
		  break;	// Leave loop when the half mass radius is reached
		}
	    }
	  xhalfmassrad[j] = rr_list[i].xproj;	// This will obviously be imprecise to one particle spacing, before real halfmass radius or after - depending on which is closer 
	}			// Loop over particles in halo
      xdisp[j] = sqrt(xdisp[j] / sum[j]);	// mass weighted average velocity dispersion within projected halfmass radius of type 'j' 
    }				// Loop over particle types

  // xz plane; y LoS 
  qsort(rr_list, num, sizeof(sort_r2list), subfind_compare_yproj_rotcurve);
  for(j = 0; j < 6; j++)
    yhalfmassrad[j] = sum[j] = ydisp[j] = 0.0;
  for(j = 0; j < 6; j++)
    {
      if(masstype[j] == 0.)
	continue;		// No reason to continue loop over particle types if this is zero
      for(i = 0; i < num; i++)
	{
	  if(j != rr_list[i].type)
	    continue;		// Only loop when particles are of the correct type
	  p = rr_list[i].index;
#ifdef PERIODIC
	  ddxx = NEAREST(P[p].Pos[1] - s[1]);
#else
	  ddxx = P[p].Pos[1] - s[1];
#endif
	  dx[1] = All.Time * ddxx;
	  dv[1] = vel_to_phys * (P[p].Vel[1] - v[1]);
	  dv[1] += H_of_a * dx[1];
	  ydisp[j] += P[p].Mass * dv[1] * dv[1];
	  sum[j] += rr_list[i].mass;	// Begin summating the masses in halo of type 'j' from centre out 
	  if(sum[j] > (masstype[j] / 2.))
	    {
	      if((sum[j] - (masstype[j] / 2.)) + ((sum[j] - rr_list[i].mass) - masstype[j] / 2.) > 0)	// The new particle takes the average mass further from mean than not adding it
		{
		  sum[j] -= rr_list[i].mass;	// Remove this particle from mass total; halfmassrad is still at previous particle value
		  break;
		}
	      else
		{
		  yhalfmassrad[j] = rr_list[i].yproj;
		  break;	// Leave loop when the half mass radius is reached
		}
	    }
	  yhalfmassrad[j] = rr_list[i].yproj;	// This will obviously be imprecise to one particle spacing, before real halfmass radius or after - depending on which is closer 
	}			// Loop over particles in halo
      ydisp[j] = sqrt(ydisp[j] / sum[j]);	// mass weighted average velocity dispersion within projected halfmass radius of type 'j' 
    }				// Loop over particle types

  // xy plane; z LoS 
  qsort(rr_list, num, sizeof(sort_r2list), subfind_compare_zproj_rotcurve);
  for(j = 0; j < 6; j++)
    zhalfmassrad[j] = sum[j] = zdisp[j] = 0.0;
  for(j = 0; j < 6; j++)
    {
      if(masstype[j] == 0.)
	continue;		// No reason to continue loop over particle types if this is zero
      for(i = 0; i < num; i++)
	{
	  if(j != rr_list[i].type)
	    continue;		// Only loop when particles are of the correct type
	  p = rr_list[i].index;
#ifdef PERIODIC
	  ddxx = NEAREST(P[p].Pos[2] - s[2]);
#else
	  ddxx = P[p].Pos[2] - s[2];
#endif
	  dx[2] = All.Time * ddxx;
	  dv[2] = vel_to_phys * (P[p].Vel[2] - v[2]);
	  dv[2] += H_of_a * dx[2];
	  zdisp[j] += P[p].Mass * dv[2] * dv[2];
	  sum[j] += rr_list[i].mass;	// Begin summating the masses in halo of type 'j' from centre out 
	  if(sum[j] > (masstype[j] / 2.))
	    {
	      if((sum[j] - (masstype[j] / 2.)) + ((sum[j] - rr_list[i].mass) - masstype[j] / 2.) > 0)	// The new particle takes the average mass further from mean than not adding it
		{
		  sum[j] -= rr_list[i].mass;	// Remove this particle from mass total; halfmassrad is still at previous particle value
		  break;
		}
	      else
		{
		  zhalfmassrad[j] = rr_list[i].zproj;
		  break;	// Leave loop when the half mass radius is reached
		}
	    }
	  zhalfmassrad[j] = rr_list[i].zproj;	// This will obviously be imprecise to one particle spacing, before real halfmass radius or after - depending on which is closer 
	}			// Loop over particles in halo
      zdisp[j] = sqrt(zdisp[j] / sum[j]);	// mass weighted average velocity dispersion within projected halfmass radius of type 'j' 
    }				// Loop over particle types

  *halfproj_stellarveldisp = 0;
  for(j = 0; j < 6; j++)
    {
      disp[j] = 0;
      if(masstype[j] == 0.)
	continue;
      proj_halfmassrad[j] = (xhalfmassrad[j] + yhalfmassrad[j] + zhalfmassrad[j]) / 3.;
      disp[j] = (xdisp[j] + ydisp[j] + zdisp[j]) / 3.;
      if(j == 4)
	*halfproj_stellarveldisp = disp[4];
    }				// Assign 4th coloumn of disp to stellarveldisp      

  // Summation of the masses to determine half-mass radii
  qsort(rr_list, num, sizeof(sort_r2list), subfind_compare_dist_rotcurve);
  for(j = 0; j < 6; j++)
    sum[j] = 0.0;
  for(j = 0; j < 6; j++)
    {
      halfmassrad[j] = 0.0;
      if(masstype[j] == 0.0)
	continue;		// No sense scanning halo when there are no particles of type 'j' present
      for(i = 0; i < num; i++)
	{
	  if(rr_list[i].type != j)
	    continue;		// Only take into account particles of type 'j'
	  sum[j] += rr_list[i].mass;	// Begin summating the masses in halo of type 'j' from centre out 
	  if(sum[j] > (masstype[j] / 2.))
	    {
	      if((sum[j] - (masstype[j] / 2.)) + ((sum[j] - rr_list[i].mass) - masstype[j] / 2.) > 0)	// The new particle takes the average mass further from mean than not adding it
		{
		  sum[j] -= rr_list[i].mass;	// Remove this particle from mass total; halfmassrad is still at previous particle value
		  break;
		}
	      else
		{
		  halfmassrad[j] = rr_list[i].r;
		  break;	// Leave loop when the half mass radius is reached
		}
	    }
	  halfmassrad[j] = rr_list[i].r;	// This will obviously be imprecise to one particle spacing, before real halfmass radius or after - depending on which is closer 
	}
    }

  /* compute cumulative mass */
  for(i = 1; i < num; i++)
    rr_list[i].mass = rr_list[i - 1].mass + rr_list[i].mass;

  for(i = num - 1, max = 0, maxrad = 0; i > 5; i--)
    if(rr_list[i].mass / rr_list[i].r > max)
      {
	max = rr_list[i].mass / rr_list[i].r;
	maxrad = rr_list[i].r;
      }

  *vmax = sqrt(All.G * max);
  *vmaxrad = maxrad;

  myfree(rr_list);
}

int subfind_compare_serial_candidates_boundlength(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->bound_length > ((struct cand_dat *) b)->bound_length)
    return -1;

  if(((struct cand_dat *) a)->bound_length < ((struct cand_dat *) b)->bound_length)
    return +1;

  if(((struct cand_dat *) a)->rank < ((struct cand_dat *) b)->rank)
    return -1;

  if(((struct cand_dat *) a)->rank > ((struct cand_dat *) b)->rank)
    return +1;

  return 0;
}

int subfind_compare_serial_candidates_rank(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->rank < ((struct cand_dat *) b)->rank)
    return -1;

  if(((struct cand_dat *) a)->rank > ((struct cand_dat *) b)->rank)
    return +1;

  if(((struct cand_dat *) a)->len > ((struct cand_dat *) b)->len)
    return -1;

  if(((struct cand_dat *) a)->len < ((struct cand_dat *) b)->len)
    return +1;

  return 0;
}

int subfind_compare_serial_candidates_subnr(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->subnr < ((struct cand_dat *) b)->subnr)
    return -1;

  if(((struct cand_dat *) a)->subnr > ((struct cand_dat *) b)->subnr)
    return +1;

  return 0;
}


#endif
