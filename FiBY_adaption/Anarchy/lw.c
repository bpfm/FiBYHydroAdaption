#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#if defined(LW_BACKGROUND) || defined(LW_LOCAL)

#include "allvars.h"
#include "proto.h"
#include "bg_proto.h"


#define LW_STELLAR_AGE_GYR 0.005


static struct lwdata_in
{
  MyDouble Pos[3];
  MyDouble Mass;
  MyFloat Metallicity;
}
 *LWDataIn, *LWDataGet;


#ifdef LW_LOCAL
int TotalCandidates;
#endif

#ifdef LW_BACKGROUND
double popII_TotSFR, popIII_TotSFR;

MyFloat set_dissociating_background_coefficients_h2(void);
MyFloat set_dissociating_background_coefficients_hm(void);

MyFloat get_h2_shielding_factor(int i);
#endif


int ngb_treefind_variable_lw(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			     int *nexport, int *nsend_local);


void lw(int mode)
{
  int i, j, k, ngrp;
  int sendTask, recvTask;

  MyFloat f_shield;


  /*
   * here compute the contribution of the background Lyman-Werner radiation
   * to the dissociation coefficients k_diss,H2 and k_diss,H-
   */

#ifdef LW_BACKGROUND
  int bin;
  MyFloat LWBackground_Kdiss_H2, LWBackground_Kdiss_Hm;
  double popii_sfr, popiii_sfr;

  /* POPIII and POPII global SFR */
  if(mode == 0)
    {
      for(bin = 0, popii_sfr = 0; bin < TIMEBINS; bin++)
	if(TimeBinCount[bin])
	  popii_sfr += TimeBinSfr[bin];
#ifdef BG_POPIII
      for(bin = 0, popiii_sfr = 0; bin < TIMEBINS; bin++)
	if(TimeBinCount[bin])
	  popiii_sfr += TimeBinPOPIIISfr[bin];
#endif

      popii_sfr -= popiii_sfr;

      MPI_Allreduce(&popii_sfr, &popII_TotSFR, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#ifdef BG_POPIII
      MPI_Allreduce(&popiii_sfr, &popIII_TotSFR, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

      /* convert to a SFR per comoving Mpc^3 */
      popII_TotSFR /= pow(All.BoxSize * All.UnitLength_in_cm / All.HubbleParam / CM_PER_MPC, 3);
      popIII_TotSFR /= pow(All.BoxSize * All.UnitLength_in_cm / All.HubbleParam / CM_PER_MPC, 3);

      LWBackground_Kdiss_H2 = set_dissociating_background_coefficients_h2();
      LWBackground_Kdiss_Hm = set_dissociating_background_coefficients_hm();
    }
#endif /* LW_BACKGROUND */


  /*
   * here compute the contribution of the local Lyman-Werner radiation
   * to the dissociation coefficients k_diss,H2 and k_diss,H-
   */

#ifdef LW_LOCAL
  int local_send_count, total_send_count;
  double age_of_star_in_Gyr;

  /* allocate to the total number of stars, eventually less memory is needed */
  LWDataIn = (struct lwdata_in *) mymalloc_msg(N_star * sizeof(struct lwdata_in), "lwdata_in");

  /* select candidates */
  for(i = 0, local_send_count = 0; i < N_star; i++)
    {
      age_of_star_in_Gyr = bg_get_elapsed_time(StarP[i].StarBirthTime, All.Time, 1);

      if(age_of_star_in_Gyr <= LW_STELLAR_AGE_GYR)
	{
	  j = StarP[i].PID; /* the corresponding P index */

	  for(k = 0; k < 3; k++)
	    LWDataIn[local_send_count].Pos[k] = P[j].Pos[k];

	  LWDataIn[local_send_count].Mass = StarP[i].InitialMass;

#ifdef BG_POPIII
#ifdef BG_METALSMOOTHING
	  LWDataIn[local_send_count].Metallicity = StarP[i].MetallicitySmoothed;
#else
	  LWDataIn[local_send_count].Metallicity = StarP[i].Metallicity;
#endif
#endif
	  local_send_count++;
	}
    }

  /* compute the total number of candidates */
  MPI_Allreduce(&local_send_count, &TotalCandidates, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(TotalCandidates > 0)
    {
      /* distribute the number of candidates in each processor */
      MPI_Allgather(&local_send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

      /* just keep the 'standard' Gadget way of doing it */
      for(i = 0, Recv_offset[0] = 0; i < NTask; i++)
	{
	  Send_count[i] = local_send_count;
	  
	  if(i > 0)
	    {
	      Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
	    }
	}

      /* allocate receive array */
      LWDataGet = (struct lwdata_in *) mymalloc_msg((TotalCandidates + 1) * sizeof(struct lwdata_in), "lwdata_in");

      /* copy local candidates to global array */
      for(i = 0, j = Recv_offset[ThisTask]; i < local_send_count; i++)
	{
	  for(k = 0; k < 3; k++)
	    LWDataGet[j].Pos[k] = LWDataIn[i].Pos[k];
	  
	  LWDataGet[j].Mass = LWDataIn[i].Mass;
	  LWDataGet[j].Metallicity = LWDataIn[i].Metallicity;
	  
	  j++;
	}

      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0) /* get the particles */
	      MPI_Sendrecv(LWDataIn,
			   Send_count[recvTask] * sizeof(struct lwdata_in), MPI_BYTE,
			   recvTask, TAG_HYDRO_A,
			   &LWDataGet[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct lwdata_in), MPI_BYTE,
			   recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

      if(mode == 0)
	lw_evaluate();
      else
	lw_evaluate_global();

      myfree(LWDataGet);
    }

  myfree(LWDataIn);
#endif /* LW_LOCAL */


  /*
   * now account for shielding and eventually add the background
   * contribution when it is larger than the local
   */

#if defined(LW_BACKGROUND) || defined(LW_LOCAL)
  if(mode == 0)
    {
      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	if(P[i].Type == 0 && SphP[i].OnEOS != 1)
	  {
	    f_shield = get_h2_shielding_factor(i);

#ifdef LW_BACKGROUND
	    SphP[i].Kdiss_H2 += LWBackground_Kdiss_H2;
	    SphP[i].Kdiss_Hm += LWBackground_Kdiss_Hm;
#endif

	    SphP[i].Kdiss_H2 *= f_shield;
	  }
    }
#endif
}


#ifdef LW_LOCAL
void lw_evaluate(void)
{
  int i, j;
  MyDouble *pos;
  double dx, dy, dz;
  double r2, r2inv;

  MyFloat lwlocal_kdiss_h2, lwlocal_kdiss_hm;
  MyFloat lwlocal_radiation_popii, lwlocal_radiation_popiii;
  double mass, metallicity;


  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
	SphP[i].Kdiss_H2 = 0;
	SphP[i].Kdiss_Hm = 0;

	SphP[i].LWRadiation_popII = 0;
	SphP[i].LWRadiation_popIII = 0;

	if(SphP[i].OnEOS != 1)
	  {
	    pos = P[i].Pos;

	    lwlocal_kdiss_h2 = 0;
	    lwlocal_kdiss_hm = 0;

	    lwlocal_radiation_popii = 0;
	    lwlocal_radiation_popiii = 0;

	    for(j = 0; j < TotalCandidates; j++)
	      {
		dx = pos[0] - LWDataGet[j].Pos[0];
		dy = pos[1] - LWDataGet[j].Pos[1];
		dz = pos[2] - LWDataGet[j].Pos[2];
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

		if(r2 > 0)
		  {
		    mass = LWDataGet[j].Mass *
		      All.UnitMass_in_g / SOLAR_MASS / All.HubbleParam; /* in solar masses */

		    metallicity = LWDataGet[j].Metallicity;

		    r2inv = (1 / r2) *
		      All.HubbleParam * CM_PER_MPC / (1000 * All.Time * All.UnitLength_in_cm) *
		      All.HubbleParam * CM_PER_MPC / (1000 * All.Time * All.UnitLength_in_cm); /* in physical kpc */

#ifdef BG_POPIII /* if popIII stars are present, discriminate by metallicity */
		    if(metallicity >= All.POPIII_MetallicityThreshold)
		      {
#endif
			/* popII contribution */
			lwlocal_kdiss_h2 += (MyFloat) (0.67 * mass * r2inv);
			lwlocal_kdiss_hm += (MyFloat) (4e3 * mass * r2inv);

			lwlocal_radiation_popii += (MyFloat) (0.002 * 1.5 * mass * r2inv);
#ifdef BG_POPIII
		      }
		    else
		      {
			/* popIII contribution */
			lwlocal_kdiss_h2 += (MyFloat) (mass * r2inv);
			lwlocal_kdiss_hm += (MyFloat) (mass * r2inv);

			lwlocal_radiation_popiii += (MyFloat) (0.002 * 7.5 * mass * r2inv);
		      }
#endif
		  }
	      }

	    SphP[i].Kdiss_H2 = 1.36e-14 * lwlocal_kdiss_h2;
	    SphP[i].Kdiss_Hm = 1.50e-12 * lwlocal_kdiss_hm;

	    SphP[i].LWRadiation_popII = lwlocal_radiation_popii;
	    SphP[i].LWRadiation_popIII = lwlocal_radiation_popiii;
	  }
      }
}

void lw_evaluate_global(void)
{
  int i, j;
  MyDouble *pos;
  double dx, dy, dz;
  double r2, r2inv;

  MyFloat lwlocal_radiation_popii, lwlocal_radiation_popiii;
  double mass, metallicity;


  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 0)
      {
	SphP[i].LWRadiation_popII = 0;
	SphP[i].LWRadiation_popIII = 0;

	if(SphP[i].OnEOS != 1)
	  {
	    pos = P[i].Pos;

	    lwlocal_radiation_popii = 0;
	    lwlocal_radiation_popiii = 0;

	    for(j = 0; j < TotalCandidates; j++)
	      {
		dx = pos[0] - LWDataGet[j].Pos[0];
		dy = pos[1] - LWDataGet[j].Pos[1];
		dz = pos[2] - LWDataGet[j].Pos[2];
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

		if(r2 > 0)
		  {
		    mass = LWDataGet[j].Mass *
		      All.UnitMass_in_g / SOLAR_MASS / All.HubbleParam; /* in solar masses */

		    metallicity = LWDataGet[j].Metallicity;

		    r2inv = (1 / r2) *
		      All.HubbleParam * CM_PER_MPC / (1000 * All.Time * All.UnitLength_in_cm) *
		      All.HubbleParam * CM_PER_MPC / (1000 * All.Time * All.UnitLength_in_cm); /* in physical kpc */

#ifdef BG_POPIII /* if popIII stars are present, discriminate by metallicity */
		    if(metallicity >= All.POPIII_MetallicityThreshold)
		      {
#endif
			/* popII contribution */
			lwlocal_radiation_popii += (MyFloat) (0.002 * 1.5 * mass * r2inv);
#ifdef BG_POPIII
		      }
		    else
		      {
			/* popIII contribution */
			lwlocal_radiation_popiii += (MyFloat) (0.002 * 7.5 * mass * r2inv);
		      }
#endif
		  }
	      }

	    SphP[i].LWRadiation_popII = lwlocal_radiation_popii;
	    SphP[i].LWRadiation_popIII = lwlocal_radiation_popiii;
	  }
      }
}
#endif /* LW_LOCAL */

MyFloat get_h2_shielding_factor(int i)
{
  double h2_shielding_factor, h2_column, x, b5;
  double rho, temp, x_H_tot, n_H_tot, n_H2;
  double a3inv;

  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  /* convert to physical density */
  rho = SphP[i].d.Density * a3inv;
  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

  /* get temperature of the gas */
  temp = bg_get_temperature(i);

  /* get total mass fractions of hydrogen */
  x_H_tot = SphP[i].Metals[element_index("Hydrogen")] / P[i].Mass;

  /* total hydrogen number density */
  n_H_tot = x_H_tot * rho / PROTONMASS;

  /* H2 number density */
  n_H2 = SphP[i].x_H2 * rho / (2 * PROTONMASS);

/*   /\* H2 column density - Jarrett's recipe *\/ */
/*   h2_column = 2.0e19 * n_H2 / n_H_tot * sqrt(n_H_tot * temp); /\* N_H2  *\/ */

/*   /\* H2 shielding factor *\/ */
/*   if(h2_column < 1.0e14) */
/*     h2_shielding_factor = 1; */
/*   else */
/*     h2_shielding_factor = pow(1.0e14 / h2_column, 0.75); */

  /* from Wolcott-Green et al 2011 */
  x = 4.0e4 * n_H2 / n_H_tot * sqrt(n_H_tot * temp); /* N_H2 / [5e14 cm^(-2)]  */

  b5 = 2.9 * sqrt(temp * 1e-3);

  h2_shielding_factor = 0.965 / pow(1 + x / b5, 1.1) +
    0.035 / sqrt(1 + x) * exp(-8.5e-4 * sqrt(1 + x));

  return (MyFloat) h2_shielding_factor;
}

#ifdef LW_BACKGROUND
MyFloat set_dissociating_background_coefficients_h2(void)
{
  double ascale;

  if(All.ComovingIntegrationOn)
    ascale = All.Time;
  else
    ascale = 1;

  return (MyFloat) ( 1.4e-09 * (popIII_TotSFR + 0.67 * popII_TotSFR) * pow(16.0 * ascale, -3) );
}

MyFloat set_dissociating_background_coefficients_hm(void)
{
  double ascale;

  if(All.ComovingIntegrationOn)
    ascale = All.Time;
  else
    ascale = 1;

  return (MyFloat) ( 1.5e-08 * (popIII_TotSFR + 4e3 * popII_TotSFR) * pow(16.0 * ascale, -3) );
}
#endif /* LW_BACKGROUND */

#endif /* LW */
