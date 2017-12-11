#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
 
#include "allvars.h"
#include "proto.h"

/*! \file blackhole.c
 *  \brief routines for gas accretion onto black holes, and black hole mergers
 */

#ifdef BLACK_HOLES


static struct blackholedata_in
{
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyFloat Density;
  MyFloat Mdot;
  MyFloat Dt;
  MyFloat Hsml;
  MyFloat Mass;
  MyFloat BH_Mass;
  MyFloat Csnd;
#ifdef BH_THERMALFEEDBACK
  MyFloat Energy;
#endif
  MyIDType ID;
  int Index;
  int NodeList[NODELISTLENGTH];
}
 *BlackholeDataIn, *BlackholeDataGet;

static struct blackholedata_out
{
  MyLongDouble Mass;
  MyLongDouble BH_Mass;
  MyLongDouble AccretedMomentum[3];
#ifdef BH_THERMALFEEDBACK
  MyFloat Energy;
#endif
#if defined(REPOSITION_ON_POTMIN) && defined(COMPUTE_POTENTIAL_ENERGY)
  MyFloat BH_MinPotPos[3];
  MyFloat BH_MinPot;
#endif
}
 *BlackholeDataResult, *BlackholeDataOut;



static double hubble_a, ascale;

static double u_to_temp_fac;

static int N_gas_swallowed, N_BH_swallowed;

void blackhole_accretion(void)
{
  int i, j, k, n, bin;
  int ndone_flag, ndone;
  int ngrp, sendTask, recvTask, place, nexport, nimport, dummy;
  int Ntot_gas_swallowed, Ntot_BH_swallowed;
  double mdot, rho, bhvel, soundspeed, meddington, dt, mdot_in_msun_per_year;
  double mass_real, total_mass_real, medd, total_mdoteddington;
  double mass_holes, total_mass_holes, total_mdot;
  int bh_count;
  MyFloat accretion_prefactor;
  MPI_Status status;

  if(ThisTask == 0)
    {
      printf("Beginning black-hole accretion\n");
      fflush(stdout);
    }

  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1
    * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  CPU_Step[CPU_MISC] += measure_time();

  if(All.ComovingIntegrationOn)
    {
      ascale = All.Time;
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
    hubble_a = ascale = 1;

  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinActive[n])
	{
	  TimeBin_BH_mass[n] = 0;
	  TimeBin_BH_dynamicalmass[n] = 0;
	  TimeBin_BH_Mdot[n] = 0;
	  TimeBin_BH_Medd[n] = 0;
	}
    }

  /* Let's first compute the Mdot values */

  bh_count = 0;

  accretion_prefactor = 4. * M_PI * All.G * All.G * All.BlackHoleAccretionFactor;
  for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
    if(P[n].Type == 5)
      {

        bh_count++;
	mdot = 0;		/* if no accretion model is enabled, we have mdot=0 */

	rho = P[n].b1.BH_Density;

	bhvel = sqrt(pow(P[n].Vel[0] - P[n].b3.BH_SurroundingGasVel[0], 2) +
		     pow(P[n].Vel[1] - P[n].b3.BH_SurroundingGasVel[1], 2) +
		     pow(P[n].Vel[2] - P[n].b3.BH_SurroundingGasVel[2], 2));

	if(All.ComovingIntegrationOn)
	  {
	    bhvel /= All.Time;
	    rho /= pow(All.Time, 3);
	  }

	soundspeed = sqrt(GAMMA * P[n].b2.BH_Entropy * pow(rho, GAMMA_MINUS1));

	if (rho > All.PhysDensThresh) {
	  mdot = pow(rho / All.PhysDensThresh,All.BlackHoleAccretionSlope);
	} else {
	  mdot = 1.0;
	}

	mdot *= accretion_prefactor * P[n].BH_Mass * P[n].BH_Mass * rho  / pow((pow(soundspeed, 2) + pow(bhvel, 2)), 1.5);

	meddington = (4 * M_PI * GRAVITY * PROTONMASS / (All.BlackHoleRadiativeEfficiency * C * THOMPSON)) * P[n].BH_Mass
	  * All.UnitTime_in_s;

	if(mdot > All.BlackHoleEddingtonFactor * meddington)
	  mdot = All.BlackHoleEddingtonFactor * meddington;

	P[n].BH_Mdot = mdot;

#ifdef VERBOSE_LOGFILES
	  //Output black hole properties to logfile:
	if(P[n].BH_Mass > 0)
	  {
	      dt = (P[n].TimeBin ? (1 << P[n].TimeBin) : 0) * All.Timebase_interval / hubble_a;
	      fprintf(FdBlackHolesDetails, "%u %g %g %g %g %g %g %g %g %g\n",
		      P[n].ID, All.Time, P[n].Mass, P[n].BH_Mass, mdot, rho,soundspeed, bhvel, P[n].BH_Energy, dt);
	  }
#endif

	dt = (P[n].TimeBin ? (1 << P[n].TimeBin) : 0) * All.Timebase_interval / hubble_a;

#ifdef BH_THERMALFEEDBACK
	P[n].BH_Energy += All.BlackHoleFeedbackFactor * All.BlackHoleRadiativeEfficiency * mdot * dt *
	  pow(C / All.UnitVelocity_in_cm_per_s, 2);
#endif

#ifdef BH_DRAG
	/* add a drag force for the black-holes,
	   accounting for the accretion */
	double fac;

	if(P[n].BH_Mass > 0)
	  {
	    fac = meddington * dt / P[n].BH_Mass;

	    if(fac > 1)
	      fac = 1;

	    if(dt > 0)
	      for(k = 0; k < 3; k++)
		P[n].g.GravAccel[k] +=
		  -ascale * ascale * fac / dt * (P[n].Vel[k] - P[n].b3.BH_SurroundingGasVel[k]) / ascale;
	  }
#endif

	P[n].BH_Mass += (1 - All.BlackHoleRadiativeEfficiency) * P[n].BH_Mdot * dt;
      }

  /* Now let's invoke the functions that stochastically swallow gas
   * and deal with black hole mergers.
   */

  if(ThisTask == 0)
    {
      printf("Start swallowing of gas particles and black holes\n");
      fflush(stdout);
    }

  N_gas_swallowed = N_BH_swallowed = 0;

  /* allocate buffers to arrange communication */
  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct blackholedata_in) +
					     sizeof(struct blackholedata_out) +
					     sizemax(sizeof(struct blackholedata_in),
						     sizeof(struct blackholedata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));
  
  /** Determine which particles may be swalled by whom */

  i = FirstActiveParticle;	/* first particle for this task */

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
	  if(blackhole_evaluate(i, 0, &nexport, Send_count) < 0)
	    break;
      
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
      
      BlackholeDataGet = (struct blackholedata_in *) mymalloc(nimport * sizeof(struct blackholedata_in));
      BlackholeDataIn = (struct blackholedata_in *) mymalloc(nexport * sizeof(struct blackholedata_in));
      
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  
	  for(k = 0; k < 3; k++)
	    {
	      BlackholeDataIn[j].Pos[k] = P[place].Pos[k];
	      BlackholeDataIn[j].Vel[k] = P[place].Vel[k];
	    }
	  BlackholeDataIn[j].Hsml = PPP[place].Hsml;
	  BlackholeDataIn[j].Mass = P[place].Mass;
	  BlackholeDataIn[j].BH_Mass = P[place].BH_Mass;
	  BlackholeDataIn[j].Density = P[place].b1.BH_Density;
	  BlackholeDataIn[j].Mdot = P[place].BH_Mdot;
	  BlackholeDataIn[j].Csnd =
	    sqrt(GAMMA * P[place].b2.BH_Entropy *
		 pow(P[place].b1.BH_Density / (ascale * ascale * ascale), GAMMA_MINUS1));
	  BlackholeDataIn[j].Dt =
	    (P[place].TimeBin ? (1 << P[place].TimeBin) : 0) * All.Timebase_interval / hubble_a;
	  BlackholeDataIn[j].ID = P[place].ID;
	  
	  memcpy(BlackholeDataIn[j].NodeList,
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
		  MPI_Sendrecv(&BlackholeDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &BlackholeDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
		}
	    }
	}
      
      
      myfree(BlackholeDataIn);
      BlackholeDataResult = (struct blackholedata_out *) mymalloc(nimport * sizeof(struct blackholedata_out));
      BlackholeDataOut = (struct blackholedata_out *) mymalloc(nexport * sizeof(struct blackholedata_out));
      
      
      /* now do the particles that were sent to us */
      
      for(j = 0; j < nimport; j++)
	blackhole_evaluate(j, 1, &dummy, &dummy);
      
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
		  MPI_Sendrecv(&BlackholeDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholedata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B,
			       &BlackholeDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholedata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
		}
	    }
	  
	}
      
      /* add the result to the particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;
#if defined(REPOSITION_ON_POTMIN) && defined(COMPUTE_POTENTIAL_ENERGY)
	  if(P[place].BH_MinPot > BlackholeDataOut[j].BH_MinPot)
	    {
	      P[place].BH_MinPot = BlackholeDataOut[j].BH_MinPot;
	      for(k = 0; k < 3; k++)
		P[place].BH_MinPotPos[k] = BlackholeDataOut[j].BH_MinPotPos[k];
	    }
#endif
	   
	}
      
      myfree(BlackholeDataOut);
      myfree(BlackholeDataResult);
      myfree(BlackholeDataGet);
    }
  while(ndone < NTask);
  
  /* Now do the swallowing of particles */

  i = FirstActiveParticle;	/* first particle for this task */

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */

      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	  if(P[i].Type == 5) {
	      if(P[i].SwallowID == 0)
		  if(blackhole_evaluate_swallow(i, 0, &nexport, Send_count) < 0)
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

      BlackholeDataGet = (struct blackholedata_in *) mymalloc(nimport * sizeof(struct blackholedata_in));
      BlackholeDataIn = (struct blackholedata_in *) mymalloc(nexport * sizeof(struct blackholedata_in));

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    BlackholeDataIn[j].Pos[k] = P[place].Pos[k];

#ifdef BH_THERMALFEEDBACK
	  BlackholeDataIn[j].Energy = P[place].BH_Energy;
#endif
	  BlackholeDataIn[j].Hsml = PPP[place].Hsml;
	  BlackholeDataIn[j].BH_Mass = P[place].BH_Mass;
	  BlackholeDataIn[j].ID = P[place].ID;

	  memcpy(BlackholeDataIn[j].NodeList,
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
		  MPI_Sendrecv(&BlackholeDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &BlackholeDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
		}
	    }
	}


      myfree(BlackholeDataIn);
      BlackholeDataResult = (struct blackholedata_out *) mymalloc(nimport * sizeof(struct blackholedata_out));
      BlackholeDataOut = (struct blackholedata_out *) mymalloc(nexport * sizeof(struct blackholedata_out));


      /* now do the particles that were sent to us */

      for(j = 0; j < nimport; j++)
	blackhole_evaluate_swallow(j, 1, &dummy, &dummy);

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
		  MPI_Sendrecv(&BlackholeDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholedata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B,
			       &BlackholeDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholedata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
		}
	    }

	}

      /* add the result to the particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  P[place].b4.dBH_accreted_Mass += BlackholeDataOut[j].Mass;
	  P[place].b5.dBH_accreted_BHMass += BlackholeDataOut[j].BH_Mass;
#ifdef BH_THERMALFEEDBACK
	  P[place].b7.dBH_accreted_BHEnergy += BlackholeDataOut[j].Energy;
#endif
	  for(k = 0; k < 3; k++)
	    P[place].b6.dBH_accreted_momentum[k] += BlackholeDataOut[j].AccretedMomentum[k];
	}

      myfree(BlackholeDataOut);
      myfree(BlackholeDataResult);
      myfree(BlackholeDataGet);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);


  MPI_Reduce(&N_gas_swallowed, &Ntot_gas_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&N_BH_swallowed, &Ntot_BH_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("Accretion done: %d gas particles swallowed, %d BH particles swallowed\n",
	     Ntot_gas_swallowed, Ntot_BH_swallowed);
      fflush(stdout);
    }

#if defined(REPOSITION_ON_POTMIN) && defined(COMPUTE_POTENTIAL_ENERGY)
  for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
    if(P[n].Type == 5 && P[n].BH_Mass < 10.0*All.OrigGasMass)
      for(k = 0; k < 3; k++) {
	P[n].Pos[k] = P[n].BH_MinPotPos[k];
      }
#endif

  for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
    if(P[n].Type == 5)
      {
	P[n].b4.BH_accreted_Mass = FLT(P[n].b4.dBH_accreted_Mass);
	P[n].b5.BH_accreted_BHMass = FLT(P[n].b5.dBH_accreted_BHMass);
#ifdef BH_THERMALFEEDBACK
	P[n].b7.BH_accreted_BHEnergy = FLT(P[n].b7.dBH_accreted_BHEnergy);
#endif
	for(k = 0; k < 3; k++)
	  P[n].b6.BH_accreted_momentum[k] = FLT(P[n].b6.dBH_accreted_momentum[k]);

	if(P[n].b4.BH_accreted_Mass > 0)
	  {
	    for(k = 0; k < 3; k++)
	      P[n].Vel[k] =
		(P[n].Vel[k] * P[n].Mass + P[n].b6.BH_accreted_momentum[k]) /
		(P[n].Mass + P[n].b4.BH_accreted_Mass);
 
	    P[n].Mass += P[n].b4.BH_accreted_Mass;
	    P[n].BH_Mass += P[n].b5.BH_accreted_BHMass;
#ifdef BH_THERMALFEEDBACK
	    P[n].BH_Energy += P[n].b7.BH_accreted_BHEnergy;
#endif
	    P[n].b4.BH_accreted_Mass = 0;
	  }

	bin = P[n].TimeBin;
	TimeBin_BH_mass[bin] += P[n].BH_Mass;
	TimeBin_BH_dynamicalmass[bin] += P[n].Mass;
	TimeBin_BH_Mdot[bin] += P[n].BH_Mdot;
	TimeBin_BH_Medd[bin] += P[n].BH_Mass;
      }

  mdot = 0;
  mass_holes = 0;
  mass_real = 0;
  medd = 0;

  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinCount[bin])
      {
	mass_holes += TimeBin_BH_mass[bin];
	mass_real += TimeBin_BH_dynamicalmass[bin];
	mdot += TimeBin_BH_Mdot[bin];
	medd += TimeBin_BH_Medd[bin];
      }

  MPI_Reduce(&mass_holes, &total_mass_holes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mass_real, &total_mass_real, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mdot, &total_mdot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&medd, &total_mdoteddington, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      /* convert to solar masses per yr */
      mdot_in_msun_per_year =
	total_mdot * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      total_mdoteddington *= ((4 * M_PI * GRAVITY * PROTONMASS) /  (0.1 * C * THOMPSON)) * All.UnitTime_in_s;

      fprintf(FdBlackHoles, "%g %d %g %g %g %g %g\n",
	      All.Time, All.TotBHs, total_mass_holes, total_mdot, mdot_in_msun_per_year,
	      total_mass_real, total_mdoteddington);
      fflush(FdBlackHoles);
    }

  /* Remove zero mass particles */
  rearrange_particle_sequence();

#ifdef VERBOSE_LOGFILES
  fflush(FdBlackHolesDetails);
#endif

  CPU_Step[CPU_BLACKHOLES] += measure_time();
}






int blackhole_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
  int startnode, numngb, j, k, n, index, id, listindex = 0;
  MyDouble *pos, *velocity;
  MyFloat h_i, dt, mdot, rho, mass, bh_mass, csnd;
  double dx, dy, dz, h_i2, r2, r, u, hinv, hinv3, wk, vrel;
#ifdef SWALLOWGAS
  double p, w;
#endif
#ifdef BH_THERMALFEEDBACK
  double reservoir;
  double reservoir_used;
#endif
#if defined(REPOSITION_ON_POTMIN) && defined(COMPUTE_POTENTIAL_ENERGY)
  MyFloat minpotpos[3] = { 0, 0, 0 }, minpot = 1.0e30;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      rho = P[target].b1.BH_Density;
      mdot = P[target].BH_Mdot;
      dt = (P[target].TimeBin ? (1 << P[target].TimeBin) : 0) * All.Timebase_interval / hubble_a;
      h_i = PPP[target].Hsml;
      mass = P[target].Mass;
      bh_mass = P[target].BH_Mass;
#ifdef BH_THERMALFEEDBACK
      reservoir = P[target].BH_Energy;
      reservoir_used = 0.0;
#endif
      velocity = P[target].Vel;
      csnd =
	sqrt(GAMMA * P[target].b2.BH_Entropy *
	     pow(P[target].b1.BH_Density / (ascale * ascale * ascale), GAMMA_MINUS1));
      index = target;
      id = P[target].ID;
    }
  else
    {
      pos = BlackholeDataGet[target].Pos;
      rho = BlackholeDataGet[target].Density;
      mdot = BlackholeDataGet[target].Mdot;
      dt = BlackholeDataGet[target].Dt;
      h_i = BlackholeDataGet[target].Hsml;
      mass = BlackholeDataGet[target].Mass;
      bh_mass = BlackholeDataGet[target].BH_Mass;
#ifdef BH_THERMALFEEDBACK
      reservoir = BlackholeDataGet[target].Energy;
      reservoir_used = 0.0;
#endif      
      velocity = BlackholeDataGet[target].Vel;
      csnd = BlackholeDataGet[target].Csnd;
      index = BlackholeDataGet[target].Index;
      id = BlackholeDataGet[target].ID;
    }

  h_i2 = h_i * h_i;

  /* Now start the actual SPH computation for this particle */
  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = BlackholeDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }


  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb = ngb_treefind_blackhole(pos, h_i, target, &startnode, mode, nexport, nSend_local);

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = Ngblist[n];

	      if(P[j].Mass > 0)
		{
		  if(mass > 0)
		    {
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

		      if(r2 < h_i2)
			{

#if defined(REPOSITION_ON_POTMIN) && defined(COMPUTE_POTENTIAL_ENERGY)
			  /* if this option is switched on, we may also encounter dark matter particles or stars */
			  if(P[j].p.Potential < minpot)
			    {
			      /* compute relative velocities */

			      for(k = 0, vrel = 0; k < 3; k++)
				vrel += (P[j].Vel[k] - velocity[k]) * (P[j].Vel[k] - velocity[k]);

			      vrel = sqrt(vrel) / ascale;

			      if(vrel <= 0.25 * csnd)
				{
				  minpot = P[j].p.Potential;
				  for(k = 0; k < 3; k++)
				    minpotpos[k] = P[j].Pos[k];
				}
			    }
#endif

			  if(P[j].Type == 5)	/* we have a black hole merger */
			    {
			      if(r2 > 0)
				{
				  /* compute relative velocity of BHs */

				  for(k = 0, vrel = 0; k < 3; k++)
				    vrel += (P[j].Vel[k] - velocity[k]) * (P[j].Vel[k] - velocity[k]);

				  vrel = sqrt(vrel) / ascale;

				  if (pow(vrel,2) < 2 * All.G * mass / h_i) {

					 if(P[j].SwallowID < id)
						P[j].SwallowID = id;
				  } else {
				    //Want to merge but cannot because velocity is too high
				    //printf("Want to merge but can't %f %f\n",pow(vrel,2),All.G * mass / h_i);
				  }
				}
			    }
			  if(P[j].Type == 0)
			    {

			      /* here we have a gas particle */
			      r = sqrt(r2);
			      hinv = 1 / h_i;
#ifndef  TWODIMS
			      hinv3 = hinv * hinv * hinv;
#else
			      hinv3 = hinv * hinv / boxSize_Z;
#endif

			      u = r * hinv;

			      if(u < 0.5)
				wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
			      else
				wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

#ifdef SWALLOWGAS
			      /* compute accretion probability */

			      if((bh_mass - mass) > 0) {
				p = (bh_mass - mass) * wk / rho;
			      }else
				p = 0;

			      /* compute random number, uniform in [0,1] */
			      w = get_random_number(P[j].ID);
			      if(w < p)
				{
				  if(P[j].SwallowID < id) {
				    P[j].SwallowID = id;
				  }
				}
#endif			      
			    }
			}
		    }
		}
	    }
	}
      
      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = BlackholeDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }
  

  /* Now collect the result at the right place */
  if(mode == 0)
    {
#if defined(REPOSITION_ON_POTMIN) && defined(COMPUTE_POTENTIAL_ENERGY)
      for(k = 0; k < 3; k++)
	P[target].BH_MinPotPos[k] = minpotpos[k];
      P[target].BH_MinPot = minpot;
#endif
    }
  else
    {
#if defined(REPOSITION_ON_POTMIN) && defined(COMPUTE_POTENTIAL_ENERGY)
      for(k = 0; k < 3; k++)
	BlackholeDataResult[target].BH_MinPotPos[k] = minpotpos[k];
      BlackholeDataResult[target].BH_MinPot = minpot;
#endif
    }

  return 0;
}


int blackhole_evaluate_swallow(int target, int mode, int *nexport, int *nSend_local)
{
  int startnode, numngb, j, k, n, id, listindex = 0;
  MyLongDouble accreted_mass, accreted_BH_mass, accreted_momentum[3], accreted_BH_Energy;
  MyLongDouble accreted_cumlseedmass, accreted_cumlaccrmass;
  MyDouble *pos;
  MyFloat h_i, bh_mass;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h_i = PPP[target].Hsml;
      id = P[target].ID;
      bh_mass = P[target].BH_Mass;
    }
  else
    {
      pos = BlackholeDataGet[target].Pos;
      h_i = BlackholeDataGet[target].Hsml;
      id = BlackholeDataGet[target].ID;
      bh_mass = BlackholeDataGet[target].BH_Mass;
    }

  accreted_mass = 0;
  accreted_BH_mass = 0;
  accreted_BH_Energy = 0;
  accreted_momentum[0] = accreted_momentum[1] = accreted_momentum[2] = 0;
  accreted_cumlseedmass = 0;
  accreted_cumlaccrmass = 0;

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = BlackholeDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb = ngb_treefind_blackhole(pos, h_i, target, &startnode, mode, nexport, nSend_local);

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = Ngblist[n];

	      if(P[j].SwallowID == id)
		{
		  if(P[j].Type == 5)	/* we have a black hole merger */
		    {
		      accreted_mass += FLT(P[j].Mass);
		      accreted_BH_mass += FLT(P[j].BH_Mass);
#ifdef BH_THERMALFEEDBACK
		      accreted_BH_Energy += FLT(P[j].BH_Energy);
#endif
		      for(k = 0; k < 3; k++)
				  accreted_momentum[k] += FLT(P[j].Mass * P[j].Vel[k]);

		      P[j].Mass = 0;
		      P[j].BH_Mass = 0;

		      N_BH_swallowed++;
		    }
		}

	      if(P[j].Type == 0)	/* we have a gas particle */
		{
		  if(P[j].SwallowID == id)
		    {
		      accreted_mass += FLT(P[j].Mass * (1.0 - All.BlackHoleFeedbackFactor));
		      for(k = 0; k < 3; k++)
			accreted_momentum[k] += FLT(P[j].Mass * P[j].Vel[k]);
		      
		      P[j].Mass = 0;
		      
		      N_gas_swallowed++;
		    }
		}
	    }
	}
      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = BlackholeDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  /* Now collect the result at the right place */
  if(mode == 0)
    {
      P[target].b4.dBH_accreted_Mass = accreted_mass;
      P[target].b5.dBH_accreted_BHMass = accreted_BH_mass;
      P[target].b7.dBH_accreted_BHEnergy = accreted_BH_Energy;
      for(k = 0; k < 3; k++)
	P[target].b6.dBH_accreted_momentum[k] = accreted_momentum[k];
    }
  else
    {
      BlackholeDataResult[target].Mass = accreted_mass;
      BlackholeDataResult[target].BH_Mass = accreted_BH_mass;
#ifdef BH_THERMALFEEDBACK
      BlackholeDataResult[target].Energy = accreted_BH_Energy;
#endif
      for(k = 0; k < 3; k++)
	BlackholeDataResult[target].AccretedMomentum[k] = accreted_momentum[k];
    }


  return 0;
}




int ngb_treefind_blackhole(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			   int *nexport, int *nsend_local)
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

#ifndef REPOSITION_ON_POTMIN
	  if(P[p].Type != 0 && P[p].Type != 5)
	    continue;
#endif
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
		endrun(12312);

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


#endif
