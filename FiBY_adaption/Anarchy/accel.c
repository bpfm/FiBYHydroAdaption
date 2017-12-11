#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#ifdef BG_SFR
#include "bg_proto.h"
#endif

#ifdef BG_MOL_NETWORK
#include "bg_molecules.h"
#endif

/*! \file accel.c
 *  \brief driver routines to carry out force computation
 */

/*! This routine computes the accelerations for all active particles.  First, the gravitational forces are
 * computed. This also reconstructs the tree, if needed, otherwise the drift/kick operations have updated the
 * tree to make it fullu usable at the current time.
 *
 * If gas particles are presented, the `interior' of the local domain is determined. This region is guaranteed
 * to contain only particles local to the processor. This information will be used to reduce communication in
 * the hydro part.  The density for active SPH particles is computed next. If the number of neighbours should
 * be outside the allowed bounds, it will be readjusted by the function ensure_neighbours(), and for those
 * particle, the densities are recomputed accordingly. Finally, the hydrodynamical forces are added.
 */
void compute_accelerations(int mode)
{
  if(ThisTask == 0)
    {
      printf("Start force computation...\n");
      fflush(stdout);
    }

  CPU_Step[CPU_MISC] += measure_time();

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)
    {
      long_range_force();

      CPU_Step[CPU_MESH] += measure_time();
    }
#endif

#ifndef ONLY_PM
  gravity_tree();		/* computes gravity accel. */

  if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0)
    gravity_tree();		/* For the first timestep, we redo it
				 * to allow usage of relative opening
				 * criterion for consistent accuracy. 
				 */
#endif


  if(All.TotN_gas > 0)
    {
      if(ThisTask == 0)
	{
	  printf("Start density computation...\n");
	  fflush(stdout);
	}


      /***** density *****/
      density();		/* computes density, and pressure */



      /***** update smoothing lengths in tree *****/
      force_update_hmax();


      if(ThisTask == 0)
	{
	  printf("Start hydro-force computation...\n");
	  fflush(stdout);
	}


      /***** hydro forces *****/
      hydro_force();		/* adds hydrodynamical accelerations 
				   and computes du/dt  */


      /***** black holes *****/
#ifdef BLACK_HOLES
      blackhole_accretion();

      CPU_Step[CPU_BLACKHOLES] += measure_time();
#endif


#ifdef BG_SFR
      bg_enrich();		/* do chemical enrichment and eventually
				 * kinetic feedback in blue-gene model */

      CPU_Step[CPU_ENRICH] += measure_time();
#endif

#if defined(BLACK_HOLES) && defined (BH_THERMALFEEDBACK)
      bg_bh_thermal_feedback(); /* do agn thermal feedback */

      CPU_Step[CPU_BLACKHOLES] += measure_time();
#endif

#ifdef BG_MOL_NETWORK
#if defined(LW_BACKGROUND) || defined(LW_LOCAL)
      lw(0);
#endif
      bg_molecules();

      CPU_Step[CPU_MOLECULES] += measure_time();
#endif

#if defined(BG_SFR) || defined(BG_COOLING)
      cooling_and_starformation();	/* do radiative cooling and star formation */

      CPU_Step[CPU_COOLINGSFR] += measure_time();
#endif

#if defined(BG_SNII_THERMAL_FEEDBACK) // || (defined(BG_POPIII) && defined(BG_POPIII_THERMAL_FEEDBACK))
      bg_thermal_feedback();  /* do SNII/PISN thermal feedback */
#endif

#ifdef BG_THERMAL_FEEDBACK_MAKE_PARTICLES_ACTIVE
  reconstruct_timebins();
#endif

/* #if defined(BG_POPIII) && defined(BG_POPIII_BH_SEEDS) */
/*       popiii_black_holes(); */

/*       CPU_Step[CPU_BLACKHOLES] += measure_time(); */
/* #endif */

#if defined(BG_SNII_THERMAL_FEEDBACK) // || (defined(BG_POPIII) && defined(BG_POPIII_THERMAL_FEEDBACK)) || (defined(BLACK_HOLES) && defined (BH_THERMALFEEDBACK))
      sigvel();
#endif
    }



  if(ThisTask == 0)
    {
      printf("force computation done.\n");
      fflush(stdout);
    }
}
