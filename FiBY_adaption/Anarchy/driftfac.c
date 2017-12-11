#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"


static double logTimeBegin;

static double logTimeMax;


double drift_integ(double a, void *param)
{
  double h;

  h = All.Omega0 / (a * a * a) + (1 - All.Omega0 - All.OmegaLambda) / (a * a)
#ifdef DARKENERGY
    + DarkEnergy_a(a);
#else
    + All.OmegaLambda;
#endif
  h = All.Hubble * sqrt(h);

  return 1 / (h * a * a * a);
}

double gravkick_integ(double a, void *param)
{
  double h;

  h = All.Omega0 / (a * a * a) + (1 - All.Omega0 - All.OmegaLambda) / (a * a)
#ifdef DARKENERGY
    + DarkEnergy_a(a);
#else
    + All.OmegaLambda;
#endif
  h = All.Hubble * sqrt(h);

  return 1 / (h * a * a);
}


double hydrokick_integ(double a, void *param)
{
  double h;

  h = All.Omega0 / (a * a * a) + (1 - All.Omega0 - All.OmegaLambda) / (a * a)
#ifdef DARKENERGY
    + DarkEnergy_a(a);
#else
    + All.OmegaLambda;
#endif
  h = All.Hubble * sqrt(h);

  return 1 / (h * pow(a, 3 * GAMMA_MINUS1) * a);
}

double growthfactor_integ(double a, void *param)
{
  double s;

  s = All.Omega0 + (1 - All.Omega0 - All.OmegaLambda) * a
#ifdef DARKENERGY
    + DarkEnergy_a(a) * a * a * a;
#else
    + All.OmegaLambda * a * a * a;
#endif
  s = sqrt(s);

  return pow(sqrt(a) / s, 3);
}


void init_drift_table(void)
{
#define WORKSIZE 100000
  int i;

  double result, abserr;

  gsl_function F;

  gsl_integration_workspace *workspace;

  logTimeBegin = log(All.TimeBegin);
  logTimeMax = log(All.TimeMax);

  workspace = gsl_integration_workspace_alloc(WORKSIZE);

  for(i = 0; i < DRIFT_TABLE_LENGTH; i++)
    {
      F.function = &drift_integ;
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), All.Hubble,	/* note: absolute error just a dummy */
			  1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      DriftTable[i] = result;


      F.function = &gravkick_integ;
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), All.Hubble,	/* note: absolute error just a dummy */
			  1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      GravKickTable[i] = result;


      F.function = &hydrokick_integ;
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), All.Hubble,	/* note: absolute error just a dummy */
			  1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      HydroKickTable[i] = result;
    }

  gsl_integration_workspace_free(workspace);
}


/*! This function integrates the cosmological prefactor for a drift
 *   step between time0 and time1. The value returned is
 *  \f[ \int_{a_0}^{a_1} \frac{{\rm d}a}{H(a)}
 *  \f]
 *  
 *  A lookup-table is used for reasons of speed. 
 */
double get_drift_factor(int time0, int time1)
{
  double a1, a2, df1, df2, u1, u2;

  int i1, i2;

  static int last_time0 = -1, last_time1 = -1;

  static double last_value;

  if(time0 == last_time0 && time1 == last_time1)
    return last_value;

  /* note: will only be called for cosmological integration */

  a1 = logTimeBegin + time0 * All.Timebase_interval;
  a2 = logTimeBegin + time1 * All.Timebase_interval;

  u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int) u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  if(i1 <= 1)
    df1 = u1 * DriftTable[0];
  else
    df1 = DriftTable[i1 - 1] + (DriftTable[i1] - DriftTable[i1 - 1]) * (u1 - i1);


  u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i2 = (int) u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  if(i2 <= 1)
    df2 = u2 * DriftTable[0];
  else
    df2 = DriftTable[i2 - 1] + (DriftTable[i2] - DriftTable[i2 - 1]) * (u2 - i2);

  last_time0 = time0;
  last_time1 = time1;

  return last_value = (df2 - df1);
}


double get_gravkick_factor(int time0, int time1)
{
  double a1, a2, df1, df2, u1, u2;

  int i1, i2;

  static int last_time0 = -1, last_time1 = -1;

  static double last_value;

  if(time0 == last_time0 && time1 == last_time1)
    return last_value;

  /* note: will only be called for cosmological integration */

  a1 = logTimeBegin + time0 * All.Timebase_interval;
  a2 = logTimeBegin + time1 * All.Timebase_interval;

  u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int) u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  if(i1 <= 1)
    df1 = u1 * GravKickTable[0];
  else
    df1 = GravKickTable[i1 - 1] + (GravKickTable[i1] - GravKickTable[i1 - 1]) * (u1 - i1);


  u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i2 = (int) u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  if(i2 <= 1)
    df2 = u2 * GravKickTable[0];
  else
    df2 = GravKickTable[i2 - 1] + (GravKickTable[i2] - GravKickTable[i2 - 1]) * (u2 - i2);

  last_time0 = time0;
  last_time1 = time1;

  return last_value = (df2 - df1);
}

double get_hydrokick_factor(int time0, int time1)
{
  double a1, a2, df1, df2, u1, u2;

  int i1, i2;

  static int last_time0 = -1, last_time1 = -1;

  static double last_value;

  if(time0 == last_time0 && time1 == last_time1)
    return last_value;

  /* note: will only be called for cosmological integration */

  a1 = logTimeBegin + time0 * All.Timebase_interval;
  a2 = logTimeBegin + time1 * All.Timebase_interval;

  u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int) u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  if(i1 <= 1)
    df1 = u1 * HydroKickTable[0];
  else
    df1 = HydroKickTable[i1 - 1] + (HydroKickTable[i1] - HydroKickTable[i1 - 1]) * (u1 - i1);


  u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i2 = (int) u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  if(i2 <= 1)
    df2 = u2 * HydroKickTable[0];
  else
    df2 = HydroKickTable[i2 - 1] + (HydroKickTable[i2] - HydroKickTable[i2 - 1]) * (u2 - i2);

  last_time0 = time0;
  last_time1 = time1;

  return last_value = (df2 - df1);
}
