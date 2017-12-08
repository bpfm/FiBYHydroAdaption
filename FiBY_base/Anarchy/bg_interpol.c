#if defined(BG_COOLING) || defined(BG_STELLAR_EVOLUTION)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "proto.h"
#include "bg_proto.h"
#include "bg_cooling.h"

#define EPS 1.e-4



/*
 * ----------------------------------------------------------------------
 * This routine returns the position i of a value x in a 1D table and the
 * displacement dx needed for the interpolation.  The table is assumed to
 * be evenly spaced.
 * ----------------------------------------------------------------------
 */

void get_index_1d(float *table, int ntable, double x, int *i, float *dx)
{
  float dxm1;

  dxm1 = (float) (ntable - 1) / (table[ntable - 1] - table[0]);

  if((float) x <= table[0] + EPS)
    {
      *i = 0;
      *dx = 0;
    }
  else if((float) x >= table[ntable - 1] - EPS)
    {
      *i = ntable - 2;
      *dx = 1;
    }
  else
    {
      *i = (int) floor(((float) x - table[0]) * dxm1);
      *dx = ((float) x - table[*i]) * dxm1;
    }
}


/*
 * ----------------------------------------------------------------------
 * This routine performs a linear interpolation
 * ----------------------------------------------------------------------
 */

float interpol_1d(float *table, int i, float dx)
{
  float result;

  result = (1 - dx) * table[i] + dx * table[i + 1];

  return result;
}


/*
 * ----------------------------------------------------------------------
 * This routine performs a linear interpolation
 * ----------------------------------------------------------------------
 */

double interpol_1d_dbl(double *table, int i, float dx)
{
  double result;

  result = (1 - dx) * table[i] + dx * table[i + 1];

  return result;
}


/*
 * ----------------------------------------------------------------------
 * This routine performs a bi-linear interpolation
 * ----------------------------------------------------------------------
 */

float interpol_2d(float **table, int i, int j, float dx, float dy)
{
  float result;

  result =
    (1 - dx) * (1 - dy) * table[i][j] +
    (1 - dx) * dy * table[i][j + 1] +
    dx * (1 - dy) * table[i + 1][j] +
    dx * dy * table[i + 1][j + 1];

  return result;
}


/*
 * ----------------------------------------------------------------------
 * This routine performs a bi-linear interpolation
 * ----------------------------------------------------------------------
 */

double interpol_2d_dbl(double **table, int i, int j, double dx, double dy)
{
  double result;

  result =
    (1 - dx) * (1 - dy) * table[i][j] +
    (1 - dx) * dy * table[i][j + 1] +
    dx * (1 - dy) * table[i + 1][j] +
    dx * dy * table[i + 1][j + 1];

  return result;
}


/*
 * ----------------------------------------------------------------------
 * This routine performs a tri-linear interpolation
 * ----------------------------------------------------------------------
 */

float interpol_3d(float ***table, int i, int j, int k, float dx, float dy, float dz)
{
  float result;

  result =
    (1 - dx) * (1 - dy) * (1 - dz) * table[i][j][k] +
    (1 - dx) * (1 - dy) * dz * table[i][j][k + 1] +
    (1 - dx) * dy * (1 - dz) * table[i][j + 1][k] +
    (1 - dx) * dy * dz * table[i][j + 1][k + 1] +
    dx * (1 - dy) * (1 - dz) * table[i + 1][j][k] +
    dx * (1 - dy) * dz * table[i + 1][j][k + 1] +
    dx * dy * (1 - dz) * table[i + 1][j + 1][k] +
    dx * dy * dz * table[i + 1][j + 1][k + 1];

  return result;
}

float interpol_3d_metals(float ***table, int i, int j, int k, float dy, float dz)
{
  float result;

  result =
    (1 - dy) * (1 - dz) * table[i][j][k] +
    (1 - dy) * dz * table[i][j][k + 1] +
    dy * (1 - dz) * table[i][j + 1][k] +
    dy * dz * table[i][j + 1][k + 1];

  return result;
}


/*
 * ----------------------------------------------------------------------
 * This routine performs a quadri-linear interpolation
 * ----------------------------------------------------------------------
 */

float interpol_4d(float ****table, int i, int ip, int j, int jp, int k, int kp,
		  int l, int lp, float dx, float dy, float dz, float dw)
{
  float result;

  result =
    (1 - dx) * (1 - dy) * (1 - dz) * (1 - dw) * table[i][j][k][l] +
    (1 - dx) * (1 - dy) * (1 - dz) * dw * table[i][j][k][lp] +
    (1 - dx) * (1 - dy) * dz * (1 - dw) * table[i][j][kp][l] +
    (1 - dx) * (1 - dy) * dz * dw * table[i][j][kp][lp] +
    (1 - dx) * dy * (1 - dz) * (1 - dw) * table[i][jp][k][l] +
    (1 - dx) * dy * (1 - dz) * dw * table[i][jp][k][lp] +
    (1 - dx) * dy * dz * (1 - dw) * table[i][jp][kp][l] +
    (1 - dx) * dy * dz * dw * table[i][jp][kp][lp] +
    dx * (1 - dy) * (1 - dz) * (1 - dw) * table[ip][j][k][l] +
    dx * (1 - dy) * (1 - dz) * dw * table[ip][j][k][lp] +
    dx * (1 - dy) * dz * (1 - dw) * table[ip][j][kp][l] +
    dx * (1 - dy) * dz * dw * table[ip][j][kp][lp] +
    dx * dy * (1 - dz) * (1 - dw) * table[ip][jp][k][l] +
    dx * dy * (1 - dz) * dw * table[ip][jp][k][lp] +
    dx * dy * dz * (1 - dw) * table[ip][jp][kp][l] +
    dx * dy * dz * dw * table[ip][jp][kp][lp];

  return result;
}

#endif
