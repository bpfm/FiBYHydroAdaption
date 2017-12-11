#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "proto.h"
#include "kernel.h"


/* This function returns the kernel value at distance u */

double kernel_wk(double u, double hinv3)
{
#ifdef SPH_KERNEL_GADGET
  if(u < 0.5)
    return hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
  else
    return hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
#endif

#ifdef SPH_KERNEL_CUBIC
  if(u < 0.5)
    return hinv3 *
        (KERNEL_COEFF_1 * (1.0 - u) * (1.0 - u) * (1.0 - u) -
            KERNEL_COEFF_2 * (0.5 - u) * (0.5 - u) * (0.5 - u));
  else
    return hinv3 *
        KERNEL_COEFF_1 * (1.0 - u) * (1.0 - u) * (1.0 - u);
#endif

#ifdef SPH_KERNEL_QUARTIC
  if(u < 0.2)
    return hinv3 *
        (KERNEL_COEFF_1 * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) -
            KERNEL_COEFF_2 * (0.6 - u) * (0.6 - u) * (0.6 - u) * (0.6 - u) +
            KERNEL_COEFF_3 * (0.2 - u) * (0.2 - u) * (0.2 - u) * (0.2 - u));
  else
    if(u >= 0.2 && u < 0.6)
      return hinv3 *
          (KERNEL_COEFF_1 * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) -
              KERNEL_COEFF_2 * (0.6 - u) * (0.6 - u) * (0.6 - u) * (0.6 - u));
    else
      return hinv3 *
          (KERNEL_COEFF_1 * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u));
#endif

#ifdef SPH_KERNEL_QUINTIC
  if(u < ONETHIRD)
    return hinv3 *
        (KERNEL_COEFF_1 * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) -
            KERNEL_COEFF_2 * (TWOTHIRDS - u) * (TWOTHIRDS - u) * (TWOTHIRDS - u) * (TWOTHIRDS - u) * (TWOTHIRDS - u) +
            KERNEL_COEFF_3 * (ONETHIRD - u) * (ONETHIRD - u) * (ONETHIRD - u) * (ONETHIRD - u) * (ONETHIRD - u));
  else
    if(u >= ONETHIRD && u < TWOTHIRDS)
      return hinv3 *
          (KERNEL_COEFF_1 * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) -
              KERNEL_COEFF_2 * (TWOTHIRDS - u) * (TWOTHIRDS - u) * (TWOTHIRDS - u) * (TWOTHIRDS - u) * (TWOTHIRDS - u));
    else
      return hinv3 *
          (KERNEL_COEFF_1 * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u));
#endif

#ifdef SPH_KERNEL_C2
    return KERNEL_COEFF_1 * hinv3 *
        (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) *
        (1.0 + 4 * u);
#endif

#ifdef SPH_KERNEL_C4
    return KERNEL_COEFF_1 * hinv3 *
        (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) *
        (3.0 + 18 * u + 35 * u * u);
#endif

#ifdef SPH_KERNEL_C6
  return KERNEL_COEFF_1 * hinv3 *
      (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) *
      (1.0 + 8 * u + 25 * u * u + 32 * u * u * u);
#endif
}


/* This function returns the kernel derivative at distance u */

double kernel_dwk(double u, double hinv4)
{
#ifdef SPH_KERNEL_GADGET
  if(u < 0.5)
    return hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
  else
    return hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
#endif

#ifdef SPH_KERNEL_CUBIC
  if(u < 0.5)
    return -hinv4 *
        (KERNEL_COEFF_3 * (1.0 - u) * (1.0 - u) -
            KERNEL_COEFF_4 * (0.5 - u) * (0.5 - u));
  else
    return -hinv4 *
        KERNEL_COEFF_3 * (1.0 - u) * (1.0 - u);
#endif

#ifdef SPH_KERNEL_QUARTIC
  if(u < 0.2)
    return -hinv4 *
        (KERNEL_COEFF_4 * (1.0 - u) * (1.0 - u) * (1.0 - u) -
            KERNEL_COEFF_5 * (0.6 - u) * (0.6 - u) * (0.6 - u) +
            KERNEL_COEFF_6 * (0.2 - u) * (0.2 - u) * (0.2 - u));
  else
    if(u >= 0.2 && u < 0.6)
      return -hinv4 *
          (KERNEL_COEFF_4 * (1.0 - u) * (1.0 - u) * (1.0 - u) -
              KERNEL_COEFF_5 * (0.6 - u) * (0.6 - u) * (0.6 - u));
    else
      return -hinv4 *
          (KERNEL_COEFF_4 * (1.0 - u) * (1.0 - u) * (1.0 - u));
#endif

#ifdef SPH_KERNEL_QUINTIC
  if(u < ONETHIRD)
    return -hinv4 *
        (KERNEL_COEFF_4 * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) -
            KERNEL_COEFF_5 * (TWOTHIRDS - u) * (TWOTHIRDS - u) * (TWOTHIRDS - u) * (TWOTHIRDS - u) +
            KERNEL_COEFF_6 * (ONETHIRD - u) * (ONETHIRD - u) * (ONETHIRD - u) * (ONETHIRD - u));
  else
    if(u >= ONETHIRD && u < TWOTHIRDS)
      return -hinv4 *
          (KERNEL_COEFF_4 * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) -
              KERNEL_COEFF_5 * (TWOTHIRDS - u) * (TWOTHIRDS - u) * (TWOTHIRDS - u) * (TWOTHIRDS - u));
    else
      return -hinv4 *
          (KERNEL_COEFF_4 * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u));
#endif

#ifdef SPH_KERNEL_C2
    return -20 * KERNEL_COEFF_1 * hinv4 *
        (1.0 - u) * (1.0 - u) * (1.0 - u) * u;
#endif

#ifdef SPH_KERNEL_C4
    return KERNEL_COEFF_1 * hinv4 *
        (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) *
        ((1.0 - u) * (18.0 + 70 * u) - (18.0 + 108 * u + 210 * u * u));
#endif

#ifdef SPH_KERNEL_C6
  return KERNEL_COEFF_1 * hinv4 *
      (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) *
      ((1.0 - u) * (8 + 50 * u + 96 * u * u) - 8 * (1.0 + 8 * u + 25 * u * u + 32 * u * u * u));
#endif
}

