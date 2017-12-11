/* Using Dark Energy instead of a Cosmological constant can be archived by
 * replacing Lambda by Lambda * a^(-3*(1+w)) in the Hubble function.
 * So easy to see that w = -1 gives back a standard Cosmological Constant !
 * Also w = -1/3 gives Lambda / a^2 which then cancel within the Hubble
 * function and is then equal to the dynamics of a universe with Lambda = 0 !
 *
 * For a time varying w once has to replace Lambda * a^(-3*(1+w)) by
 * Lambda * exp(Integral(a,1,3*(1+w)/a))
 *
 * Dark Energy does not alter the powerspectrum of initial conditions.
 * To get the same cluster for various values or functions of w, once
 * has do assign a new redshift to the initial cond. to match the
 * linear growth factors, so g(z=0)/g(z_ini) == g_w(z=0)/g_w(z_ini^new)
 * Also the initial velocities field has to be scaled by 
 * (Hubble_w(z_ini^new)*Omega_w(z_ini^new)^0.6)/(Hubble(z_ini)*Omega(z_ini)^0.6)
 * where _w means the according functions including the terms for
 * Dark Energy.
 */

#ifdef DARKENERGY

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#ifdef TIMEDEPDE

#define ANZ_W_A_IN 1000
#define ANZ_W_A 10000

static MyFloat atab[ANZ_W_A_IN];

static MyFloat wtab[ANZ_W_A_IN];

static MyFloat intwtab[ANZ_W_A + 1];

/* set-up table with Exp(-Integral_a^1 (q+w(a'))/a' da'
 * needed in hubble function for time dependent w
 */

void fwa_init(void)
{
  int count = 0, i;

  char buf[200], buf1[200], buf2[200];

  FILE *fd;

  MyFloat a_first, w_first, a, w, sum;


  if(ThisTask == 0)
    {
      printf("initialize time dependent w ...\n");
      fflush(stdout);
    }

  if((fd = fopen(All.DarkEnergyFile, "r")))
    {
      if(ThisTask == 0)
	{
	  printf("\nreading w of a from file `%s'\n", All.DarkEnergyFile);
	  fflush(stdout);
	}
      atab[0] = -1.0;		/* we have to extrapolate wtab[0] later ! */
      count = 1;
      while(!feof(fd) && count < ANZ_W_A_IN)
	{
	  if(fgets(buf, 200, fd))
	    {
	      if(sscanf(buf, "%s%s", buf1, buf2) < 2)
		{
		  if(ThisTask == 0)
		    {
		      printf("Wrong syntax in file `%s', line %d\n", All.DarkEnergyFile, count);
		      fflush(stdout);
		    }
		  endrun(0);
		}
	      a = 1. / atof(buf1);
	      if(a == 0.0 && count == 1)
		count--;	/* w(0) present in file, so fill the first element ! */
	      atab[count] = a;
	      wtab[count] = atof(buf2);
	      count++;
	    }
	}
      fclose(fd);
      if(count >= ANZ_W_A_IN - 1)
	{
	  if(ThisTask == 0)
	    {
	      printf("File `%s' contains to many datapoints, increase ANZ_W_A_IN !\n", All.DarkEnergyFile);
	      fflush(stdout);
	    }
	  endrun(0);
	}
      if(count <= 2)
	{
	  if(ThisTask == 0)
	    {
	      printf("File `%s' has to less Data Points (%d) !\n", All.DarkEnergyFile, count);
	      fflush(stdout);
	    }
	  endrun(0);
	}

      if(atab[0] < 0.)		/* We still have to extrapolate w to a = 0 (w[0]) */
	{
	  atab[0] = 0.;
	  wtab[0] = wtab[1] - (wtab[2] - wtab[1]) / (atab[2] - atab[1]) * (atab[1] - atab[0]);
	}

/* Calculate w(1) if needed */
      if(atab[count - 1] < 1.)
	{
	  atab[count] = 1.0;
	  wtab[count] = wtab[count - 1] + (wtab[count - 1] - wtab[count - 2])
	    / (atab[count - 1] - atab[count - 2]) * (1. - atab[count - 1]);
/*            if(ThisTask ==0) 
              {
                printf("%d %f %f %f\n",count,atab[count-2],atab[count-1],atab[count]);
                printf("%d %f %f %f\n",count,wtab[count-2],wtab[count-1],wtab[count]);
              }*/
	  count++;
	}

/* Now calculated the integral (starting from last to first to save Time !
 * Explicit asume that a[0]=0. and a[count-1]=1. , which is enshured by   
 * The loading precedure !                                                  */

      count--;			/* Place count on last entry in table */
      intwtab[ANZ_W_A] = All.OmegaLambda;	/* Todays Omega_DarkEnergy */
      a_first = atab[count];	/* Startinv value should be 1.0 ! */
      w_first = wtab[count];	/* Starting value from table ! */
      sum = 0.0;		/* set int to 0.0 */
      for(i = ANZ_W_A - 1; i >= 1; i--)
	{
	  a = (FLOAT) i / (FLOAT) ANZ_W_A;
	  if(count > 1)		/* Still inside the table */
	    {
	      while(atab[count - 1] > a && count > 0)
		{
		  sum += 0.5 * ((1. + w_first) / a_first + (1. + wtab[count - 1]) / atab[count - 1])
		    * (a_first - atab[count - 1]);
		  count--;
		  a_first = atab[count];
		  w_first = wtab[count];
		}
	      w = w_first - (wtab[count] - wtab[count - 1]) / (atab[count] - atab[count - 1]) * (a_first - a);
	      sum += 0.5 * ((1. + w_first) / a_first + (1. + w) / a) * (a_first - a);
	      w_first = w;
	      a_first = a;
	    }
	  else
	    {
	      w = w_first - (wtab[count] - wtab[count - 1]) / (atab[count] - atab[count - 1]) * (a_first - a);
	      sum += 0.5 * ((1. + w_first) / a_first + (1. + w) / a) * (a_first - a);
	      w_first = w;
	      a_first = a;
	    }
	  intwtab[i] = All.OmegaLambda * exp(3. * sum);
	}
      /* artificially define value for a=0 */
      intwtab[0] = intwtab[1];
    }
  else
    {
      if(ThisTask == 0)
	{
	  printf("\nFile `%s' not found !\n", All.DarkEnergyFile);
	  fflush(stdout);
	}
      endrun(0);
    }

  if(ThisTask == 0)
    {
      printf("Integrating w(a) finisched.\n");
      fflush(stdout);
    }
}

/* This function tzhe integral w(a) therm for the actual time
 * needed in the hubble function.
 */
double INLINE_FUNC fwa(double a)
{
  int ai;

  double fwa = 0.0;

  ai = a * ANZ_W_A;
  if(ai >= ANZ_W_A)
    {
      fwa = intwtab[ANZ_W_A];
    }
  else
    {
      fwa = intwtab[ai] + (intwtab[ai + 1] - intwtab[ai]) * (a * (double) ANZ_W_A - (double) ai);
    }
/*   if(ThisTask==0) printf("%f %f %f %f\n",a,intwtab[ai],fwa,intwtab[ai+1]);*/
  return (fwa);
}

#endif



double DarkEnergy_a(double a)	/* only needed for comoving integration */
{
#ifdef TIMEDEPDE
  return fwa(a);
#else
  return (All.OmegaLambda * pow(a, -3. * (1 + All.DarkEnergyParam)));
#endif
}


double DarkEnergy_t(double Time)	/* only needed for physical integration */
{
  return All.DarkEnergyParam;
}

#endif
