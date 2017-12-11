#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef BG_EXTRA_ARRAYS

#include "allvars.h"
#include "proto.h"
#include "bg_proto.h"
#include "bg_vars.h"


void bg_compute_extra_arrays(void)
{
  int i;

  double temperature;

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0)
	if(SphP[i].OnEOS <= 0)
	  {
	    if(SphP[i].Entropy > SphP[i].MaximumEntropy)
	      {
		SphP[i].MaximumEntropy = SphP[i].Entropy;
		SphP[i].TimeMaximumEntropy = (MyFloat) All.Time;
	      }

	    temperature = bg_get_temperature(i);

	    if((MyFloat) temperature > SphP[i].MaximumTemperature)
	      {
		SphP[i].MaximumTemperature = (MyFloat) temperature;
		SphP[i].TimeMaximumTemperature = (MyFloat) All.Time;
	      }
	  }
    }
}

#endif
