#include <stdio.h>

#include "allvars.h"

#ifdef SUBFIND
#include "subfind.h"
#include "fof.h"


int Ncollective;

int MaxNsubgroups;

int Nsubgroups;

int TotNsubgroups;

/* #ifdef IGNOREFUZZ */
/* float sizelimit; */
/* int skipped; */
/* #endif */

struct subgroup_properties *SubGroup;


struct nearest_r2_data *R2Loc;

struct nearest_ngb_data *NgbLoc;


struct r2data *R2list;

struct nearest_ngb_data *NgbLoc;

struct nearest_r2_data *R2Loc;

double *Dist2list;

struct my_serial_index *Ind;

int *Rank;

int *SkipGrUnbind;

#endif
