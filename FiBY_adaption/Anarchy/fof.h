#ifndef FOF_H
#define FOF_H

#include "allvars.h"

int fof_compare_FOF_PList_MinID(const void *a, const void *b);
int fof_compare_FOF_GList_MinID(const void *a, const void *b);
int fof_compare_FOF_GList_MinIDTask(const void *a, const void *b);
int fof_compare_FOF_GList_LocCountTaskDiffMinID(const void *a, const void *b);
int fof_compare_FOF_GList_ExtCountMinID(const void *a, const void *b);
int fof_compare_Group_GrNr(const void *a, const void *b);
int fof_compare_Group_MinIDTask(const void *a, const void *b);
int fof_compare_Group_MinID(const void *a, const void *b);
int fof_compare_ID_list_GrNrID(const void *a, const void *b);
int fof_compare_ID_list_GrNrTypeID(const void *a, const void *b);
int fof_compare_Group_MinIDTask_MinID(const void *a, const void *b);
int fof_compare_Group_Len(const void *a, const void *b);
int fof_compare_Star_list_GrNrID(const void *a, const void *b);

void fof_compute_group_properties(int gr, int start, int len);
void fof_exchange_group_data(void);
void fof_finish_group_properties(void);

#ifdef HAVE_HDF5
#include <hdf5.h>
void fof_write_hdf5_double_field(double *data, char *tag, int elems, char *description, int blocknr, hid_t grp);
void fof_write_hdf5_int_field(int *data, char *tag, int elems, char *description, hid_t grp);

void fof_write_hdf5_myIDType_field(MyIDType *data, char *tag, int elems, char *description, hid_t grp);
#endif

#ifdef INCL_DMONLY
int list_DM_then_baryons(const void *a, const void *b);
void increase_DM_mass(int size);
#endif

void fof_prepare_id_list(int *ids);
void fof_prepare_star_list(float *stellarage, float *metallicity, float *alphaen, float *initialmass);

extern int Ngroups, TotNgroups;
extern long long TotNids;

extern struct group_properties
{
  int Len;
  int Offset;
  int OffsetType[6];
  MyIDType MinID;
  MyIDType MinIDTask;
  int GrNr;
  int LenType[6];
  double MassType[6];
  double Mass;
  double CM[3];
  double Vel[3];
  MyDouble FirstPos[3];
  double Sfr;
#ifdef BLACK_HOLES
  double BH_Mass;
  double BH_Mdot;
  double MaxDens;
  int index_maxdens, task_maxdens;
#endif

#ifdef BG_SFR /* BG_SFR */
  double InitialMassWeightedStellarAge, InitialMassWeightedStellarBirthZ;
  double MassWeightedTemperature_nsf, MassWeightedTemperature_sf;
  double MassWeightedEntropy_nsf, MassWeightedEntropy_sf;
  double SFMass, NSFMass, StarMass, StarInitialMass;

#ifdef BG_STELLAR_EVOLUTION /* BG_STELLAR_EVOLUTION */
  double SFMetals[BG_NELEMENTS], NSFMetals[BG_NELEMENTS], StarMetals[BG_NELEMENTS];
  double SFMetalMass, NSFMetalMass, StarMetalMass;

#ifdef BG_METALSMOOTHING
  double SFMetalsSmoothed[BG_NELEMENTS], NSFMetalsSmoothed[BG_NELEMENTS], StarMetalsSmoothed[BG_NELEMENTS];
  double SFMetalMassSmoothed, NSFMetalMassSmoothed, StarMetalMassSmoothed;
#endif
#ifdef BG_SNIA_IRON
#ifdef BG_METALSMOOTHING
  double SFIronFromSNIaSmoothed, NSFIronFromSNIaSmoothed, StarIronFromSNIaSmoothed;
#endif
  double SFIronFromSNIa, NSFIronFromSNIa, StarIronFromSNIa;
#endif
#ifdef BG_Z_WEIGHTED_POTENTIAL
  double MetallicityWeightedPotential_sf;
  double MetallicityWeightedPotential_nsf;
  double MetallicityWeightedPotential_stars;
#endif
#ifdef BG_Z_WEIGHTED_REDSHIFT
  double MetallicityWeightedRedshift_sf;
  double MetallicityWeightedRedshift_nsf;
  double MetallicityWeightedRedshift_stars;
#endif
#endif /* BG_STELLAR_EVOLUTION */
#ifdef BG_EXTRA_ARRAYS
  double MassWeightedMaxTemp_sf, MassWeightedMaxTempAExp_sf;
  double MassWeightedMaxEntr_sf, MassWeightedMaxEntrAExp_sf;
  double MassWeightedMaxTemp_nsf, MassWeightedMaxTempAExp_nsf;
  double MassWeightedMaxEntr_nsf, MassWeightedMaxEntrAExp_nsf;
  double MassWeightedMaxTemp_stars, MassWeightedMaxTempAExp_stars;
  double MassWeightedMaxEntr_stars, MassWeightedMaxEntrAExp_stars;
#endif
#ifdef EVALPOTENTIAL
  double MassWeightedPotential_sf;
  double MassWeightedPotential_nsf;
  double MassWeightedPotential_stars;
#endif

#endif /* BG_SFR */
} *Group;

#endif
