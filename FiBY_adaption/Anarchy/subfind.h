#ifndef SUBFIND_H
#define SUBFIND_H





typedef struct
{
  MyIDType index; 
  double r;
  double xproj; // projected radius for yz plane
  double yproj; // projected radius for xz plane
  double zproj; // projected radius for xy plane
  double mass;
  int type;
}
sort_r2list;

void subfind(int num);
#ifdef CONTAMINATION
int subfind_contamination_treefind(MyDouble *searchcenter, MyFloat hsml, int target, int *startnode,
                                      int mode, int *nexport, int *nsend_local, double *Mass);
int subfind_contamination_evaluate(int target, int mode, int *nexport, int *nsend_local);
void subfind_contamination(void);
#endif

int subfind_force_treeevaluate_potential(int target, int mode, int *nexport, int *nsend_local);
void subfind_density(void);
void subfind_overdensity(void);
int subfind_overdensity_evaluate(int target, int mode, int *nexport, int *nsend_local);
double subfind_ovderdens_treefind(MyDouble *searchcenter, MyFloat hsml, int target, int *startnode,
				  int mode, int *nexport, int *nsend_local);
void subfind_save_densities(int num);
void subfind_save_local_densities(int num);
void subfind_setup_smoothinglengths(void);
int subfind_density_evaluate(int target, int mode, int *nexport, int *nsend_local);
void subfind_save_local_catalogue(int num);
void subfind_save_local_catalogue_hdf5(int num);
void subfind_save_final(int num);
int subfind_linkngb_evaluate(int target, int mode, int *nexport, int *nsend_local);
int subfind_ngb_treefind_linkngb(MyDouble *searchcenter, double hsml, int target, int *startnode, int mode,
                                 double *hmax, int *nexport, int *nsend_local);
int subfind_ngb_treefind_nearesttwo(MyDouble *searchcenter, double hsml, int target, int *startnode, int mode,
                                    double *hmax, int *nexport, int *nsend_local);
void subfind_distribute_particles(int mode);
void subfind_unbind_independent_ones(int count);
void subfind_distribute_groups(void);
void subfind_potential_compute(int num, struct unbind_data * d, int phase, double weakly_bound_limit);

void subfind_process_group_collectively(void);
int subfind_col_unbind(struct unbind_data *d, int num);
void subfind_col_determine_sub_halo_properties(struct unbind_data *d, int num, double *mass,
					       double *pos, double *vel, double *cm, double *veldisp,
					       double *vmax, double *vmaxrad, double *spin,
					       MyIDType *mostboundid, double *halfmassrad,
					       double *proj_halfmassrad);
/*
void subfind_col_determine_sub_halo_properties(struct unbind_data *d, int num, double *mass, double *masstype,
					       double *pos, double *vel, double *cm, double *veldisp, double *stellarveldisp,
					       double *vmax, double *vmaxrad, double *spin, MyIDType *mostboundid, double *halfmassrad);
*/
void subfind_col_determine_R200(double hmr, double center[3],
				double *m_Mean200, double *r_Mean200,
				double *m_Crit200, double *r_Crit200, double *m_TopHat, double *r_TopHat);
void subfind_find_linkngb(void);
int subfind_loctree_treebuild(int npart, struct unbind_data *mp);
void subfind_loctree_update_node_recursive(int no, int sib, int father);
double subfind_loctree_treeevaluate_potential(int target);
double subfind_loctree_treeevaluate_accurate_potential(int target);
void subfind_loctree_copyExtent(void);
double subfind_locngb_treefind(MyDouble *xyz, int desngb, double hguess);
void subfind_loctree_findExtent(int npart, struct unbind_data *mp);
int subfind_locngb_treefind_variable(MyDouble *searchcenter, double hguess);
size_t subfind_loctree_treeallocate(int maxnodes, int maxpart);
void subfind_loctree_treefree(void);
void subfind_find_nearesttwo(void);
int subfind_nearesttwo_evaluate(int target, int mode, int *nexport, int *nsend_local);
int subfind_process_group_serial(int gr, int offset, int *iter, int *skipunbind);
int subfind_unbind(struct unbind_data *ud, int len, int *lentype, int *iter, int *skipunbind);
void subfind_determine_sub_halo_properties(struct unbind_data *ud, int num, double *mass, double *masstype,
					   double *pos, double *vel, double *cm, double *veldisp, double *stellarveldisp, double *halfproj_stellarveldisp,
					   double *vmax, double *vmaxrad, double *spin, MyIDType *mostboundid, double *halfmassrad, double *proj_halfmassrad);
int subfind_compare_P_GrNr_DM_Density(const void *a, const void *b);
int subfind_compare_Ind_GrNr_DM_Density(const void *a, const void *b);
int subfind_compare_P_GrNrGrNr(const void *a, const void *b);
int subfind_locngb_compare_key(const void *a, const void *b);
int subfind_compare_serial_candidates_subnr(const void *a, const void *b);
int subfind_compare_serial_candidates_rank(const void *a, const void *b);
int subfind_compare_dens(const void *a, const void *b);
int subfind_compare_energy(const void *a, const void *b);
int subfind_compare_grp_particles(const void *a, const void *b);
int subfind_compare_candidates_boundlength(const void *a, const void *b);
int subfind_compare_candidates_nsubs(const void *a, const void *b);
int subfind_compare_serial_candidates_boundlength(const void *a, const void *b);
int subfind_compare_P_submark(const void *a, const void *b);
int subfind_compare_dist_rotcurve(const void *a, const void *b);
int subfind_compare_xproj_rotcurve(const void *a, const void *b);
int subfind_compare_yproj_rotcurve(const void *a, const void *b);
int subfind_compare_zproj_rotcurve(const void *a, const void *b);
int subfind_compare_binding_energy(const void *a, const void *b);
int subfind_compare_unbind_data_Potential(const void *a, const void *b);
int subfind_compare_unbind_data_seqnr(const void *a, const void *b);
int subfind_compare_densities(const void *a, const void *b);
int subfind_compare_candidates_rank(const void *a, const void *b);
int subfind_ngb_compare_dist(const void *a, const void *b);
int subfind_compare_hsml_data(const void *a, const void *b);
int subfind_compare_ID_list(const void *a, const void *b);
int subfind_compare_SubGroup_GrNr_SubNr(const void *a, const void *b);
int subfind_compare_candidates_subnr(const void *a, const void *b);
void subfind_poll_for_requests(void);
long long subfind_distlinklist_setrank_and_get_next(long long index, long long *rank);
long long subfind_distlinklist_get_rank(long long index);
void subfind_distlinklist_set_next(long long index, long long next);
void subfind_distlinklist_add_particle(long long index);
void subfind_distlinklist_add_bound_particles(long long index, int nsub);
void subfind_distlinklist_mark_particle(long long index, int target, int submark);
long long subfind_distlinklist_get_next(long long index);
long long subfind_distlinklist_get_head(long long index);
void subfind_distlinklist_set_headandnext(long long index, long long head, long long next);
void subfind_distlinklist_set_tailandlen(long long index, long long tail, int len);
void subfind_distlinklist_get_tailandlen(long long index, long long *tail, int *len);
void subfind_distlinklist_set_all(long long index, long long head, long long tail, int len, long long next);
int subfind_distlinklist_get_ngb_count(long long index, long long *ngb_index1, long long *ngb_index2);
long long subfind_distlinklist_set_head_get_next(long long index, long long head);

void initialize_SubGroup(int size);


#ifdef HAVE_HDF5
#include <hdf5.h>
void subfind_save_local_densities_hdf5(int num);

void subfind_write_hdf5_double_field(double *data, char *tag, int elems, char *description, int blocknr, hid_t grp);
void subfind_write_hdf5_int_field(int *data, char *tag, int elems, char *description, hid_t grp);

void subfind_write_hdf5_MyIDType_field(MyIDType *data, char *tag, int elems, char *description, hid_t grp);
void subfind_write_hdf5_MyIDType_particle(MyIDType *data, char *tag, int elems, char *description, hid_t grp);

void subfind_write_hdf5_MyDouble_particle(MyDouble *data, char *tag, int elems, char *description, int blocknr, hid_t grp);
void subfind_write_hdf5_int_particle(int *data, char *tag, int elems, char *description, int blocknr, hid_t grp);
#endif


extern int Ncollective;
extern int MaxNsubgroups;
extern int Nsubgroups;
extern int TotNsubgroups;

/* #ifdef IGNOREFUZZ */
/* extern float sizelimit; */
/* extern int skipped; */
/* #endif */

extern struct subgroup_properties
{
  int Len;
  int LenType[6];
  int GrNr;
  int SubNr;
  int SubParent;    
  int Offset;
  int OffsetType[6];
  MyIDType SubMostBoundID;
  double Mass;
  double MassType[6];
  double SubVelDisp;
  double SubStellarVelDisp;
  double SubStellarVelDispHalfProj;
  double SubVmax;
  double SubVmaxRad;
  double SubHalfMass[6];
  double SubHalfMassProj[6];
  double Pos[3];
  double CM[3];
  double Vel[3];
  double Spin[3];
} *SubGroup;


extern struct nearest_r2_data
{
  double dist[2];
}
*R2Loc;

extern struct nearest_ngb_data
{
  long long index[2];
  int count;
}
*NgbLoc;


extern struct r2data
{
  MyFloat r2;
  int   index;
}
*R2list;

extern double *Dist2list;


extern struct my_serial_index
{
  int Index;
  double Density;
  int GrNr;
}
*Ind;

extern int *Rank;

extern int *SkipGrUnbind;

#endif
