#include "allvars.h"


double bg_time_integ(double a, void *param);

void bg_compute_extra_arrays(void);

void binning(void);

void get_group_name(enum iofields blocknr, char *buf);

double bg_get_temperature(int i);
double bg_get_temperaturep(struct sph_particle_data *sphp, struct particle_data *pp);

double bg_get_entropyp(struct sph_particle_data *sphp);
double bg_get_entropy(int i);

#ifdef BG_STELLAR_EVOLUTION
void output_metallicity_data(void);
#endif

/*
 * ----------------------------------------------------------------------
 * Functions declarations for BlueGene/L interpolation routines
 * ----------------------------------------------------------------------
 */

#if defined(BG_STELLAR_EVOLUTION) || defined(BG_COOLING)

void get_index_1d(float *table, int ntable, double x, int *i, float *dx);

float interpol_1d(float *table, int i, float dx);

double interpol_1d_dbl(double *table, int i, float dx);

float interpol_2d(float **table, int i, int j, float dx, float dy);

double interpol_2d_dbl(double **table, int i, int j, double dx, double dy);

float interpol_3d(float ***table, int i, int j, int k,
		  float dx, float dy, float dz);

float interpol_3d_metals(float ***table, int i, int j, int k, float dy, float dz);

float interpol_4d_old(float ****table, int i, int j, int k, int l,
		      float dx, float dy, float dz, float dw);

float interpol_4d(float ****table, int i, int ip, int j, int jp,
		  int k, int kp, int l, int lp,
		  float dx, float dy, float dz, float dw);

int element_index(char *element_name);

int element_present(char *element_name);

int get_element_index(char **table, int size, char *element_name);

#endif

/*
 * ----------------------------------------------------------------------
 * Functions declarations for BlueGene/L cooling routines
 * ----------------------------------------------------------------------
 */

#if defined(BG_COOLING) || defined(BG_STELLAR_EVOLUTION)
void InitChemistry(void);
#endif

#ifdef BG_COOLING

void InitCool(void);

void LoadCoolingTables(int z_index);

void GetCollisTable(void);
#ifdef BG_COOLING_SHIELDING
void LoadCollisionalCoolingTables();
void GetCollisTableShielding();
#endif

void GetCoolingTables(int highz_table_index, int lowz_table_index);

#ifdef BG_COOLING_OLD_TABLES
void BroadcastCoolingTables(int mode);
#else
void BroadcastCoolingTables(void);
#endif

double CoolingRate(double u, float inn_h, double inZ, float d_z, double particle_Z[], double n_H2, double n_HD);

#ifdef BG_MOL_COOLING
double DoCooling(double u_old, double rho, double dt, double dz, double inz,
		 float d_z, double particle_Z[], double x_H2, double x_HD);
#else
double DoCooling(double u_old, double rho, double dt, double dz, double inz,
		 float d_z, double particle_Z[]);
#endif

void allocate_header_arrays(void);

int set_cooling_SolarAbundances(double *ElementAbundance,
				double *cooling_ElementAbundance);

double Reionization(double z, double dz, float n_h);

void set_solar_metallicity(void);

#endif

void get_redshift_index(float z, int *z_index, float *dz);

double convert_u_to_temp(float d_z, double u, double inn_h, double inHe);

double convert_u_to_temp_old(float d_z, double u, double inn_h, double inHe);


/*
 * ----------------------------------------------------------------------
 * Functions declarations for BlueGene/L stellar evolution routines
 * ----------------------------------------------------------------------
 */

#ifdef BG_SFR
int bg_enrich_evaluate(int target, int mode, int *nexport, int *nsend_local);

void bg_thermal_feedback(void);

void bg_ngbmass(void);

int bg_ngbmass_evaluate(int target, int mode, int *nexport, int *nsend_local);

int bg_thermal_feedback_ngb_treefind(MyDouble searchcenter[3], MyFloat hsml,
				     int target, int *startnode, int mode,
				     int *nexport, int *nsend_local);

void bg_popiii_thermal_feedback(void);

int bg_popiii_thermal_feedback_evaluate(int target, int mode, int *nexport,
					int *nsend_local);

void bg_snii_thermal_feedback(void);

int bg_snii_thermal_feedback_evaluate(int target, int mode, int *nexport,
				      int *nsend_local);

#ifdef BG_DOUBLE_IMF
void bg_snii_kinetic_feedback(int j, unsigned int id, MyDouble pos[3],
			      double mass_star, double ngb_mass,
			      double birth_density);
#else
void bg_snii_kinetic_feedback(int j, unsigned int id, MyDouble pos[3],
			      double mass_star, double ngb_mass);
#endif

double bg_snii_kinetic_feedback_second_loop(void);

int bg_enrich_ngb_treefind(MyDouble searchcenter[3], MyFloat hsml, int target, 
                           int *startnode, int mode, int *nexport, int *nsend_local);

void bg_enrich(void);

int  bg_solidangle_compare_key(const void *a, const void *b);

int  bg_enrich_evaluate_weightsum(int target, int mode, int *nexport, int *nsend_local);

void bg_enrich_evaluate_ngbmass(int target, int mode, int *nexport, int *nsend_local);

void bg_enrich_determine_weights(void);

double bg_evaluate_solid_angle(int target, double r, double h);

void output_wind_data(void);
#endif

#ifdef BG_STELLAR_EVOLUTION

void init_imf(void);
#ifdef BG_POPIII
void init_popiii_imf(void);
#endif
void init_yields(void);

void compute_yields(int);

void compute_ejecta(void);

int read_yield_tables(void);

int bcast_yield_table_dim(void);

int bcast_yield_tables(void);

int allocate_yield_tables(void);

int bg_metals_compare_key(const void *a, const void *b);

double integrate_imf(double log_min_dying_mass, double log_max_dying_mass,
		     double m2, int mode);

void determine_imf_bins(double log_min_dying_mass, double log_max_dying_mass,
			int *ilow, int *ihigh);
#ifdef BG_POPIII
double integrate_popiii_imf(double log_min_dying_mass, double log_max_dying_mass,
			    double m2, int mode);

void determine_popiii_imf_bins(double log_min_dying_mass, double log_max_dying_mass,
			       int *ilow, int *ihigh);
#endif
void get_particle_metal_content(int iPart, double *initial_metals);

void set_particle_metal_content(int iPart, double mass_released);

#ifdef BG_SNIA_IRON
#ifdef BG_DOUBLE_IMF
void evolve_SNIa(double log_min_dying_mass, double log_max_dying_mass,
		 double dt_in_Gyr, double age_of_star_in_Gyr, double metallicity,
		 MyFloat *metals_released, MyFloat *metal_mass_released,
		 MyFloat *iron_from_snia, int mode, double phys_dens);
/* MyFloat *iron_from_snia, double *Number_of_SNIa, int mode, double phys_dens); */
#else
void evolve_SNIa(double log_min_dying_mass, double log_max_dying_mass,
		 double dt_in_Gyr, double age_of_star_in_Gyr, double metallicity,
		 MyFloat *metals_released, MyFloat *metal_mass_released,
		 MyFloat *iron_from_snia, int mode);
/* MyFloat *iron_from_snia, double *Number_of_SNIa, int mode); */
#endif
#else
#ifdef BG_DOUBLE_IMF
void evolve_SNIa(double log_min_dying_mass, double log_max_dying_mass,
		 double dt_in_Gyr, double age_of_star_in_Gyr, double metallicity,
		 MyFloat *metals_released, MyFloat *metal_mass_released,
		 int mode, double phys_dens);
/* double *Number_of_SNIa, int mode, double phys_dens); */
#else
void evolve_SNIa(double log_min_dying_mass, double log_max_dying_mass,
		 double dt_in_Gyr, double age_of_star_in_Gyr, double metallicity,
		 MyFloat *metals_released, MyFloat *metal_mass_released,
		 int mode);
/* double *Number_of_SNIa, int mode); */
#endif
#endif

#ifdef BG_DOUBLE_IMF
void evolve_SNII(double log_min_mass, double log_max_mass, double log_metallicity,
		 double *initial_metals, MyFloat *metals_released,
		 MyFloat *metal_mass_released, MyFloat *dust_mass_released, double phys_dens);
/* MyFloat *metal_mass_released, double *Number_of_SNII, double phys_dens); */
#else
void evolve_SNII(double log_min_mass, double log_max_mass, double log_metallicity,
		 double *initial_metals, MyFloat *metals_released,
		 MyFloat *metal_mass_released, MyFloat *dust_mass_released);
/* MyFloat *metal_mass_released, double *Number_of_SNII); */
#endif

#ifdef BG_DOUBLE_IMF
void evolve_AGB(double log_min_mass, double log_max_mass, double log_metallicity,
		double *initial_metals, MyFloat *metals_released,
		MyFloat *metal_mass_released, MyFloat *dust_mass_released, double phys_dens);
#else
void evolve_AGB(double log_min_mass, double log_max_mass, double log_metallicity,
		double *initial_metals, MyFloat *metals_released,
		MyFloat *metal_mass_released, MyFloat *dust_mass_released);
#endif

#ifdef BG_POPIII
void evolve_POPIII(double log_min_mass, double log_max_mass, double log_metallicity,
		   double *initial_metals, MyFloat *metals_released,
		   MyFloat *metal_mass_released, MyFloat *dust_mass_released);

void popiii_black_holes();
#endif

#ifdef BG_SNIA_IRON
void bg_stellar_evolution(double age_of_star_in_Gyr, double dtime_in_Gyr,
			  int iPart, MyFloat *metals_released, MyFloat *metal_mass_released,
			  MyFloat *dust_mass_released, MyFloat *iron_from_snia, MyFloat *energy_released);
/* double *Number_of_SNIa); */
#else
void bg_stellar_evolution(double age_of_star_in_Gyr, double dtime_in_Gyr,
			  int iPart, MyFloat *metals_released, MyFloat *metal_mass_released,
			  MyFloat *dust_mass_released, MyFloat *energy_released);
/* MyFloat *energy_released, double *Number_of_SNIa); */
#endif

double dying_mass_msun(double time_in_Gyr, double metallicity);

double lifetime_in_Gyr(double mass, double metallicity);

#ifdef BG_SNIA_IRON
#ifdef BG_DOUBLE_IMF
void stellar_evolution(double age_of_star_in_Gyr, double dtime_in_Gyr, double metallicity,
		       double *initial_metals, MyFloat *metals_released,
		       MyFloat *metal_mass_released, MyFloat *dust_mass_released,
		       MyFloat *iron_from_snia, double phys_dens);
/* double *Number_of_SNIa, double *Number_of_SNII, double phys_dens); */
#else
void stellar_evolution(double age_of_star_in_Gyr, double dtime_in_Gyr, double metallicity,
		       double *initial_metals, MyFloat *metals_released,
		       MyFloat *metal_mass_released, MyFloat *dust_mass_released,
		       MyFloat *iron_from_snia);
/* double *Number_of_SNIa, double *Number_of_SNII); */
#endif
#else
#ifdef BG_DOUBLE_IMF
void stellar_evolution(double age_of_star_in_Gyr, double dtime_in_Gyr, double metallicity,
		       double *initial_metals, MyFloat *metals_released,
		       MyFloat *metal_mass_released, MyFloat *dust_mass_released,
		       double phys_dens);
/* MyFloat *metal_mass_released, double *Number_of_SNIa, double *Number_of_SNII, double phys_dens); */
#else
void stellar_evolution(double age_of_star_in_Gyr, double dtime_in_Gyr, double metallicity,
		       double *initial_metals, MyFloat *metals_released,
		       MyFloat *metal_mass_released, MyFloat *dust_mass_released);
/* MyFloat *metal_mass_released, double *Number_of_SNIa, double *Number_of_SNII); */
#endif
#endif

void test_yields(void);

void output_snia_data(void);

#endif

#ifdef BLACK_HOLES
void bg_bh_thermal_feedback(void);

int bg_bh_thermal_feedback_ngb_treefind(MyDouble searchcenter[3],
					  MyFloat hsml, int target,
					  int *startnode, int mode,
					  int *nexport, int *nsend_local);

int bg_bh_thermal_feedback_evaluate(int target, int mode, int *nexport,
				      int *nsend_local);

#endif
