/* conversion from eV to K */
#define EV_IN_K 11604.50520

/* the foolowing coul be put in as run parameters */
#define EPS 0.1 /* set the chemical timestep as EPS * n / (dn/dt) */
#define MIN_CHEM_TSTEP_YR 1000
#define MOL_NETWORK_TEMP_THRESH 1.3e4

/* definition of rate coefficient tables */
#ifdef BG_MOL_NETWORK_TABULATED_RATES
#define MOL_NETWORK_TABLE_SIZE 500

#define MOL_NETWORK_MIN_TEMP 1.0
#define MOL_NETWORK_MAX_TEMP 3e4

#define MOL_NETWORK_MIN_TEMP_LOG10 (log10(MOL_NETWORK_MIN_TEMP))
#define MOL_NETWORK_MAX_TEMP_LOG10 (log10(MOL_NETWORK_MAX_TEMP))

#define MOL_NETWORK_DLOGTEMP ((MOL_NETWORK_MAX_TEMP_LOG10-MOL_NETWORK_MIN_TEMP_LOG10)/(MOL_NETWORK_TABLE_SIZE-1))
#define MOL_NETWORK_DLOGTEMP_M1 (1.0/MOL_NETWORK_DLOGTEMP)
#endif

/* minimum allowed fraction */
#define MOL_NETWORK_SMALL_NUMBER 1e-30

/* tables */
#ifdef BG_MOL_NETWORK_TABULATED_RATES
extern double *ttab;

extern double *k1tab;
extern double *k2tab;
extern double *k3tab;
extern double *k4tab;
extern double *k5tab;
extern double *k6tab;
extern double *k7tab;
extern double *k8tab;
extern double *k9tab;

extern double *k10tab;
extern double *k11tab;
extern double *k12tab;
extern double *k13tab;
extern double *k14tab;
extern double *k15tab;
extern double *k16tab;
extern double *k17tab;
extern double *k18tab;
extern double *k19tab;

extern double *k20tab;
extern double *k21tab;
extern double *k22tab;
extern double *k23tab;
extern double *k24tab;

extern double *d1tab;
extern double *d2tab;
extern double *d3tab;
extern double *d4tab;
extern double *d5tab;
extern double *d6tab;
#endif

/* functions */
void bg_molecules();

void bg_molecules_evaluate(int i);

double prod_H(double n_H, double n_Hp, double n_Hm, double n_H2, double n_H2p, double n_He, double n_D, double n_e);
double dest_H(double n_H, double n_Hp, double n_Hm, double n_H2, double n_H2p, double n_Hep, double n_Dp, double n_HD, double n_e);

double prod_Hp(double n_H, double n_H2, double n_H2p, double n_Hep, double n_Dp, double n_e);
double dest_Hp(double n_H, double n_Hm, double n_H2, double n_He, double n_D, double n_HD, double n_e);

double prod_Hm(double n_H, double n_e);
double dest_Hm(double n_H, double n_Hp, double n_H2p, double n_e);

double prod_H2(double n_H, double n_Hp, double n_Hm, double n_H2, double n_H2p, double n_HD);
double dest_H2(double n_H, double n_Hp, double n_H2, double n_D, double n_Dp, double n_e);

double prod_H2p(double n_H, double n_Hp, double n_Hm, double n_H2);
double dest_H2p(double n_H, double n_Hm, double n_e);

double prod_He(double n_H, double n_Hep, double n_e);
double dest_He(double n_Hp, double n_e);

double prod_Hep(double n_Hp, double n_He, double n_Hepp, double n_e);
double dest_Hep(double n_H, double n_e);

double prod_Hepp(double n_Hep, double n_e);
double dest_Hepp(double n_e);

double prod_D(double n_H, double n_Dp, double n_HD);
double dest_D(double n_Hp, double n_H2);

double prod_Dp(double n_Hp, double n_D, double n_HD);
double dest_Dp(double n_H, double n_H2);

double prod_HD(double n_H2, double n_D, double n_Dp);
double dest_HD(double n_H, double n_Hp);

double prod_e(double n_H, double n_Hp, double n_Hm, double n_H2, double n_He, double n_Hep, double n_e);
double dest_e(double n_H, double n_Hp, double n_Hm, double n_H2, double n_H2p, double n_He, double n_Hep, double n_Hepp, double n_e);

double free_electrons(double n_Hp, double n_Hm, double n_H2p, double n_Hep, double n_Hepp, double n_Dp);

double k1(double temp);
double k2(double temp);
double k3(double temp);
double k4(double temp);
double k5(double temp);
double k6(double temp);
double k7(double temp);
double k8(double temp);
double k9(double temp);
double k10(void);
double k11(double temp);
double k12(double temp);
double k13(double temp);
double k14(double temp);
double k15(double temp);
double k16(double temp);
double k17(double temp);
double k18(double temp);
double k19(double temp);
double k20(double temp);
double k21(double temp);
double k22(double temp);
double k23(double temp);
double k24(double temp);

double d1(double temp);
double d2(double temp);
double d3(double temp);
double d4(double temp);
double d5(double temp);
double d6(double temp);

#ifdef BG_MOL_NETWORK_TABULATED_RATES
double k_table_interpol(double *table, double temp);
void init_rate_coefficient_tables(void);
void output_rate_tables(void);
#endif

void set_rate_coefficients(double temp);

void output_rates(void);


