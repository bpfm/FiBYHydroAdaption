#define N_MASS_BIN 200 /* IMF number of bins */
#define POPIII_N_MASS_BIN 200 /* IMF number of bins for popIII */

#define MIN_METAL -20  /* minimum metallicity */

extern double NumberOfSNIa, NumberOfSNII;

extern double *imf_by_number, *imf_mass_bin, *imf_mass_bin_log10;
#ifdef BG_DOUBLE_IMF
extern double *imf_by_number1;
#endif
extern double *yield_mass_bin, *stellar_yield, *integrand;

#ifdef BG_POPIII
extern double *popiii_imf_by_number, *popiii_imf_mass_bin, *popiii_imf_mass_bin_log10;
extern double *popiii_yield_mass_bin, *popiii_stellar_yield, *popiii_integrand;
#endif

extern struct YieldTypeIa
{
  int N_ELEMENTS;               /* number of yields */
  char **ElementName;           /* name of species */
  double *Yield;                /* Yield (in solar masses) */
  double *SPH;                  /* Yields for corresponding SPH elements (in solar masses) */
  double TotalMetals_SPH;       /* Total metal mass */
} yieldsSNIa;

extern struct YieldTypeII_and_AGB
{
  int N_ELEMENTS, N_MASS, N_Z;  /* number of elements, mass, and initial metallicity bins */
  char **ElementName;           /* name of element */
  double *Mass;                 /* mass bins[N_MASS] */
  double *Metallicity;          /* metallicity bins[N_Z] (total metallicity/total mass) */
  double **Ejecta;              /* size [N_MASS*N_Z] with index_3d(imass,iz)  CHANGE MADE HERE! */
  double **Ejecta_SPH;          /* size [N_MASS*N_Z] with index_3d(imass,iz)  CHANGE MADE HERE! */
  double ***Yield;              /* size [N_ELEMENTS*N_MASS*N_Z] with index_3d(iel,imass,iz) */
  double ***SPH;                /* Yields for corresponding SPH elements (in solar masses) */
  double **TotalMetals;         /* Total metal mass */
  double **TotalMetals_SPH;     /* Total metal mass */
} yieldsSNII, yieldsAGB;

#ifdef BG_DUST
extern struct DustTypeII_and_AGB
{
  int N_MASS, N_Z;              /* number of mass and initial metallicity bins */
  double *Mass;                 /* mass bins[N_MASS] */
  double *Metallicity;          /* metallicity bins[N_Z] (total metallicity/total mass) */
  double **TotalDust;           /* Total dust mass */
  double **TotalDust_SPH;       /* Total dust mass */
} dustSNII, dustAGB;
#endif

extern struct Lifetime_Table
{
  int N_MASS, N_Z;              /* number of elements, mass, and initial metallicity bins */
  double *Mass;                 /* mass bins[N_MASS]     */
  double *Metallicity;          /* metallicity bins[N_Z] (total metallicity/total mass) */
  double **Dyingtime;           /* size [N_MASS*N_Z] with index_2d(imass,iz) */
} Lifetimes;

#ifdef BG_POPIII
extern struct YieldPOPIII
{
  int N_ELEMENTS, N_MASS;       /* number of elements and mass bins */
  char **ElementName;           /* name of element */
  double *Mass;                 /* mass bins[N_MASS] */
  double *Ejecta;               /* size [N_MASS] */
  double *Ejecta_SPH;           /* size [N_MASS] */
  double **Yield;               /* size [N_ELEMENTS*N_MASS] */
  double **SPH;                 /* Yields for corresponding SPH elements (in solar masses) */
  double *TotalMetals;          /* Total metal mass */
  double *TotalMetals_SPH;      /* Total metal mass */
} yieldsPOPIII;

extern struct POPIIILifetime_Table
{
  int N_MASS;                   /* number of elements, mass, and initial metallicity bins */
  double *Mass;                 /* mass bins[N_MASS] */
  double *Dyingtime;            /* size [N_MASS] */
} POPIIILifetimes;

#ifdef BG_DUST
extern struct DustPOPIII
{
  int N_MASS;                   /* number of mass bins */
  double *Mass;                 /* mass bins[N_MASS] */
  double *TotalDust;            /* Total dust mass */
  double *TotalDust_SPH;        /* Total dust mass */
} dustPOPIII;
#endif
#endif


/* names of yield tables */
#define AGB_yieldname  "AGB.hdf5"
#define SNIa_yieldname "SNIa.hdf5"
#define SNII_yieldname "SNII.hdf5"
#define Lifetime_name  "Lifetimes.hdf5"

#ifdef BG_DUST
#define AGB_dustname    "AGB_dust.hdf5"
#define SNII_dustname   "SNII_dust.hdf5"
#endif

#ifdef BG_POPIII
#define POPIII_yieldname "POPIII.hdf5"
#define POPIIILifetime_name  "POPIII_Lifetimes.hdf5"
#ifdef BG_DUST
#define POPIII_dustname "POPIII_dust.hdf5"
#endif
#endif

/* constants for SNIa */
#define A_SNIa_Lia   0.05
#define gamma_SNIa   2.
#define A_SNIa_scan  4.4e-2 /* not used */
#define B_SNIa_scan  2.6 /* not used */

#define TYPEIA_MIN_MASS 3.0 /* caution! Ia rates have been calibrated with 3 solar mass! */
#define TYPEIA_MAX_MASS 8.0 /* caution! Ia rates have been calibrated with 8 solar mass! */
#define TAU_SNIa_MAN_G  3.3 /* (tau for the gaussian SNIa rate) */
#define SIGMA_SNIa_G    0.66 /* (sigma for the gaussian SNIa rate) */

#define TAU_SNIa_MAN_E 2.0 /* (tau for the efolding SNIa rate) */

#define TAU_SNIa_MAN_EC 3.0 /* (tau for the efolding part of the combination rate) */
#define TAU_SNIa_MAN_GC 0.05 /* (tau for the gaussian part combination rate) */
#define SIGMA_SNIa_GC   0.01 /* (sigma for the gaussian part combination rate) */
