/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define statements
 *     - add #include "allvars.h" 
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#include "tags.h"


#define  GADGETVERSION   "3.0-BG"      /*!< code version string */
#define  ANARCHYVERSION  "0.0-Cosmo"   /*!< code version string */


#define  TIMEBINS        29

#define  TIMEBASE        (1<<TIMEBINS) /*!< The simulated timespan is mapped onto the integer interval [0,TIMESPAN],
                                        *   where TIMESPAN needs to be a power of 2. Note that (1<<28) corresponds
                                        *   to 2^29
                                        */
#define  TISTEP(x) ((x) ? (1 << (x)) : 0)

/* #ifdef TIMESTEP_LIMITER */
/* #define DT_LIMITER_FACTOR 4 */
/* #define DT_LIMITER_FACTOR_LOG2 2 */

/* #define BITFLAG_DT_FEEDBACK   0 */
/* #define BITFLAG_DT_LIMITER    1 */
/* #endif */

#ifndef  MULTIPLEDOMAINS
#define  MULTIPLEDOMAINS     1
#endif

#ifndef  TOPNODEFACTOR
#define  TOPNODEFACTOR       2.5
#endif

#define  NODELISTLENGTH      8

typedef unsigned long long peanokey;

#define  BITS_PER_DIMENSION 21	/* for Peano-Hilbert order. Note: Maximum is 10 to fit in 32-bit integer ! */
#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))


#define  PI            3.1415927  /*!< this is pi */

#ifdef USER_GAMMA
#define  GAMMA            USER_GAMMA
#else
#define  GAMMA            (5.0/3)	/*!< adiabatic index of simulated gas */
#endif

#define  GAMMA_INV        (1.0/GAMMA)
#define  GAMMA_MINUS1     (GAMMA-1)
#define  GAMMA_MINUS1_INV (1.0/GAMMA_MINUS1)


#define  HYDROGEN_MASSFRAC 0.76	/*!< mass fraction of hydrogen, relevant only for radiative cooling */


#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#define  RNDTABLE 8192

/* ... often used physical constants (cgs units) */

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  CM_PER_KM   1e5
#define  PROTONMASS  1.6726e-24
#define  ELECTRONMASS 9.10953e-28
#define  THOMPSON     6.6524587e-25
#define  ELECTRONCHARGE  4.8032e-10
#define  HUBBLE          3.2407789e-18	/* in h/sec */
#define  LYMAN_ALPHA      1215.6e-8	/* 1215.6 Angstroem */
#define  LYMAN_ALPHA_HeII  303.8e-8	/* 303.8 Angstroem */
#define  OSCILLATOR_STRENGTH       0.41615
#define  OSCILLATOR_STRENGTH_HeII  0.41615
#define  T_CMB0      2.728	/* present-day CMB temperature */
#define  ZSOLAR      0.012663729 /* fraction of solar mass in metals */
#ifndef M_PI
#define M_PI 3.1415926535897932385
#endif

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7


#ifndef ASMTH
/*! ASMTH gives the scale of the short-range/long-range force split in units of FFT-mesh cells */
#define ASMTH 1.25
#endif

#ifndef RCUT
/*! RCUT gives the maximum distance (in units of the scale used for the force split) out to which short-range
 * forces are evaluated in the short-range tree walk.
 */
#define RCUT  4.5		
#endif

#ifndef HPM_SMTH
#define HPM_SMTH ASMTH
#endif


#define MAXLEN_OUTPUTLIST 350	/*!< maxmimum number of entries in output list */

#define DRIFT_TABLE_LENGTH  1000	/*!< length of the lookup table used to hold the drift and kick factors */


#define MAXITER 150

#define LINKLENGTH 0.2
#define GROUP_MIN_LEN 25

#define MINRESTFAC 0.05

/* size of description string in io.c */
#define NDESC 200


#ifndef LONGIDS
typedef unsigned int MyIDType;
#else
typedef unsigned long long MyIDType;
#endif

#ifndef DOUBLEPRECISION     /* default is single-precision */
typedef float  MyFloat;
typedef float  MyDouble;
#endif

#if (DOUBLEPRECISION == 1)  /* everything double-precision */
typedef double  MyFloat;
typedef double  MyDouble;
#endif

#if (DOUBLEPRECISION == 2)   /* mixed precision */
typedef float   MyFloat;
typedef double  MyDouble;
#endif

#ifdef FLTROUNDOFFREDUCTION  
#define FLT(x) ((MyFloat)(x))
#ifdef SOFTDOUBLEDOUBLE      /* this requires a C++ compilation */
#include "dd.h"
typedef dd MyLongDouble;
#else
typedef long double MyLongDouble;
#endif
#else  /* not enabled */
#define FLT(x) (x)
typedef MyFloat MyLongDouble;
#endif  /* end FLTROUNDOFFREDUCTION */

struct unbind_data
{
  int index;
};


#define CPU_ALL            0
#define CPU_TREEWALK1      1
#define CPU_TREEWALK2      2
#define CPU_TREEWAIT1      3
#define CPU_TREEWAIT2      4
#define CPU_TREESEND       5
#define CPU_TREERECV       6
#define CPU_TREEMISC       7
#define CPU_TREEBUILD      8
#define CPU_TREEUPDATE     9
#define CPU_TREEHMAXUPDATE 10
#define CPU_DOMAIN         11
#define CPU_DENSCOMPUTE    12
#define CPU_DENSWAIT       13
#define CPU_DENSCOMM       14
#define CPU_DENSMISC       15
#define CPU_HYDCOMPUTE     16
#define CPU_HYDWAIT        17
#define CPU_HYDCOMM        18
#define CPU_HYDMISC        19
#define CPU_DRIFT          20
#define CPU_TIMELINE       21
#define CPU_POTENTIAL      22
#define CPU_MESH           23
#define CPU_PEANO          24
#define CPU_COOLINGSFR     25        
#define CPU_SNAPSHOT       26
#define CPU_FOF            27
#define CPU_BLACKHOLES     28
#define CPU_MISC           29
#define CPU_LINEOFSIGHT    30
#define CPU_BINNING        31
#define CPU_ENRICH         32
#define CPU_ENRICHSTEVOL   33
#define CPU_COOLING        34
#define CPU_SFR            35
#define CPU_MOLECULES      36
#define CPU_PARTS          37  /* this gives the number of parts above (must be last) */

#define CPU_STRING_LEN 120

extern double CPU_Step[CPU_PARTS];
extern char CPU_Symbol[CPU_PARTS];
extern char CPU_SymbolImbalance[CPU_PARTS];
extern char CPU_String[CPU_STRING_LEN+1];


#ifndef  TWODIMS
#define  NUMDIMS 3                                      /*!< For 3D-normalized kernel */
/* #define  KERNEL_COEFF_1  2.546479089470                 /\*!< Coefficients for SPH spline kernel and its derivative *\/  */
/* #define  KERNEL_COEFF_2  15.278874536822 */
/* #define  KERNEL_COEFF_3  45.836623610466 */
/* #define  KERNEL_COEFF_4  30.557749073644 */
/* #define  KERNEL_COEFF_5  5.092958178941 */
/* #define  KERNEL_COEFF_6  (-15.278874536822) */
/* #define  NORM_COEFF      4.188790204786                 /\*!< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 *\/  */
#else
#define  NUMDIMS 2                                      /*!< For 2D-normalized kernel */
/* #define  KERNEL_COEFF_1  (5.0/7*2.546479089470)         /\*!< Coefficients for SPH spline kernel and its derivative *\/  */
/* #define  KERNEL_COEFF_2  (5.0/7*15.278874536822) */
/* #define  KERNEL_COEFF_3  (5.0/7*45.836623610466) */
/* #define  KERNEL_COEFF_4  (5.0/7*30.557749073644) */
/* #define  KERNEL_COEFF_5  (5.0/7*5.092958178941) */
/* #define  KERNEL_COEFF_6  (5.0/7*(-15.278874536822)) */
/* #define  NORM_COEFF      M_PI                           /\*!< Coefficient for kernel normalization. *\/ */
#endif

#define NUMDIMS_INV (1.0/NUMDIMS)


#if defined(BG_SFR) || defined(BG_COOLING) || defined(BG_STELLAR_EVOLUTION)

#define STR_LENGTH     100
#define BG_NELEMENTS     9

#ifdef BG_COOLING_OLD_TABLES
#define EL_NAME_LENGTH  20
#define TABLE_INDEX_HHE_COOLING 0
#define TABLE_INDEX_Z_COOLING   1
#define TABLE_INDEX_THERMTOTEMP 2
#else
#define EL_NAME_LENGTH  10
#endif

#endif


#if defined (BLACK_HOLES) || defined (BG_SFR)
#define PPP P
#else
#define PPP SphP
#endif


#define DMAX(a,b) (dmax1=(a),dmax2=(b),(dmax1>dmax2)?dmax1:dmax2)
#define DMIN(a,b) (dmin1=(a),dmin2=(b),(dmin1<dmin2)?dmin1:dmin2)
#define IMAX(a,b) (imax1=(a),imax2=(b),(imax1>imax2)?imax1:imax2)
#define IMIN(a,b) (imin1=(a),imin2=(b),(imin1<imin2)?imin1:imin2)

#ifdef PERIODIC
extern MyDouble boxSize, boxHalf;
#ifdef LONG_X
extern MyDouble boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
extern MyDouble boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
extern MyDouble boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif

extern size_t AllocatedBytes;
extern size_t FreeBytes;

#ifdef PERIODIC
#define NGB_PERIODIC_LONG_X(x) (xtmp=fabs(x),(xtmp>boxHalf_X)?(boxSize_X-xtmp):xtmp)
#define NGB_PERIODIC_LONG_Y(x) (xtmp=fabs(x),(xtmp>boxHalf_Y)?(boxSize_Y-xtmp):xtmp)
#define NGB_PERIODIC_LONG_Z(x) (xtmp=fabs(x),(xtmp>boxHalf_Z)?(boxSize_Z-xtmp):xtmp)
#else
#define NGB_PERIODIC_LONG_X(x) fabs(x)
#define NGB_PERIODIC_LONG_Y(x) fabs(x)
#define NGB_PERIODIC_LONG_Z(x) fabs(x)
#endif

#define FACT1 0.366025403785	/* FACT1 = 0.5 * (sqrt(3)-1) */



/*********************************************************/
/*  Global variables                                     */
/*********************************************************/

/* #ifdef FORCE_DOMAIN_DECOMPOSITION */
/* extern double PartsPerSecond; */
/* #endif */

extern int FirstActiveParticle;
extern int *NextActiveParticle;

extern int TimeBinCount[TIMEBINS];
extern int TimeBinCountSph[TIMEBINS];
extern int TimeBinActive[TIMEBINS];

extern int FirstInTimeBin[TIMEBINS];
extern int LastInTimeBin[TIMEBINS];
extern int *NextInTimeBin;
extern int *PrevInTimeBin;

#ifdef BG_SFR
extern double TimeBinSfr[TIMEBINS];
#ifdef BG_POPIII
extern double TimeBinPOPIIISfr[TIMEBINS];
#endif
#endif

#ifdef BLACK_HOLES
extern double TimeBin_BH_mass[TIMEBINS];
extern double TimeBin_BH_dynamicalmass[TIMEBINS];
extern double TimeBin_BH_Mdot[TIMEBINS];
extern double TimeBin_BH_Medd[TIMEBINS];
#endif

extern int ThisTask;		/*!< the number of the local processor  */
extern int NTask;		/*!< number of processors */
extern int PTask;		/*!< note: NTask = 2^PTask */

extern double CPUThisRun;	/*!< Sums CPU time of current process */

extern int NumForceUpdate;	/*!< number of active particles on local processor in current timestep  */
extern long long GlobNumForceUpdate;

extern int NumSphUpdate;	/*!< number of active SPH particles on local processor in current timestep  */

extern int MaxTopNodes;	        /*!< Maximum number of nodes in the top-level tree used for domain decomposition */

extern int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from
				   initial conditions, 1 is resuming a run from a set of restart files, while 2
				   marks a restart from a snapshot file. */

extern int RestartSnapNum;

extern int *Exportflag;	        /*!< Buffer used for flagging whether a particle needs to be exported to another process */
extern int *Exportnodecount;
extern int *Exportindex;

extern int *Send_offset, *Send_count, *Recv_count, *Recv_offset, *Sendcount_matrix;


extern int Flag_FullStep;       /*!< Flag used to signal that the current step involves all particles */


extern int TreeReconstructFlag;  
                
extern double WallclockTime;    /*!< This holds the last wallclock time measurement for timings measurements */


extern int GlobFlag;

extern int NumPart;		/*!< number of particles on the LOCAL processor */
extern int N_gas;		/*!< number of gas particles on the LOCAL processor  */
extern long long Ntype[6];      /*!< total number of particles of each type */
extern int NtypeLocal[6];       /*!< local number of particles of each type */
#if defined(BG_SFR)
extern int N_star;             /*!< number of star particles on the LOCAL processor  */
#endif

extern gsl_rng *random_generator;	/*!< the random number generator used */


#ifdef BG_SFR
extern int Stars_converted;	/*!< current number of star particles in gas particle block */
#endif


extern double TimeOfLastTreeConstruction;	/*!< holds what it says */

extern int *Ngblist;		/*!< Buffer to hold indices of neighbours retrieved by the neighbour search
				   routines */


extern double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
extern int *DomainStartList, *DomainEndList;



extern double *DomainWork;
extern int *DomainCount;
extern int *DomainCountSph;
extern int *DomainTask;
extern int *DomainNodeIndex;
extern int *DomainList, DomainNumChanged;

#if defined(BG_SFR)
extern int *DomainCountStar;
#endif

#ifdef SUBFIND
extern int GrNr;
extern int NumPartGroup;
#endif


extern peanokey *Key, *KeySorted;

extern struct topnode_data
{
  peanokey Size;
  peanokey StartKey;
  long long Count;
  MyFloat GravCost;
  int Daughter;
  int Pstart;
  int Blocks;
  int Leaf;
} *TopNodes;

extern int NTopnodes, NTopleaves;




extern double RndTable[RNDTABLE];


/* variables for input/output , usually only used on process 0 */


extern char ParameterFile[100];  /*!< file name of parameterfile used for starting the simulation */

extern FILE *FdInfo,   /*!< file handle for info.txt log-file. */
  *FdEnergy,           /*!< file handle for energy.txt log-file. */
  *FdTimings,          /*!< file handle for timings.txt log-file. */
  *FdBalance,          /*!< file handle for balance.txt log-file. */
  *FdCPU;              /*!< file handle for cpu.txt log-file. */

#ifdef BG_SFR
extern FILE *FdSfr;    /*!< file handle for sfr.txt log-file. */
#ifdef BG_POPIII
extern FILE *FdPOPIIISfr;    /*!< file handle for popiii_sfr.txt log-file. */
#endif
#endif

#ifdef BG_STELLAR_EVOLUTION
extern FILE *FdMetGas;    /*!< file handle for metals_gas.txt log-file. */
extern FILE *FdMetStars;  /*!< file handle for metals_star.txt log-file. */
extern FILE *FdMetSF;     /*!< file handle for metals_sf.txt log-file. */
extern FILE *FdMetTot;    /*!< file handle for metals_tot.txt log-file. */
extern FILE *FdSNIa;      /*!< file handle for SNIa.txt log-file. */
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK
extern FILE *FdKineticFeedback;      /*!< file handle for wind.txt log-file. */
#endif

#ifdef BG_SNII_THERMAL_FEEDBACK
extern FILE *FdThermalFeedback;      /*!< file handle for thermal_feedback.txt log-file. */
#endif

#ifdef BLACK_HOLES
extern FILE *FdBlackHoles;    /*!< file handle for blackholes.txt log-file. */
#ifdef VERBOSE_LOGFILES
extern FILE *FdBlackHolesDetails;
extern FILE *FdBlackHolesFeedback;
#endif

#ifdef FOF
extern FILE *FdBlackHolesSeeds;
#endif

#ifdef BH_DEBUG
extern FILE *FdBlackHolesEnergy;
#endif
#endif



#ifdef BG_SNII_KINETIC_FEEDBACK
extern MyFloat ExpectedMassLoading;
extern MyFloat ActualMassLoading;
extern MyFloat ActualMassLoading2ndLoop;
#endif



/*! table for the cosmological drift factors */
extern double DriftTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for gravitational forces */
extern double GravKickTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for hydrodynmical forces */
extern double HydroKickTable[DRIFT_TABLE_LENGTH];

extern void *CommBuffer;	/*!< points to communication buffer, which is used at a few places */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
extern struct global_data_all_processes
{
  long long TotNumPart;		/*!<  total particle numbers (global value) */
  long long TotN_gas;		/*!<  total gas particle number (global value) */

#ifdef BLACK_HOLES
  int TotBHs;
#endif

  int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one
				   processor. */
  int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one
				     processor. */

#ifdef BG_SFR
  int MaxPartStar;
#endif

  int ICFormat;			/*!< selects different versions of IC file-format */

  int SnapFormat;		/*!< selects different versions of snapshot file-formats */

  int NumFilesPerSnapshot;      /*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel; /*!< maximum number of files that may be written simultaneously when
                                   writing/reading restart-files, or when writing snapshot files */ 

  int BufferSize;		/*!< size of communication buffer in MB */
  int BunchSize;		/*!< number of particles fitting into the buffer in the parallel tree algorithm  */

  int DoDynamicUpdate;

#ifdef BG_SFR
  int Generations;
  double SF_EOSGammaEffective;
  double SF_EOSMinPhysDens_HpCM3;
  double EOSPhysDensThresh;
  double SF_EOSMinOverDens;
  double EOSOverDensThresh;
  double SF_EOSEnergyAtThreshold_ERG;
  double SF_EOSTempThreshMargin_DEX;
  double SF_SchmidtLawCoeff_MSUNpYRpKPC2;
  double SF_SchmidtLawExponent;
  double SF_SchmidtLawCoeff_GpSpCM2;
#ifdef BG_DOUBLE_IMF
  double SF_SchmidtLawCoeff1_MSUNpYRpKPC2;
  double SF_SchmidtLawCoeff1_GpSpCM2;
#endif
#endif

#ifdef BG_SFR
#ifdef BG_SNII_KINETIC_FEEDBACK
  int    SNII_WindOn;
  int    SNII_WindIsotropicOn;
  double SNII_WindSpeed_KMpS;
  double SNII_WindMassLoading;

#ifdef BG_DOUBLE_IMF
  double SNII_WindMassLoading1;
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING
double SNII_WindDecouplingTime_YR;
#endif
#endif /* BG_SNII_KINETIC_FEEDBACK */

#if defined(BG_SNII_KINETIC_FEEDBACK) || defined(BG_SNII_THERMAL_FEEDBACK)
  double SNII_WindDelay_YR;
#endif
#endif /* BG_SFR */

#if defined(BG_COOLING) || defined(BG_STELLAR_EVOLUTION)
  double InitAbundance_Hydrogen;
  double InitAbundance_Helium;
  double InitAbundance_Carbon;
  double InitAbundance_Nitrogen;
  double InitAbundance_Oxygen;
  double InitAbundance_Neon;
  double InitAbundance_Magnesium;
  double InitAbundance_Silicon;
  double InitAbundance_Iron;

  double CalciumOverSilicon;
  double SulphurOverSilicon;
#endif

#ifdef BG_COOLING
  int    MetDepCoolingOn;
  double REION_H_ZSigma;
  double REION_H_ZCenter;
  double REION_H_Heating_EVpH;
  double REION_H_Heating_ERGpG;
  double REION_He_ZSigma;
  double REION_He_ZCenter;
  double REION_He_Heating_EVpH;
  double REION_He_Heating_ERGpG;
  char   CoolTablePath[STR_LENGTH];

#ifdef BG_COOLING_SHIELDING
  double UVBG_PhysDensThresh_HpCM3;
  double UVBG_PhysDensThresh;
#endif
#endif

#ifdef BG_STELLAR_EVOLUTION
  /* IMF definition */
  char   IMF_Model[STR_LENGTH];
  char   IMF_LifetimeModel[STR_LENGTH];
  int    IMF_LifetimeMode;
#ifdef BG_DOUBLE_IMF
  double IMF_Exponent1;
  double IMF_PhysDensThresh_HpCM3;
  double IMF_PhysDensThresh;
#endif
  double IMF_Exponent;
  double IMF_MinMass_MSUN;
  double IMF_MaxMass_MSUN;
  /* SNIa parameters and flags */
  char SNIa_Model[STR_LENGTH];
  double SNIa_Efficiency_fracwd;
  int SNIa_Mode;
  int SNIa_MassTransferOn;
  int SNIa_EnergyTransferOn;
  double SNIa_Energy_ERG;
  /* SNII parameters and flags */
  int SNII_MassTransferOn;
  double SNII_Factor_Hydrogen;
  double SNII_Factor_Helium;
  double SNII_Factor_Carbon;
  double SNII_Factor_Nitrogen;
  double SNII_Factor_Oxygen;
  double SNII_Factor_Neon;
  double SNII_Factor_Magnesium;
  double SNII_Factor_Silicon;
  double SNII_Factor_Iron;
  /* AGB flag */
  int AGB_MassTransferOn;
  /* Yield tables path */
  char YieldTablePath[STR_LENGTH];
#endif

#ifdef BG_SFR
#ifdef BG_SNII_THERMAL_FEEDBACK
  int SNII_EnergyTransferOn;

  double SNII_EnergyPerUnitMass_ERGperG;
  double SNII_EnergyPerUnitMass;
  double SNII_EnergyFraction;
  double SNII_Energy_ERG;
  double SNII_AvailableEnergyPerUnitMass;
#endif

#if defined(BG_STELLAR_EVOLUTION) || defined(BG_SNII_THERMAL_FEEDBACK)
  double SNII_MinMass_MSUN;
  double SNII_MaxMass_MSUN;
#endif

#ifdef BG_POPIII
#ifdef BG_POPIII_THERMAL_FEEDBACK
  int POPIII_EnergyTransferOn;

  double POPIII_EnergyPerUnitMass_ERGperG;
  double POPIII_EnergyPerUnitMass;
  double POPIII_EnergyFraction;
  double POPIII_AvailableEnergyPerUnitMass;

  double POPIII_SNII_Energy_ERG;
  double POPIII_PISN_Energy_ERG;

  double POPIII_PISN_MinMass_MSUN;
  double POPIII_PISN_MaxMass_MSUN;
#endif
  double POPIII_SNII_MinMass_MSUN;
  double POPIII_SNII_MaxMass_MSUN;
#endif

#endif /* BG_SFR */


  double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
				   NOT be balanced.  Each processor allocates memory for PartAllocFactor times
				   the average number of particles to allow for that */

  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

  double TopNodeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

#ifdef BG_SFR
  double StarAllocFactor;
#endif

#ifdef BG_THERMAL_FEEDBACK_MAKE_PARTICLES_ACTIVE
  short LowestOccupiedTimeBin;
#endif

  /* some SPH parameters */

  int DesNumNgb;                /*!< Desired number of SPH neighbours */
  double MaxNumNgbDeviation;    /*!< Maximum allowed deviation neighbour number */

  double ArtBulkViscConst;      /*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasU;		/*!< the same, but converted to thermal energy per unit mass */
  double MinEgySpec;            /*!< the minimum allowed temperature expressed as energy per unit mass */

#if defined(BG_SFR) || defined(BG_COOLING)
  double InitGasU_ERG;
  double MinGasU_ERG;
#endif

  /* some force counters  */

  long long TotNumOfForces;	/*!< counts total number of force computations  */

  long long NumForcesSinceLastDomainDecomp;	/*!< count particle updates since last domain decomposition */

  /* some variable for dynamic work-load adjustment based on CPU measurements */

  double Cadj_Cost;
  double Cadj_Cpu;

  /* system of units  */

  double UnitTime_in_s,		/*!< factor to convert internal time unit to seconds/h */
    UnitMass_in_g,		/*!< factor to convert internal mass unit to grams/h */
    UnitVelocity_in_cm_per_s,	/*!< factor to convert intqernal velocity unit to cm/sec */
    UnitLength_in_cm,		/*!< factor to convert internal length unit to cm/h */
    UnitPressure_in_cgs,	/*!< factor to convert internal pressure unit to cgs units (little 'h' still
				   around!) */
    UnitDensity_in_cgs,		/*!< factor to convert internal length unit to g/cm^3*h^2 */
    UnitCoolingRate_in_cgs,	/*!< factor to convert internal cooling rate to cgs units */
    UnitEnergy_in_cgs,		/*!< factor to convert internal energy to cgs units */
    UnitTime_in_Megayears,	/*!< factor to convert internal time to megayears/h */
    GravityConstantInternal,	/*!< If set to zero in the parameterfile, the internal value of the
				   gravitational constant is set to the Newtonian value based on the system of
				   units specified. Otherwise the value provided is taken as internal gravity
				   constant G. */
    G;				/*!< Gravity-constant in internal units */

  /* Cosmology */

  double Hubble;		/*!< Hubble-constant in internal units */
  double Omega0,		/*!< matter density in units of the critical density (at z=0) */
    OmegaLambda,		/*!< vaccum energy density relative to crictical density (at z=0) */
    OmegaBaryon,		/*!< baryon density in units of the critical density (at z=0) */
    HubbleParam;		/*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute
				 * physical values for cooling physics
				 */

  double BoxSize;		/*!< Boxsize in case periodic boundary conditions are used */

  /* Code options */

  int ComovingIntegrationOn;	/*!< flags that comoving integration is enabled */
  int PeriodicBoundariesOn;	/*!< flags that periodic boundaries are enabled */
  int ResubmitOn;		/*!< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;	/*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative
				   criterion */
  int TypeOfTimestepCriterion;	/*!< gives type of timestep criterion (only 0 supported right now - unlike
				   gadget-1.1) */
  int OutputListOn;		/*!< flags that output times are listed in a specified file */
  int CoolingOn;		/*!< flags that cooling is enabled */
  int StarformationOn;		/*!< flags that star formation is enabled */
  int StellarEvolutionOn;       /*!< flags that stellar evolution is enabled */


  /* parameters determining output frequency */

  int SnapshotFileCount;	/*!< number of snapshot that is written next */
  double TimeBetSnapshot,	/*!< simulation time interval between snapshot files */
    TimeOfFirstSnapshot,	/*!< simulation time of first snapshot files */
    CpuTimeBetRestartFile,	/*!< cpu-time between regularly generated restart files */
    TimeLastRestartFile,	/*!< cpu-time when last restart-file was written */
    TimeBetStatistics,		/*!< simulation time interval between computations of energy statistics */
    TimeLastStatistics;		/*!< simulation time when the energy statistics was computed the last time */
  int NumCurrentTiStep;		/*!< counts the number of system steps taken up to this point */

  /* Current time of the simulation, global step, and end of simulation */

  double Time,			/*!< current time of the simulation */
    TimeBegin,			/*!< time of initial conditions of the simulation */
    TimeStep,			/*!< difference between current times of previous and current timestep */
    TimeMax;			/*!< marks the point of time until the simulation is to be evolved */

  /* variables for organizing discrete timeline */

  double Timebase_interval;	/*!< factor to convert from floating point time interval to integer timeline */
  int Ti_Current;		/*!< current time on integer timeline */
  int Ti_nextoutput;		/*!< next output time on integer timeline */


#ifdef PMGRID
  int PM_Ti_endstep, PM_Ti_begstep;
  double Asmth[2], Rcut[2];
  double Corner[2][3], UpperCorner[2][3], Xmintot[2][3], Xmaxtot[2][3];
  double TotalMeshSize[2];
#endif

  int Ti_nextlineofsight;
#ifdef OUTPUTLINEOFSIGHT
  double TimeOfFirstLineOfSight;
  double TimeBetLineOfSight;
#endif
#ifdef ADAPTIVE_OUTPUT
  int AO_SubSamplingFactor;
  double AO_TimeOfFirstOutput;
  int NeedFileRefresh;
#endif
#ifdef BG_OUTPUT_GRID
  int Ti_nextgridoutput;
  double TimeOfFirstGridOutput;
  double TimeBetGridOutput;
#endif
#ifdef FOF
  int Ti_nextfof;
  double TimeOfFirstFOF;
  double TimeBetFOF;
  int FOFFileCount;
#endif

  /* variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;
  double CPU_Sum[CPU_PARTS];    /*!< sums wallclock time/CPU consumption in whole run */
  
  /* tree code opening criterion */

  double ErrTolTheta;		/*!< BH tree opening angle */
  double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */


  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
				   timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep,       /*!< minimum allowed timestep. Normally, the simulation terminates if the
                                  timestep determined by the timestep criteria falls below this limit. */ 
         MaxSizeTimestep;       /*!< maximum allowed timestep */

  double MaxRMSDisplacementFac; /*!< this determines a global timestep criterion for cosmological simulations
                                     in comoving coordinates.  To this end, the code computes the rms velocity
                                     of all particles, and limits the timestep such that the rms displacement
                                     is a fraction of the mean particle separation (determined from the
                                     particle mass and the cosmological parameters). This parameter specifies
                                     this fraction. */



  double CourantFac;		/*!< SPH-Courant factor */


  /* frequency of tree reconstruction/domain decomposition */


  double TreeDomainUpdateFrequency;	/*!< controls frequency of domain decompositions  */


  /* gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening
   * length)
   *
   * five groups of particles are supported 0=gas,1=halo,2=disk,3=bulge,4=stars
   */
  double MinGasHsmlFractional, /*!< minimum allowed SPH smoothing length in units of SPH gravitational
                                  softening length */
    MinGasHsml;                /*!< minimum allowed SPH smoothing length */


  double SofteningGas,    /*!< for type 0 */ 
    SofteningHalo,        /*!< for type 1 */ 
    SofteningDisk,        /*!< for type 2 */ 
    SofteningBulge,       /*!< for type 3 */ 
    SofteningStars,       /*!< for type 4 */ 
    SofteningBndry;       /*!< for type 5 */ 

  double SofteningGasMaxPhys,   /*!< for type 0 */ 
    SofteningHaloMaxPhys,       /*!< for type 1 */ 
    SofteningDiskMaxPhys,       /*!< for type 2 */ 
    SofteningBulgeMaxPhys,      /*!< for type 3 */ 
    SofteningStarsMaxPhys,      /*!< for type 4 */ 
    SofteningBndryMaxPhys;      /*!< for type 5 */ 

  double SofteningTable[6];  /*!< current (comoving) gravitational softening lengths for each particle type */
  double ForceSoftening[6];  /*!< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */


  /*! If particle masses are all equal for one type, the corresponding entry in MassTable is set to this
   *  value, * allowing the size of the snapshot files to be reduced
   */
  double MassTable[6];


  /* some filenames */
  char InitCondFile[100],
    OutputDir[100],
    SnapshotFileBase[100],
    EnergyFile[100],
    CpuFile[100],
    InfoFile[100], 
    TimingsFile[100], 
    RestartFile[100], 
    ResubmitCommand[100], 
    OutputListFilename[100];

  char RunLabel[200];
  /*! table with desired output times */
  double OutputListTimes[MAXLEN_OUTPUTLIST];

  int OutputListLength; /*!< number of times stored in table of desired output times */



#if defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS_HSML)
  double ReferenceGasMass;
#endif

  double OverDensThresh;
  double PhysDensThresh;
  double MaxPhysDensThresh;
  double OrigGasMass; 
#ifdef BG_SFR
  int SF_THRESH_MaxPhysDensOn;
  double SF_THRESH_MinOverDens;
  double SF_THRESH_MinPhysDens_HpCM3;
  double SF_THRESH_MaxPhysDens_HpCM3;
  double SF_THRESH_MaxTemp_K;
#ifdef BG_POPIII
  int POPIII_MassTransferOn;
  double POPIII_MetallicityThreshold_SOLAR;
  double POPIII_MetallicityThreshold;
  double POPIII_IMF_Exponent;
  double POPIII_IMF_MinMass_MSUN;
  double POPIII_IMF_MaxMass_MSUN;
  double POPIII_MinMass_MSUN;
  double POPIII_MaxMass_MSUN;
#endif
#endif

#ifdef DARKENERGY
  double DarkEnergyParam;	/*!< fixed w for equation of state */
#ifdef TIMEDEPDE
  char DarkEnergyFile[100];	/*!< tabelized w for equation of state */
#endif
#endif

#ifdef RESCALEVINI
  double VelIniScale;		/*!< Scale the initial velocities by this amount */
#endif


#ifdef BLACK_HOLES
  double SeedBHMassOverGasMass; 
  double BlackHoleAccretionFactor; /*!< alpha_0 */
  double BlackHoleAccretionSlope; /*!< beta */
  double BlackHoleEddingtonFactor;
  double BlackHoleFeedbackFactor; /*! epsilon_f */
  double BlackHoleRadiativeEfficiency; /*! epsilon_f */
  double BlackHoleCriticalEnergy;
  double BlackHoleMinHeatTemp;
  double BlackHoleNumberOfNeighboursToHeat;

#ifdef FOF
  double MinFoFSizeForNewSeed;	/*!< Halo particle number before new seed is put in */
#endif
#endif

#ifdef BG_SFR
  long long TotN_star;
  double GasFraction;
#endif

#if defined(BG_MOL_COOLING) || defined(BG_MOL_NETWORK)
  double InitAbundance_Hp;
  double InitAbundance_Hm;
  double InitAbundance_H2;
  double InitAbundance_H2p;
  double InitAbundance_Hep;
  double InitAbundance_Hepp;
  double InitAbundance_D;
  double InitAbundance_Dp;
  double InitAbundance_HD;
  double InitNumberDensity_e;
#endif

/* #ifdef BG_DUST */
/* #if defined(BG_DUST_DESTRUCTION_SUBLIMATION) || defined(BG_DUST_DESTRUCTION_SPUTTERING) */
/*   double InitDustGrainSize_nM; */
/* #endif */
/* #endif */
}
All;




/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern struct particle_data
{
  MyDouble Pos[3];   /*!< particle position at its current time */
  MyDouble Vel[3];   /*!< particle velocity at its current time */
  MyDouble Mass;     /*!< particle mass */
  MyIDType ID;

  union
  {
    MyFloat       GravAccel[3];		/*!< particle acceleration due to gravity */
    MyLongDouble dGravAccel[3];
  } g;
#ifdef PMGRID
  MyFloat GravPM[3];		/*!< particle acceleration due to long-range PM gravity force*/
#endif
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUTPOTENTIAL)
  union
  {
    MyFloat       Potential;		/*!< gravitational potential */
    MyLongDouble dPotential;
  } p;
#endif
  MyFloat OldAcc;			/*!< magnitude of old gravitational force. Used in relative opening
                                  criterion */

#if defined(EVALPOTENTIAL) && defined(PMGRID)
  MyFloat PM_Potential;
#endif

#ifdef BLACK_HOLES
  MyFloat BH_BirthTime;
#endif

#if defined (BLACK_HOLES) || defined (BG_SFR)
  MyFloat Hsml;	

  union
  {
    MyFloat       NumNgb;
    MyLongDouble dNumNgb;
  } n;
#endif

#ifdef BG_SFR
  //#ifndef LONGIDS
  unsigned int StarID;
  //#else
  //  unsigned long long StarID;
  //#endif
#endif

#ifdef ADAPTIVE_OUTPUT
  MyFloat StepsSinceLastOutput;
#endif

#ifdef BLACK_HOLES
  int SwallowID;
  MyFloat BH_Mass;
  MyFloat BH_Mdot;
  MyFloat NgbMassSum;
  int NumNgb;
  union
  {
    MyFloat BH_Density;
    MyLongDouble dBH_Density;
  } b1;
  union
  {
    MyFloat BH_Entropy;
    MyLongDouble dBH_Entropy;
  } b2;
  union
  {
    MyFloat BH_SurroundingGasVel[3];
    MyLongDouble dBH_SurroundingGasVel[3];
  } b3;
  union
  {
    MyFloat BH_accreted_Mass;
    MyLongDouble dBH_accreted_Mass;
  } b4;
  union
  {
    MyFloat BH_accreted_BHMass;
    MyLongDouble dBH_accreted_BHMass;
  } b5;
  union
  {
    MyFloat BH_accreted_momentum[3];
    MyLongDouble dBH_accreted_momentum[3];
  } b6;
  union
  {
    MyFloat BH_accreted_BHEnergy;
    MyLongDouble dBH_accreted_BHEnergy;
  } b7;
#ifdef REPOSITION_ON_POTMIN
  MyFloat BH_MinPotPos[3];
  MyFloat BH_MinPot;
#endif
#ifdef BH_THERMALFEEDBACK
  MyFloat BH_Energy;
#endif
#endif

  float GravCost;		/*!< weight factor used for balancing the work-load */

  int Ti_begstep;		/*!< marks start of current timestep of particle on integer timeline */
  int Ti_current;		/*!< current time of the particle */

  short int Type;		/*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
  short int TimeBin;

#if defined(ORDER_SNAPSHOTS_BY_ID) && !defined(SUBFIND)
  int     GrNr;
  int     SubNr;
#endif
}
 *P,				/*!< holds particle data on local processor */
 *DomainPartBuf;		/*!< buffer for particle data used in domain decomposition */


#ifdef BG_SFR
extern struct star_particle_data
{
  MyFloat InitialMass;
  MyFloat GasDensity;
  double StarBirthTime;

#ifdef BG_SNII_KINETIC_FEEDBACK
  MyFloat WindFlag;
#endif

#if defined(BG_SNII_KINETIC_FEEDBACK) || defined(BG_SNII_THERMAL_FEEDBACK)
  MyDouble NgbMassSum;
#endif

#ifdef BG_EXTRA_ARRAYS
  MyFloat MaximumEntropy;
  MyFloat MaximumTemperature;
  MyFloat TimeMaximumEntropy;
  MyFloat TimeMaximumTemperature;
#endif

#ifdef BG_STELLAR_EVOLUTION
  MyFloat Metals[BG_NELEMENTS];
  MyFloat Metallicity;
  MyFloat SolidAngleWeightSum;
  MyFloat SNIaRate;
#ifdef BG_METALSMOOTHING
  MyFloat MetalsSmoothed[BG_NELEMENTS];
  MyFloat MetallicitySmoothed;
#ifdef BG_SNIA_IRON
  MyFloat IronFromSNIaSmoothed;
#endif
#endif
#ifdef BG_SNIA_IRON
  MyFloat IronFromSNIa;
#endif
#ifdef BG_Z_WEIGHTED_REDSHIFT
  MyFloat MetallicityWeightedRedshift;
#endif
#endif /* BG_STELLAR_EVOLUTION */

  //#ifndef LONGIDS
  unsigned int PID;
  //#else
  //  unsigned long long PID;
  //#endif
}
*StarP,
*DomainStarBuf;
#endif /* BG_SFR */


/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
extern struct sph_particle_data
{
  MyDouble Entropy;             /*!< current value of entropy (actually entropic function) of particle */
  MyFloat EntropyPred;
  MyFloat Pressure;		/*!< current pressure */
  MyFloat VelPred[3];		/*!< predicted SPH particle velocity at the current time */

  MyFloat OldDivVel;              /*!< local velocity divergence at the previous time-step */
  MyFloat ArtViscParam;
  MyFloat ArtDiffParam;

  union
  {
    MyFloat       R;
    MyLongDouble dR;
  } aky;

  union
  {
    MyFloat       Div2U;
    MyLongDouble dDiv2U;
  } bky;

#ifdef PRESSURE_ENTROPY_SPH
  MyFloat EntropyVarPred;

  union
  {
    MyFloat       WeightedDensity;
    MyLongDouble dWeightedDensity;
  } cky;

  union
  {
    MyFloat       DhsmlPressure;
    MyLongDouble dDhsmlPressure;
  } dky;
#endif

  MyFloat Metallicity;
/* #ifdef BG_DUST */
/*   MyFloat Dusticity; */
/* #ifdef BG_DUST_DESTRUCTION_SUBLIMATION */
/*   MyFloat DustGrainSize; */
/* #endif */
/* #endif */

#if defined(BG_SFR) && defined(BG_STELLAR_EVOLUTION)
  MyDouble OldMass;			/*!< particle mass */
#endif

  union
  {
    MyFloat Density;		/*!< current baryonic mass density of particle */
    MyLongDouble dDensity;
  } d;

#if !defined(BLACK_HOLES) && !defined(BG_SFR)
  MyFloat Hsml;			/*!< current smoothing length */
  union
  {
    MyFloat NumNgb;
    MyLongDouble dNumNgb;
  } n;
#endif

  union
  {
    MyFloat DtEntropy;              /*!< rate of change of entropy */
    MyLongDouble dDtEntropy;
  } e;
  union
  {
    MyFloat HydroAccel[3];		/*!< acceleration due to hydrodynamical force */
    MyLongDouble dHydroAccel[3];
  } a;

  union
  {
    MyFloat CurlVel;	 	        /*!< local velocity curl */
    MyFloat Rot[3];		        /*!< local velocity curl */
    MyLongDouble dRot[3];
  } r;

  union
  {
    MyFloat DivVel;			/*!< local velocity divergence */
    MyLongDouble dDivVel;
  } v;

  union
  {
    MyFloat DhsmlDensityFactor;     /*!< correction factor needed in the equation of motion of the conservative
                                    entropy formulation of SPH */
    MyLongDouble dDhsmlDensityFactor; 
  } h;

/* #ifdef ALTERNATIVE_VISCOUS_TIMESTEP */
/*   MyFloat MinViscousDt; */
/* #else */
  MyFloat MaxSignalVel;           /*!< maximum signal velocity */
/* #endif */

#ifdef BG_SFR
  MyFloat Sfr;
  MyFloat OnEOS;
#endif

/* #if defined(BH_THERMALFEEDBACK) || defined(BG_SNII_THERMAL_FEEDBACK) */
/*    union */
/*       { */
/*       MyFloat Injected_Energy; */
/*       MyLongDouble dInjected_Energy; */
/*       } i; */
/* #endif */

#ifdef BG_EXTRA_ARRAYS
  MyFloat MaximumEntropy;
  MyFloat MaximumTemperature;
  MyFloat TimeMaximumEntropy;
  MyFloat TimeMaximumTemperature;
#endif

#ifdef BG_STELLAR_EVOLUTION
#ifdef BG_METALSMOOTHING
  MyFloat MetalsSmoothed[BG_NELEMENTS];
  MyFloat MetallicitySmoothed;
/* #ifdef BG_DUST */
/*   MyFloat DusticitySmoothed; */
/* #endif */
#ifdef BG_SNIA_IRON
  MyFloat IronFromSNIaSmoothed;
#endif
#endif

#ifdef BG_SNIA_IRON
  MyFloat IronFromSNIa;
#endif

#ifdef BG_Z_WEIGHTED_REDSHIFT
  MyFloat MetallicityWeightedRedshift;
#endif
#endif  /* BG_STELLAR_EVOLUTION */

#if defined(BG_STELLAR_EVOLUTION) || defined(BG_COOLING)
  MyFloat Metals[BG_NELEMENTS];
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK
  MyFloat WindFlag;
#endif

#ifdef BG_MOL_NETWORK
  MyFloat x_Hp;
  MyFloat x_Hm;
  MyFloat x_H2;
  MyFloat x_H2p;
  MyFloat x_Hep;
  MyFloat x_Hepp;
  MyFloat x_Dp;
  MyFloat x_HD;
  MyFloat n_e;

#ifdef BG_METALSMOOTHING
  MyFloat x_H2_Smoothed;
  MyFloat x_H2p_Smoothed;
  MyFloat x_HD_Smoothed;
#endif

  MyFloat Chem_dt;

/* #if defined(LW_BACKGROUND) || defined(LW_LOCAL) */
/*   MyFloat Kdiss_H2; */
/*   MyFloat Kdiss_Hm; */

/*   MyFloat LWRadiation_popII; */
/*   MyFloat LWRadiation_popIII; */
/* #endif */
#endif

#ifdef TIMESTEP_LIMITER
  short FlagFeedbackDtLimiter;
  int Ti_begstep;
#endif
}
  *SphP,				/*!< holds SPH particle data on local processor */
  *DomainSphBuf;			/*!< buffer for SPH particle data in domain decomposition */


extern peanokey *DomainKeyBuf;

/* global state of system 
*/
extern struct state_of_system
{
  double Mass,
    EnergyKin,
    EnergyPot,
    EnergyInt,
    EnergyTot,
    Momentum[4],
    AngMomentum[4],
    CenterOfMass[4],
    MassComp[6],
    EnergyKinComp[6],
    EnergyPotComp[6],
    EnergyIntComp[6], 
    EnergyTotComp[6], 
    MomentumComp[6][4], 
    AngMomentumComp[6][4], 
    CenterOfMassComp[6][4];
}
SysState, SysStateAtStart, SysStateAtEnd;



/* Various structures for communication 
 */

extern struct data_index
{
  int Task;
  int Index;
  int IndexGet;
}
 *DataIndexTable;		/*!< the particles to be exported are grouped
                                  by task-number. This table allows the
                                  results to be disentangled again and to be
                                  assigned to the correct particle */

extern struct data_nodelist
{
  int NodeList[NODELISTLENGTH];
}
*DataNodeList;		

extern struct gravdata_in
{
  MyFloat Pos[3];
#ifdef UNEQUALSOFTENINGS
  int Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  MyFloat Soft;
#endif
#endif
  MyFloat OldAcc;
  int NodeList[NODELISTLENGTH];
}
 *GravDataIn,			/*!< holds particle data to be exported to other processors */
 *GravDataGet;			/*!< holds particle data imported from other processors */


extern struct gravdata_out
{
  MyLongDouble Acc[3];
#ifdef EVALPOTENTIAL
  MyLongDouble Potential;
#endif
  int Ninteractions;
}
 *GravDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


extern struct potdata_out
{
  MyLongDouble Potential;
}
 *PotDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *PotDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */




/*! Header for the standard file format.
 */
extern struct io_header
{
  int npart[6];           /*!< number of particles of each type in this file */
  double mass[6];              /*!< mass of particles of each type. If 0, then the masses are explicitly
                                 stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                 /*!< time of snapshot file */
  double redshift;             /*!< redshift of snapshot file */
  int flag_sfr;           /*!< flags whether the simulation was including star formation */
  int flag_feedback;      /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];      /*!< total number of particles of each type in this snapshot. This can be
                                        different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;       /*!< flags whether cooling was included  */
  int num_files;          /*!< number of files in multi-file snapshot */
  double BoxSize;              /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;               /*!< matter density in units of critical density */
  double OmegaLambda;          /*!< cosmological constant parameter */
  double HubbleParam;          /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;    /*!< flags whether the file contains formation times of star particles */
  int flag_metals;        /*!< flags whether the file contains metallicity values for gas and star
                                 particles */
  unsigned int npartTotalHighWord[6];   /*!< High word of the total number of particles of each type */
  int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
  char fill[60];	       /*!< fills to 256 Bytes */
}
header;				/*!< holds header for snapshot files */



enum iofields
{ IO_POS,
  IO_VEL,
  IO_ID,
  IO_MASS,
  IO_U,
  IO_CR_C0,
  IO_CR_Q0,
  IO_CR_P0,
  IO_CR_E0,
  IO_CR_n0,
  IO_CR_ThermalizationTime,
  IO_CR_DissipationTime,
  IO_RHO,
  IO_WRHO,
  IO_iMass,
  IO_Zs,
  IO_CLDX,
  IO_HSML,
  IO_SFR,
  IO_AGE,
  IO_Z,
  IO_POT,
  IO_ACCEL,
  IO_DTENTR,
  IO_TSTP,
  IO_BFLD,
  IO_DBDT,
  IO_DIVB,
  IO_ABVC,
  IO_AMDC,
  IO_PHI,
  IO_COOLRATE,
  IO_CONDRATE,
  IO_BSMTH,
  IO_DENN,
  IO_EGYPROM,
  IO_EGYCOLD,
  IO_BHMASS,
  IO_BHMDOT,
  IO_BH_ENERGY,
  IO_AO,
  IO_BH_BIRTH_TIME,
  IO_MACH,
  IO_DTENERGY,
  IO_PRESHOCK_DENSITY,
  IO_PRESHOCK_ENERGY,
  IO_PRESHOCK_XCR,
  IO_DENSITY_JUMP,
  IO_ENERGY_JUMP,
  IO_CRINJECT,
  IO_BG_TEMP,
  IO_BG_ONEOS,
  IO_BG_SNII_KINETIC_FEEDBACK,
  IO_BG_METALS,
  IO_BG_METALS_SMOOTHED,
  IO_BG_METALLICITY,
  IO_BG_METALLICITY_SMOOTHED,
  IO_BG_METALLICITY_WEIGHTED_REDSHIFT,
  IO_BG_METALLICITY_WEIGHTED_POTENTIAL,
  IO_BG_INITIAL_MASS,
  IO_BG_IRON_FROM_SNIA,
  IO_BG_IRON_FROM_SNIA_SMOOTHED,
  IO_BG_MAX_ENTROPY,
  IO_BG_MAX_TEMPERATURE,
  IO_BG_TIME_MAX_ENTROPY,
  IO_BG_TIME_MAX_TEMPERATURE,
  IO_BG_STELLAR_AGE,

#ifdef SUBFIND
  IO_SUBFIND_MASS,
  IO_SUBFIND_MASSTYPE,
  IO_SUBFIND_POS,
  IO_SUBFIND_CMPOS,
  IO_SUBFIND_HALFMASS,
  IO_SUBFIND_HALFMASSPROJ,
  IO_SUBFIND_VMAXRAD,
  IO_SUBFIND_SPIN,
  IO_SUBFIND_VMAX,
  IO_SUBFIND_CMVEL,
  IO_SUBFIND_VELDISP,
  IO_SUBFIND_STELLARVELDISP,
  IO_SUBFIND_STELLARVELDISPHALFPROJ,
  IO_SUBFIND_PARTPOS,
  IO_SUBFIND_PARTVEL,
  IO_SUBFIND_PARTMASS,
  IO_M_MEAN,
  IO_M_CRIT,
  IO_M_TOPHAT,
  IO_R_MEAN,
  IO_R_CRIT,
  IO_R_TOPHAT,
#endif

  IO_BG_e,
  IO_BG_HII,
  IO_BG_Hminus,
  IO_BG_H2I,
  IO_BG_H2II,
  IO_BG_HeII,
  IO_BG_HeIII,
  IO_BG_DII,
  IO_BG_HD,

/*   IO_BG_DUST, */
/*   IO_BG_DUST_SMOOTHED, */
/*   IO_BG_DUST_SIZE, */

/*   IO_LW_popII, */
/*   IO_LW_popIII, */

  IO_LASTENTRY			/* This should be kept - it signals the end of the list */
};



/*
 * Variables for Tree
 * ------------------
 */

extern struct NODE
{
  MyFloat len;			/*!< sidelength of treenode */
  MyFloat center[3];		/*!< geometrical center of node */
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  MyFloat maxsoft;                /*!< hold the maximum gravitational softening of particles in the 
                                     node if the ADAPTIVE_GRAVSOFT_FORGAS option is selected */
#endif
  union
  {
    int suns[8];		/*!< temporary pointers to daughter nodes */
    struct
    {
      MyFloat s[3];		/*!< center of mass of node */
      MyFloat mass;		/*!< mass of node */
      unsigned int bitflags;	/*!< flags certain node properties */
      int sibling;		/*!< this gives the next node in the walk in case the current node can be used */
      int nextnode;		/*!< this gives the next node in case the current node needs to be opened */
      int father;		/*!< this gives the parent node of each node (or -1 if we have the root node) */
    }
    d;
  }
  u;
  int Ti_current;
}
 *Nodes_base,			/*!< points to the actual memory allocted for the nodes */
 *Nodes;			/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
				   gives the first allocated node */


extern struct extNODE
{
  MyLongDouble dp[3];
#ifdef FLTROUNDOFFREDUCTION
  MyFloat s_base[3];
  MyFloat len_base;
#endif
  MyFloat vs[3];
  MyFloat vmax;
  MyFloat divVmax;
  MyFloat hmax;			/*!< maximum SPH smoothing length in node. Only used for gas particles */
  int Ti_lastkicked;
  int Flag;
}
 *Extnodes, *Extnodes_base;


extern int MaxNodes;		/*!< maximum allowed number of internal nodes */
extern int Numnodestree;	/*!< number of (internal) nodes in each tree */


extern int *Nextnode;		/*!< gives next node in tree walk  (nodes array) */
extern int *Father;		/*!< gives parent node in tree (Prenodes array) */

#ifdef STATICNFW
extern double Rs, R200;
extern double Dc;
extern double RhoCrit, V200;
extern double fac;
#endif


#endif  /* ALLVARS_H */

