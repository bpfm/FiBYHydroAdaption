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

#include "allvars.h"



#ifdef PERIODIC
MyDouble boxSize, boxHalf;

#ifdef LONG_X
MyDouble boxSize_X, boxHalf_X;
#else
#endif
#ifdef LONG_Y
MyDouble boxSize_Y, boxHalf_Y;
#else
#endif
#ifdef LONG_Z
MyDouble boxSize_Z, boxHalf_Z;
#else
#endif
#endif


/*********************************************************/
/*  Global variables                                     */
/*********************************************************/


#ifdef FORCE_DOMAIN_DECOMPOSITION
double PartsPerSecond;
#endif


int ThisTask;			/*!< the number of the local processor  */

int NTask;			/*!< number of processors */

int PTask;			/*!< note: NTask = 2^PTask */

double CPUThisRun;		/*!< Sums CPU time of current process */

int NumForceUpdate;		/*!< number of active particles on local processor in current timestep  */

long long GlobNumForceUpdate;

int NumSphUpdate;		/*!< number of active SPH particles on local processor in current timestep  */

int MaxTopNodes;		/*!< Maximum number of nodes in the top-level tree used for domain decomposition */

int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from
				   initial conditions, 1 is resuming a run from a set of restart files, while 2
				   marks a restart from a snapshot file. */

int RestartSnapNum;

int *Exportflag;		/*!< Buffer used for flagging whether a particle needs to be exported to another process */

int *Exportnodecount;

int *Exportindex;

int *Send_offset, *Send_count, *Recv_count, *Recv_offset, *Sendcount_matrix;

int FirstActiveParticle;

int *NextActiveParticle;

int TimeBinCount[TIMEBINS];

int TimeBinCountSph[TIMEBINS];

int TimeBinActive[TIMEBINS];

int FirstInTimeBin[TIMEBINS];

int LastInTimeBin[TIMEBINS];

int *NextInTimeBin;

int *PrevInTimeBin;

#ifdef BG_SFR
double TimeBinSfr[TIMEBINS];
#ifdef BG_POPIII
double TimeBinPOPIIISfr[TIMEBINS];
#endif
#endif

#ifdef BLACK_HOLES
double TimeBin_BH_mass[TIMEBINS];

double TimeBin_BH_dynamicalmass[TIMEBINS];

double TimeBin_BH_Mdot[TIMEBINS];

double TimeBin_BH_Medd[TIMEBINS];
#endif



size_t AllocatedBytes;

size_t FreeBytes;

double CPU_Step[CPU_PARTS];

char CPU_Symbol[CPU_PARTS] =
  { '-', '*', '=', ';', '<', '[', '^', ':', '.', '~', '|', '+', '"', '/', '`', ',', '>', '@', '#', '&', '$',
  ']', '(', '?', ')', '1', '2', '3', '4', '5', '6', '7', '8'
};

char CPU_SymbolImbalance[CPU_PARTS] =
  { 'a', 't', 'u', 'v', 'b', 'w', 'd', 'r', 'h', 'm', 'n', 'l', 'o', 'p', 's', 'f', 'i', 'g', 'c', 'e', 'x',
  'y', 'z', 'A', 'I', 'W', 'T', 'V', 'B', 'C', 'D', 'E', 'F'
};

char CPU_String[CPU_STRING_LEN + 1];

double WallclockTime;		/*!< This holds the last wallclock time measurement for timings measurements */

int Flag_FullStep;		/*!< Flag used to signal that the current step involves all particles */


int TreeReconstructFlag;

int GlobFlag;


int NumPart;			/*!< number of particles on the LOCAL processor */

int N_gas;			/*!< number of gas particles on the LOCAL processor  */

long long Ntype[6];		/*!< total number of particles of each type */

int NtypeLocal[6];		/*!< local number of particles of each type */

#ifdef BG_SFR
int N_star;			/*!< number of star particles on the LOCAL processor  */
#endif


gsl_rng *random_generator;	/*!< the random number generator used */


#ifdef BG_SFR
int Stars_converted;		/*!< current number of star particles in gas particle block */
#endif


double TimeOfLastTreeConstruction;	/*!< holds what it says */

int *Ngblist;			/*!< Buffer to hold indices of neighbours retrieved by the neighbour search
				   routines */


double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;

int *DomainStartList, *DomainEndList;



double *DomainWork;

int *DomainCount;

int *DomainCountSph;

int *DomainTask;

int *DomainNodeIndex;

int *DomainList, DomainNumChanged;

#ifdef BG_SFR
int *DomainCountStar;
#endif

#ifdef SUBFIND
int GrNr;

int NumPartGroup;
#endif


peanokey *Key, *KeySorted;

struct topnode_data *TopNodes;

int NTopnodes, NTopleaves;




double RndTable[RNDTABLE];


/* variables for input/output , usually only used on process 0 */


char ParameterFile[100];	/*!< file name of parameterfile used for starting the simulation */

FILE *FdInfo,			/*!< file handle for info.txt log-file. */
 *FdEnergy,			/*!< file handle for energy.txt log-file. */
 *FdTimings,			/*!< file handle for timings.txt log-file. */
 *FdBalance,			/*!< file handle for balance.txt log-file. */
 *FdCPU;			/*!< file handle for cpu.txt log-file. */

#ifdef BG_SFR
FILE *FdSfr;			/*!< file handle for sfr.txt log-file. */
#ifdef BG_POPIII
FILE *FdPOPIIISfr;              /*!< file handle for popiii_sfr.txt log-file. */
#endif
#endif

#ifdef BG_STELLAR_EVOLUTION
FILE *FdMetGas;			/*!< file handle for metals_gas.txt log-file. */

FILE *FdMetStars;		/*!< file handle for metals_star.txt log-file. */

FILE *FdMetSF;			/*!< file handle for metals_sf.txt log-file. */

FILE *FdMetTot;			/*!< file handle for metals_tot.txt log-file. */

FILE *FdSNIa;			/*!< file handle for SNIa.txt log-file. */
#endif

#ifdef BG_SNII_KINETIC_FEEDBACK
FILE *FdKineticFeedback;	/*!< file handle for wind.txt log-file. */
#endif

#ifdef BG_SNII_THERMAL_FEEDBACK
FILE *FdThermalFeedback;	/*!< file handle for thermal_feedback.txt log-file. */
#endif

#ifdef BLACK_HOLES
FILE *FdBlackHoles;		/*!< file handle for blackholes.txt log-file. */
#ifdef VERBOSE_LOGFILES
FILE *FdBlackHolesDetails;
FILE *FdBlackHolesFeedback;
#endif

#ifdef FOF
FILE *FdBlackHolesSeeds;
#endif

#ifdef BH_DEBUG
FILE *FdBlackHolesEnergy;
#endif

#endif



#ifdef BG_SNII_KINETIC_FEEDBACK
MyFloat ExpectedMassLoading;
MyFloat ActualMassLoading;
MyFloat ActualMassLoading2ndLoop;
#endif



/*! table for the cosmological drift factors */
double DriftTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for gravitational forces */
double GravKickTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for hydrodynmical forces */
double HydroKickTable[DRIFT_TABLE_LENGTH];

void *CommBuffer;		/*!< points to communication buffer, which is used at a few places */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
struct global_data_all_processes All;




/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data *P,	/*!< holds particle data on local processor */
 *DomainPartBuf;		/*!< buffer for particle data used in domain decomposition */



#ifdef BG_SFR
struct star_particle_data *StarP, *DomainStarBuf;
#endif /* BG_SFR */


/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
struct sph_particle_data *SphP,	/*!< holds SPH particle data on local processor */
 *DomainSphBuf;			/*!< buffer for SPH particle data in domain decomposition */


peanokey *DomainKeyBuf;

/* global state of system 
*/
struct state_of_system SysState, SysStateAtStart, SysStateAtEnd;



/* Various structures for communication 
 */

struct data_index *DataIndexTable;	/*!< the particles to be exported are grouped
					   by task-number. This table allows the
					   results to be disentangled again and to be
					   assigned to the correct particle */

struct data_nodelist *DataNodeList;

struct gravdata_in *GravDataIn,	/*!< holds particle data to be exported to other processors */
 *GravDataGet;			/*!< holds particle data imported from other processors */


struct gravdata_out *GravDataResult,	/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


struct potdata_out *PotDataResult,	/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *PotDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */




/*! Header for the standard file format.
 */
struct io_header header;	/*!< holds header for snapshot files */





/*
 * Variables for Tree
 * ------------------
 */

struct NODE *Nodes_base,	/*!< points to the actual memory allocted for the nodes */
 *Nodes;			/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
				   gives the first allocated node */


struct extNODE *Extnodes, *Extnodes_base;


int MaxNodes;			/*!< maximum allowed number of internal nodes */

int Numnodestree;		/*!< number of (internal) nodes in each tree */


int *Nextnode;			/*!< gives next node in tree walk  (nodes array) */

int *Father;			/*!< gives parent node in tree (Prenodes array) */

#ifdef STATICNFW
double Rs, R200;

double Dc;

double RhoCrit, V200;

double fac;
#endif
