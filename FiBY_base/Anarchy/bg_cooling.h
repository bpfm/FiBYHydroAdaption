#if defined(BG_COOLING) && !defined(BG_COOLING_OLD_TABLES)

/*
 * ----------------------------------------------------------------------
 * Cooling tables arrays
 * ----------------------------------------------------------------------
 */

extern float ****cooling_MetalsNetHeating;
extern float ****cooling_HplusHeNetHeating;
extern float ****cooling_HplusHeElectronAbundance;
extern float ****cooling_ThermalToTemp;

extern float ***cooling_SolarElectronAbundance;

#ifdef BG_COOLING_SHIELDING
extern float ***cooling_MetalsNetHeatingCollisional;
extern float ***cooling_HplusHeNetHeatingCollisional;
extern float ***cooling_HplusHeElectronAbundanceCollisional;
extern float ***cooling_ThermalToTempCollisional;

extern float **cooling_SolarElectronAbundanceCollisional;
#endif

extern float *cooling_Temp;
extern float *cooling_nH;
extern float *cooling_Therm;
extern float *cooling_HeFrac;
extern float *cooling_Redshifts;
extern float *cooling_SolarAbundances; 

extern char **cooling_ElementNames;
extern char **cooling_SolarAbundanceNames;

extern int cooling_N_He;
extern int cooling_N_Temp;
extern int cooling_N_nH;
extern int cooling_N_Elements;
extern int cooling_N_SolarAbundances;
extern int cooling_N_Redshifts;

extern int *ElementNamePointers;
extern int *SolarAbundanceNamePointers;

extern double *cooling_ElementAbundance_SOLAR;
extern double *cooling_ElementAbundance_SOLARM1;

#ifdef BG_MOL_COOLING
extern int mol_cooling_H2_N_Temp;
extern float *mol_cooling_H2_Temp;
extern float *mol_cooling_H2_Rate;

extern int mol_cooling_HD_N_Temp;
extern float *mol_cooling_HD_Temp;
extern float *mol_cooling_HD_Rate;
#endif

/*
 * ----------------------------------------------------------------------
 * Cooling routine constants
 * ----------------------------------------------------------------------
 */

#define STEFAN       7.5657e-15
#define COMP_COEFF   (4.0 * STEFAN * THOMPSON * BOLTZMANN / (ELECTRONMASS * C))
#define EV_TO_ERG    1.60217646e-12
#define COMPTON_RATE (COMP_COEFF * T_CMB0 * T_CMB0 * T_CMB0 * T_CMB0)

#endif
