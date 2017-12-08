#if defined(BG_SFR) || defined(BG_COOLING)
extern char ElementNames[BG_NELEMENTS][EL_NAME_LENGTH];
extern double SolarMetallicity;
#endif

/* constants for metal dependent density threshold */
#define DENSITY_THRESHOLD_EXP   -0.32
#define DENSITY_THRESHOLD_COEFF (pow(pow(10, 20.75) / (3.06e21 * pow(0.1, DENSITY_THRESHOLD_EXP)), 2))

/* definition for max and min operators */
#define max(x, y) ((x > y) ? x : y)
#define min(x, y) ((x < y) ? x : y)

