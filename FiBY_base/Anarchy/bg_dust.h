#define SIGMA 5.67e-5

#define THERMAL_SPUTTERING_TIME_SCALE (2000 * SEC_PER_YEAR)  /* yr */
#define THERMAL_SPUTTERING_TEMP_THRESH 1e6                   /* K */

#define SUBLIMATION_C0 2.154435e-8  /* cm [(1e-23)^(1/3)] */
#define SUBLIMATION_C1 1e15         /* s^(-1) */
#define SUBLIMATION_C2 7e4          /* K */

void dust_destruction(int j);
double get_dust_temperature(double T_gas, double n_gas, double T_rad);
double dust_temperature_evaluate(double x, void *params);
