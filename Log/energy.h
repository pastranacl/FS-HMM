#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define kBT	4.1		// Thermal energy at room temperature,  pN nm 
#define dz      0.001		// Integration step

// Statistical mechanics functions
double calc_Z(double F, double Lc, double p);
double extension(double F, double Lc, double p);
double stiffness(double z, double F, double Lc, double p);

double prob_ext_wlc(double z, double eq_z, double k);
double prob_dz_wlc(double z, double k);

double energy_z(double z, double z_eq, double k);
double energy_dz(double var_z, double k);



// Integration functions
double trapz(double a, double b, double z_eq, double k);

//double simpson(double a, double b, double z_eq, double k);
double simpson(double a, double b, double k);

double trap_int(double *x, double *y, int N); // Not used