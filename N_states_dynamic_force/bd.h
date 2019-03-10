
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


#define FILE_NAME "bd.dat"

#define MAX 2.0e21
#define RMX RAND_MAX + 1

#define pi	3.14159265359

#define dt 	0.1	  	  // Simulation time step, ms
#define T 	100 	   	  // Total simulation time, s

#define kBT	4.1               // Kinetic energy RT, pN nm
//#define F	1		  // Stretching force, pN
#define R	500               // Bead radious, nm
#define eta	0.000001          // Viscosity of the medium, pN·nm^2/ms
#define D 	kBT/(6*pi*eta*R)  // Difusion constant, nm^2/ms


#define STATES	2
#define MAX_FORCE	3	// Max force for the fec
#define MIN_FORCE	0.02	// Min force for the fec

const float Lc[STATES] = {2000, 1950};	 // Contour lenght of the first state, nm
const float Lp[STATES] = {50, 50};	 	// Persistence lenght of each state, nm	

/* Transition matrix */
const float trans_p[STATES][STATES] =  {  {1-0.00002, 0.00002},  
					  {0.00002, 1-0.00002}
				  }; 


/* Main functions */
double extension(double Force, double Lc, double p);
double stiffness_z(double z, double Lc, double p);


/* Accesory funcions */
double randu(void);
double randn(double mean, double std);
 
