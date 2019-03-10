 /**********************************************************************

  Hidden Markov Model Analysis
  
  Energy.c: This file contains the functions necessary for the 
  determination of emision probabilities in the Hidden Markov Model.
  These probabilities are calculated from the equilibrium position
  at a certain extension with a Gaussian distribution.
    
  Copyright (C) 2016, Cesar Lopez Pastrana
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
  
***********************************************************************/

#include "energy.h"




/**********************************************************************/
/*                            calc_Z 		                      */
/* This function calculates the partition function for the extension  */
/* energies the worm like chain at a given force with contour lenght  */
/* Lc and persistence length p					      */
/*							              */
/**********************************************************************/
/*
double calc_Z(double F, double Lc, double p)
{
  double z_eq   = extension(F, Lc, p);
  double k 	= stiffness(z_eq, F, Lc, p);
  double Z = 0;
  
  // Note that the with this method the energy is for a certain extension
  // is calculated twice 
  for( double i = 0; i < Lc-dz; i+= dz )
    Z += simpson(i, i + dz, z_eq, k);
  return Z;
}
*/
 
double calc_Z(double F, double Lc, double p)
{
  double Z = 0;
  double z_eq   = extension(F, Lc, p);
  double k 	= stiffness(z_eq, F, Lc, p);  

  for( double i = 0; i < Lc-dz; i+= dz )
    Z += simpson(i, i + dz, k);
  return Z;
}



/**********************************************************************/
/*                              stiffnes 		              */
/* The function returns the stifness of the DNA at a certain force    */
/* and extension. The equation comes from the dF/dz of Marko's WLC    */
/* interpolation formula. 					      */
/*							              */
/**********************************************************************/
double stiffness(double z, double F, double Lc, double p)
{
  return (kBT/(Lc*p)) * ( 1 + 1/( 2 * pow(1 - z/Lc,3) ) );
}


/**********************************************************************/
/*                              prob_ext_wlc 		              */
/* Determines the Boltzmann factor (non-normalized)of a certain       */
/* extension with stretched force F, with parameters Lc and p for     */
/* persistence and contour lenghts respectively.                      */
/*							              */
/**********************************************************************/
double prob_ext_wlc(double z, double eq_z, double k)
{
  return exp(-energy_z(z, eq_z, k)/kBT);
}


/**********************************************************************/
/*                              prob_dz_wlc 		              */
/* Determines the Boltzmann factor (non-normalized)of a certain       */
/* out of equilibrium with stretched force F, with parameters Lc and  */
/* p for persistence and contour lenghts respectively.                */
/*							              */
/**********************************************************************/
double prob_dz_wlc(double z, double k)
{
  return exp(-energy_dz(z, k)/kBT);
}



/**********************************************************************/
/*                              energy_z 		              */
/* Calculates the  energy at a certain extension out of equilibrium   */
/* based on a harmonic potential with extension z and stiffness k     */
/*       							      */
/**********************************************************************/
double energy_z(double z, double z_eq, double k)
{
  return 0.5*k*(z - z_eq)*(z - z_eq);
}


/**********************************************************************/
/*                              energy_dz 		              */
/* Calculates the  energy at a certain dz out of the equilibrium      */
/* extension. This is used to compute the partition function 	      */
/**********************************************************************/
double energy_dz(double var_z, double k)
{
  return 0.5*k*var_z*var_z;
}


/**********************************************************************/
/*                             extension 		              */
/* The function returns the expected extension for a molecule at      */
/* stretching force F, with parameters Lc and p for persistence       */
/* and contour lenghts respectively using the Moroz and Nelson WLC.   */
/*								      */
/**********************************************************************/
double extension(double F, double Lc, double p)
{
  return Lc*( 1 - 0.5*sqrt( 1/((p*F)/kBT -1/32) )   );
}




/**********************************************************************/
/*       							      */
/*       		INTEGRATION FUNCTIONS			      */
/*       							      */
/**********************************************************************/



/**********************************************************************/
/*                             trapz 	               	              */
/* Calculates the define integral between a and b of the Boltzmann    */
/* factor for the energy of the wlc.                                  */
/*       							      */
/**********************************************************************/
double trapz(double a, double b, double z_eq, double k)
{
  return (b-a)*0.5*(prob_ext_wlc(a, z_eq, k) + prob_ext_wlc(b, z_eq, k));
}



/**********************************************************************/
/*                             simpson 	               	              */
/* Calculates the define integral between a and b of the Boltzmann    */
/* factor for the energy of the wlc.                                  */
/*       							      */
/**********************************************************************/
/*double simpson(double a, double b, double z_eq, double k)
{
    return ( (b-a)/6 ) * ( prob_ext_wlc(a, z_eq, k) + prob_ext_wlc(b, z_eq, k)  + 4 * prob_ext_wlc( (a+b)/2, z_eq, k) );
}*/


double simpson(double a, double b, double k)
{
    return ( (b-a)/6 ) * ( prob_dz_wlc(a, k) + prob_dz_wlc(b, k)  + 4 * prob_dz_wlc( (a+b)/2, k) );
}

/**********************************************************************/
/*                              trap_int 		              */
/* Calculates the integral of the function values in x, and y using   */
/* the trapezoidal rule.                                              */
/*       							      */
/**********************************************************************/
double trap_int(double *x, double *y, int N)
{
  double cdf = 0;
  for(int i = 1; i < N; i++)
    cdf += 0.5 * (x[i+1] * x[i]) * (y[i-1] + y[i]);
  return cdf;
}
