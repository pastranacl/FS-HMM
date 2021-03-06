/**********************************************************************

  Brownian dynamics simulations
  
  This code performes Brownian dynamics simulations of bead in an
  harmonic potential based on the paper,
  
  "Quantitative Guidelines for Force Calibration through Spectral 
  Analysis of Magnetic Tweezers Data", Velthuis A.J.W, et al. 2010

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
  
  
  Compile with:
  
***********************************************************************/

#include "bd.h"


int main()
{
  int N_ts = (int)T*1000/dt;
  printf("Simulation steps:%d\n", N_ts);
 
  double *x = (double *)malloc( N_ts * sizeof(double *) ),
         *y = (double *)malloc( N_ts * sizeof(double *) ),
         *z = (double *)malloc( N_ts * sizeof(double *) ),
	 *t = (double *)malloc( N_ts * sizeof(double *) );
	 
  int *states = (int *)malloc( N_ts * sizeof(int *) );
  
  if(states == NULL)
  {
    printf("Not enough memory\n");
    exit(-1);
  }
  
  int state = 0;
  states[0] = 1;
  t[0] = 0;
    
  double eq_z = extension(F, Lc[state], Lp[state]),
	 eq_x = 0,
	 eq_y = 0;
	 
  x[0] = 0; 
  y[0] = 0; 
  z[0] = 0;
  
  double k_x = F/eq_z,
	 k_y = F/(eq_z + R),
	 k_z = stiffness_z(eq_z, F, Lc[state], Lp[state]);
   
	 
  // Brownian dynamics simulation
  
  double std_diff = sqrt(2*D*dt); 
  printf("Simulation start\n");
  srand(time(0));
  
  for(int i = 1; i < N_ts; i++)
  {    
      x[i] = x[i-1] - (x[i-1]*k_x*D*dt)/kBT + randn(0, std_diff);
      y[i] = y[i-1] - (y[i-1]*k_y*D*dt)/kBT + randn(0, std_diff);
      z[i] = z[i-1]  - (z[i-1] *k_z*D*dt)/kBT + + randn(0, std_diff);

      
      states[i] = state;
      t[i] = t[i-1] + dt;
     
      
      // Random change of state based on transition probabilities
      int state_trial =  rand() % STATES;
      
      if( randu() < trans_p[state][state_trial] )
      {
	state = state_trial;
	eq_z = extension(F, Lc[state], Lp[state]);
	k_z = stiffness_z(eq_z, F, Lc[state], Lp[state]);
	k_x = F/eq_z,
	k_y = F/(eq_z + R);
      }
  }
  
  
  
  // Save the data of the simulation
  //states[0] =  extension(F, Lc[0], Lp[0]);;
  
  FILE *open_file;
  open_file = fopen (FILE_NAME, "w");
  //fprintf(open_file, "t(s)\tstate\tx(nm)\ty(nm)\tz(nm)\n"); //Header
  fprintf(open_file, "%f\t%d\t%f\t%f\t%f\n", t[0],states[0],x[0],y[0],z[0]);
  
  for(int i = 0; i < N_ts; i++)
  {
    
    state = (int)states[i];
    eq_z = extension(F, Lc[state], Lp[state]);
    states[i] = eq_z;
    z[i] += eq_z;
    fprintf(open_file, "%f\t%d\t%f\t%f\t%f\n", t[i],states[i],x[i],y[i],z[i]);
  }
  fclose (open_file);
  
  free(x); 
  free(y); 
  free(z); 
  free(t);
  free(states);
  
  return 0;
}



// Returns WLC extension in nm
double extension(double Force, double Lc, double p)
{
  return Lc*( 1 - 0.5*sqrt( 1/((p*F)/kBT -1/32) )   );
}


// Returns the spring constant of a wlc polymer in pN/nm
double stiffness_z(double z, double Force, double Lc, double p)
{
  return (kBT/(Lc*p)) * ( 1 + 1/( 2 * pow(1 - z/Lc,3) ) );
}


// Generates a random number between 0 and 1
double randu(void)
{
  return (double)rand() / (double)RAND_MAX;
}

// Box-Muller polar form for gaussian distribution
double randn(double mean, double std)
{
  double x1, x2, w, y1, y2;

  do {
	x1 = 2.0 * randu() - 1.0;
	x2 = 2.0 * randu() - 1.0;
	w = x1 * x1 + x2 * x2;
  } while ( w >= 1.0 );
  w = sqrt( (-2.0 * log( w ) ) / w );
  y1 = x1 * w;
  y2 = x2 * w;
  return( mean + y1 * std );
}
 