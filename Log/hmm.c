/**********************************************************************

  Hidden Markov Model Analysis
  
  Analysis of a two-state model based in the energy of the WLC for 
  looped and unlooped states.
    
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
  gcc -lm -std=c99 hmm.c energy.c -o hmm
  
***********************************************************************/

#include "hmm.h"
#include "energy.h"
#include "files.h"

/* Physical Parameters of each state */

#define F 	8.89			// Constant force, pN

// Contour lenght of each state
const double Lc[STATES]	  	=	{106.657819034744, 120.97720620033, 135.296593365918, 149.615980531505, 163.935367697092, 178.254754862679, 192.574142028265, 206.893529193852, 221.212916359439};  

// Persistence lenght of each state
const double p[STATES]		= 	{0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};   


/* Markov Chain Parameters */

// Transition matrix
double trans_m[STATES][STATES]	=   { {0,0,0,0,0,0,0,0,0},
				      {0,0,0,0,0,0,0,0,0},
				      {0,0,0,0,0,0,0,0,0},
				      {0,0,0,0,0,0,0,0,0},
				      {0,0,0,0,0,0,0,0,0},
				      {0,0,0,0,0,0,0,0,0},
				      {0,0,0,0,0,0,0,0,0},
				      {0,0,0,0,0,0,0,0,0},
				      {0,0,0,0,0,0,0,0,0}, };    
				      
// Initial probabilities (1.0)				 
double init_p[STATES]	=   {0.111, 0.111, 0.111, 0.111, 0.111, 0.111, 0.111, 0.111, 0.111};	                 


/***********************************************************************/


int main()
{
  
  /********* Viterbi algorithm test *************/
  const int N=1428404;
  double *data = load_row("bd_26pts.dat", 4, N);
   
  double *tm = baum_welch(data, N);
 
  //Replace transition matrix
  for(int i = 0; i<STATES; i++)
  {
    for(int j = 0; j < STATES; j++)
    {
      printf("  %f  ", tm[i*STATES + j]);
      trans_m[i][j] = tm[i*STATES + j];
    }
    printf("\n");
  }
  
  
   
  double *states = viterbi(data, N);
  
  // Save the data
  FILE *open_file;
  open_file = fopen ("viterbi.dat", "w");
  
  for(int i = 0; i<N; i++)
    fprintf(open_file, "%f\n", states[i]);
  
  fclose(open_file);

  
  // Free the memory 
  free(states);
  free(data);
  
  
  return 0;
}


/**********************************************************************/
/*                              viterbi			              */
/* Implementation of the Viterbi algorithm (in log form). 	      */
/* The function returns the better state sequence after backtracking  */
/* converted into equilibrium extensions for each state.	      */
/*								      */
/* The input data are the extension and the # of data points N.	      */
/**********************************************************************/
double *viterbi (double *z, int N)
{
  
  /****************************************************************/
  // Allocate memory  
  // For delta variable, STATES-rows x 2-COLUMNS for delta
  double *delta = (double *)malloc( STATES * N * sizeof(double *) ); 
  
   // For the backstrack variable
  int *psi = (int *)malloc( STATES * N * sizeof(int *) );
  if (delta == NULL)
    printf("Not enough memory for delta at viterbi\n");
  if(psi == NULL)
    printf("Not enough memory for psi at viterbi\n");
   
  
  /****************************************************************/
  
  // Calculates equilibrium extension of each state
  double z_eq[STATES];
  for(int i = 0; i < STATES; i++) 
    z_eq[i] = extension(F, Lc[i], p[i]);
  
  // Calculates the width/stiffness of each state
  double k[STATES];
  for(int i = 0; i < STATES; i++) 
    k[i] = stiffness(z_eq[i], F, Lc[i], p[i]);
  
  // Calculates the partition function for each state
  double Z[STATES];
  for(int i = 0; i < STATES; i++) 
    Z[i] = calc_Z(F, Lc[i], p[i]);

  ///////////////////////////////////////////////////////////////////////////
  //				VITERBI ALGORITHM                	   //
  ///////////////////////////////////////////////////////////////////////////
    
  printf("- Starting Viterbi\n");
    
  /* 1.- Initialization */
  for(int i = 0; i < STATES; i++)
  {
    delta[i] = log(init_p[i]) + log(prob_ext_wlc(z[0], z_eq[i], k[i]) / Z[i]);
    psi[i] = 0;
  }   
  
  /* 2.- Recursion */
  double    trans_prob =  -DBL_MAX,
	   max_trans_p =  -DBL_MAX;
  
  for( int t = 1; t < N; t++)
  {
    for( int j = 0; j < STATES; j++)
    {
      max_trans_p = -DBL_MAX; 
      for( int i = 0; i < STATES; i++ )
      {
	trans_prob = delta[(t-1)*STATES + i] + log(trans_m[i][j]);
	// Max state 
	if( trans_prob > max_trans_p ){
	  delta[t*STATES + j] = trans_prob + log( prob_ext_wlc(z[t], z_eq[j], k[j]) / Z[j] );
	  psi[t*STATES + j] = i;
	  max_trans_p = trans_prob;
	}
      }// i 
    } 
  }
  
  /* 3.- Termination */
  // In this case there is not an special state at the end.
  int *q = (int *)malloc( N * sizeof(int *) );
  double *z_states = (double *)malloc( N * sizeof(double *) );
  
  double mx_dt = -DBL_MAX;
  
  for(int i = 0; i < STATES; i++)
  {
    if( delta[(N-1)*STATES + i] > mx_dt ){
      mx_dt = delta[(N-1)*STATES + i];
      q[N-1] = i;
    }
  }
  printf("  log(P*) = %e\n", mx_dt);
  
  // 4.- Backtracking
  for(int t = N-2; t > 0; t--){
    // Backtrack position. Assing the expected extension for the state 
    q[t] = psi[STATES*(t+1) + q[t+1]]; 
    z_states[t] = z_eq[q[t]]; 
  }
  
  free(delta);
  free(psi);
  free(q);

  return z_states;
}



/**********************************************************************/
/*                           forward_backward	                      */
/*								      */
/**********************************************************************/
double *baum_welch(double *z, int N)
{ 
  
  printf("- Starting Baum-Welch likelihood maximization.\n");
  
  // Calculates equilibrium extension of each state
  double z_eq[STATES];
  for(int i = 0; i < STATES; i++) 
    z_eq[i] = extension(F, Lc[i], p[i]);
  
  // Calculates the width/stiffness of each state
  double k[STATES];
  for(int i = 0; i < STATES; i++) 
    k[i] = stiffness(z_eq[i], F, Lc[i], p[i]);
  
  // Calculates the partition function for each state
  double Z[STATES];
  for(int i = 0; i < STATES; i++) 
    Z[i] = calc_Z(F, Lc[i], p[i]);
  
  // Creates the transition matrixes and assigns init. prob. to every value
  double *a = (double *)malloc(STATES*STATES * sizeof(double*));
  double *a_opt = (double *)malloc(STATES*STATES * sizeof(double*));
 
  for(int i = 0; i < STATES; i++)
  {
    for(int j =0; j < STATES; j++)
    {
      if(i==j)
      {
	a[i*STATES + j] = 1-STATES*INIT_PROB_FB;
	a_opt[i*STATES + j] = 1-STATES*INIT_PROB_FB;
      }
      else
      {
	a[i*STATES + j] = INIT_PROB_FB;
	a_opt[i*STATES + j] = INIT_PROB_FB;
      }
    }
  }

  // Allocates the forward and backward variables 
  double *alpha = (double *)malloc( STATES * N * sizeof(double *) ); 
  double *beta = (double *)malloc( STATES * N * sizeof(double *) );
  
  // Allocates the variable for the sum of the gamma variable
  double *sum_gamma = (double *)malloc( STATES * sizeof(double *) ); 
  
  
  // Scaling variable
  double *c = (double *)malloc( N * sizeof(double *) ); c[0]=0;
  double *s = (double *)malloc( N * sizeof(double *) ); s[0]=0;
  
  double p = 0;
  
  
  if (c == NULL)
    printf("Not enough memory for Forward-Backward algorithm\n");
  
  int its = 0;
  
  do{
    
    memcpy(a, a_opt, sizeof(a)); // Transition for the next iteration updated
    
    //////			FORWARD ALGORITHM                	   //////
    /* 1.- Initialization */
    for(int i = 0; i < STATES; i++)
    {
      alpha[i] = init_p[i] * prob_ext_wlc(z[0], z_eq[i], k[i]) / Z[i];
      c[0] += alpha[i];
    }
    
    // Scaling of forward for t=0;
    for(int i = 0; i < STATES; i++)
      alpha[i] /= c[0];
    
    /* 2.- Induction */
    for( int t=1; t < N; t++ ){
      c[t] = 0;
      for(int j = 0; j < STATES; j++)
      {
	alpha[t*STATES + j] = 0;
	for(int i = 0; i < STATES; i++) 
	{ 
	  alpha[t*STATES + j] += alpha[(t-1)*STATES + i] * a[i*STATES +j] * prob_ext_wlc(z[t], z_eq[j], k[j]) / Z[j];
	}
	c[t] += alpha[t*STATES + j]; // Scaling variable
      }
      // Scaling the forward variable
      for(int j = 0; j < STATES; j++)
	alpha[t*STATES + j] /= c[t];
    }
  
    /*  Calculate log[p(O | lambda)] */
    p = 0;
    for(int t = 0; t < N; t++ )
      p += log(c[t]);
    p = -p;
    
    //////			BACKWARD ALGORITHM                	 //////
    
   /* 1.- Initialization. For all states Beta_T = 1*/
   for( int i = 0; i < STATES; i++)
     beta[(N-1)*STATES +i] = 1;
   
    /* 	2.- Induction 	*/
    for( int t=(N-2); t >= 0; t-- )
    {
      s[t]=0;
      for(int i = 0; i < STATES; i++)
      {
	beta[t*STATES + i] = 0;
	for(int j = 0; j < STATES; j++)
	{
	  beta[t*STATES + i] += a[i*STATES +j] * beta[(t+1)*STATES + j] * prob_ext_wlc(z[t+1], z_eq[j], k[j]) / Z[j];
	}
	s[t] += beta[t*STATES + i]; // Scaling variable
      }
      // Scaling the backward variable
      for(int i=0;i<STATES;i++)
	beta[t*STATES + i] /= s[t];
    }
   
    /***************************** EXPECTATION ********************************/
    
    ////////////////////////////////////////////////////////////////////////////
    //				BAUM-WELCH PARAMETERS               	     //
    ///////////////////////////////////////////////////////////////////////////
    
   for(int i=0; i<STATES; i++)
   {
      // Calculates the sum of all epsilon_(i,j)
      for(int j=0; j<STATES; j++)
      {
	a_opt[i*STATES+j] = 0;
	for(int t=0; t<(N-1); t++)
	{
	  // This calculates the transition matrix summing all epsilon_(i,j) for all t (NON NORMALIZED)
	  a_opt[i*STATES+j] += alpha[t*STATES+i]*a[i*STATES +j]*beta[(t+1)*STATES +j]*prob_ext_wlc(z[t+1], z_eq[j], k[j]) / Z[j];
	} 
      }
      
      // Calculates the normalization factor, ie the Sum  of epsilon for every  = gamma(i)
      sum_gamma[i] = 0;
      for(int j=0; j < STATES; j++)
      {
	for(int t=0; t < (N-1); t++)
	  sum_gamma[i] += alpha[t*STATES+i]*a[i*STATES +j]*beta[(t+1)*STATES +j]*prob_ext_wlc(z[t+1], z_eq[j], k[j]) / Z[j];
      }
      
      /***************************** MAXIMIZATION ********************************/
      // Normalize the probabilities a_ij 
      // S ^T epsilon(i,j) / S^T gamma(i)
      for( int j = 0; j < STATES; j++)
	a_opt[i*STATES + j] /= sum_gamma[i];
      
   } 
   
   //printf("%f\n", matrx_diff(a, a_opt, STATES));
   its++;
   
  }while( ( THR_FOR_BACK < matrx_diff(a, a_opt, STATES) ) && ( its < MAX_ITS_FOR_BACK ) );
  
  
  
  free(alpha);
  free(beta);
  free(sum_gamma);
  free(c);
  free(s);
  
  printf("  Forward-backward completed after %d iterations with log[p(O|L)] = %e\n", its, -p);
  
  
  
  return a_opt;
}




/**********************************************************************/
/*                           forward		                      */
/* Here is wrote the forward algorithm to determine P( O | Lambda)    */
/* Because it is not scaled it underflows and for this particular     */
/* application the porpose in merely didactic. 		              */
/* 								      */
/* NOTE: This code is not correct? */
/**********************************************************************/
double forward(double *z, int T)
{
  // Allocating memory
  double *alpha = (double *)malloc( STATES * T * sizeof(double *) ); 
  if (alpha == NULL)
    printf("Not enough memory for the forward variable alpha at forward\n");

  // 0.- Calculates physical parameters
  // Calculates equilibrium extension of each state
  double z_eq[STATES];
  for(int i = 0; i < STATES; i++) 
    z_eq[i] = extension(F, Lc[i], p[i]);
  
  // Calculates the width/stiffness of each state
  double k[STATES];
  for(int i = 0; i < STATES; i++) 
    k[i] = stiffness(z_eq[i], F, Lc[i], p[i]);
  
  // Calculates the partition function for each state
  double Z[STATES];
  for(int i = 0; i < STATES; i++) 
    Z[i] = calc_Z(F, Lc[i], p[i]);

  ///////////////////////////////////////////////////////////////////////////
  //				FORWARD ALGORITHM                	   //
  ///////////////////////////////////////////////////////////////////////////
    
  printf("Starting Forward algorithm\n");
    
  /* 1.- Initialization */
  for(int i = 0; i < STATES; i++)
    alpha[i] = init_p[i] * prob_ext_wlc(z[0], z_eq[i], k[i]) / Z[i];
     
  // 2.- Induction
  for(int t = 1; t < T; t++){
    for(int j = 0; j<STATES; j++){
      for(int i = 0; i<STATES; i++){
	alpha[t*STATES + j] += alpha[(t-1)*STATES + i] * trans_m[i][j];   
      }
      alpha[t*STATES + j] *= prob_ext_wlc(z[t], z_eq[j], k[j]) / Z[j];
    }
  }
  
  // 3.- Termination
  double P = 0;
  for(int i = 0; i<STATES; i++)
    P += alpha[(T-1)*STATES + i];
  
  return P;
}


/**********************************************************************/
/*                              matrx_diff 		              */
/* This function returns the cummulative difference of all the        */
/* elements of two square matrices M and N with K x K elements        */
/**********************************************************************/
double matrx_diff(double *M, double *N, int K)
{
  double d = 0;
  for(int i = 0; i < K; i++ )
  {
    for( int j = 0; j < K; j++)
       d += fabs(M[i*K + j]-N[i*K +j]);
  }
  return d;
}