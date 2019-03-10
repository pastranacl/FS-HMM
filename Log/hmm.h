#include <stdio.h>
#include <stdlib.h>
#include <string.h> // for memcpy
#include <math.h>
#include <float.h>

#define STATES 	 		9		// Number of hidden states

#define INIT_PROB_FB		3.8462e-04	// Initial probababilites for the forward-backward
#define MAX_ITS_FOR_BACK 	50  		// Max numbef of iteractions of the forward_backward
#define THR_FOR_BACK 		0.0001		// Threshold for stopping the forward-backward optimization of the transition matrix

double *viterbi (double *extension, int N);

double forward(double *z, int T);

double *baum_welch(double *z, int N);
double matrx_diff(double *M, double *N, int K);





