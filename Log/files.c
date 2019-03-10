#include "files.h"


/**********************************************************************/
/*                              load_row 		              */
/* This function loads one of the rows obtained from the brownian     */
/* dynamics simulation and returns it				      */
/*								      */
/**********************************************************************/
double *load_row(char file_name[], int column, int MAX_ROWS)
{
  double *data = (double *)malloc(MAX_ROWS * sizeof(double *) );
  double tmp_val = 0;
  int r = 0;
  FILE* data_file;
  data_file = fopen(file_name,"r");
  
  for(int i = 0; i < MAX_ROWS; i++){
    for(int j = 0; j < 5; j++){
      fscanf(data_file,"%lf",&tmp_val);
      if( j == column ){
	data[r] = tmp_val;
	r++;
      }
    }
  }
  
  fclose(data_file);
  return data;
  
}

/**********************************************************************/
/*                              save row 		              */
/**********************************************************************/
// For the viterbin algorithm

