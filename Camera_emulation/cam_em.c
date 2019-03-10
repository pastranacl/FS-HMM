#include <stdio.h>
#include <stdlib.h> 
#include <unistd.h>


#define ACQ_FREQ 	200.0f   		// Acquistion frequency in Hertz
#define DT_CAM  	1.0/ACQ_FREQ		// Time step (s)
#define DT_CAM_MS	DT_CAM*1000   		// Time step (ms)

#define FILE_NAME	"bd.dat"
#define CAM_FILE	"bd_cam.dat"

double *load_row (char file_name[], int column, int MAX_ROWS);


/* Call as cam_em #lines */
int main (int argc, char *argv[])
{
  // The calling to the program must include the number of data points
  // The data must be provided 
  if( argc < 2 ) printf("One argument must be provided");
  
  int N =  atoi(argv[1]); // atoi converts char to int (lib unistd)
  
  double *time = load_row(FILE_NAME, 0, N);
  double *states = load_row(FILE_NAME, 1, N);
  double *x = load_row(FILE_NAME, 2, N);
  double *y = load_row(FILE_NAME, 3, N);
  double *z = load_row(FILE_NAME, 4, N);
  
  
  double sim_dt = (time[N-1]-time[N-2])/1000; // Reduction to seconds 
  
  // How many points of the simulation are blurred given the shutter time
  int pts_cam =  DT_CAM/sim_dt;
  
  
  
  int cam_data_pts = (int)N/pts_cam;
  int n = 1;
  
  double *shutter_data_t = (double *)malloc(cam_data_pts * sizeof(double));
  double *shutter_states = (double *)malloc(cam_data_pts * sizeof(double));
  double *shutter_data_x = (double *)malloc(cam_data_pts * sizeof(double));
  double *shutter_data_y = (double *)malloc(cam_data_pts * sizeof(double));
  double *shutter_data_z = (double *)malloc(cam_data_pts * sizeof(double));
  
  shutter_data_t[0] = 0;
  shutter_states[0] = 0;
  shutter_data_x[0] = 0;
  shutter_data_y[0] = 0;
  shutter_data_z[0] = z[0];
  

  
  // Blurring the image with the mean position
  int t_check = 0;
  for(int t = 0; t < N; t++ )
  {
    
    if( t_check == (pts_cam) )
    {
      //printf(" t =  %d\n", t);
      shutter_data_x[n] /= pts_cam;
      shutter_data_y[n] /= pts_cam;
      shutter_data_z[n] /= pts_cam;
      shutter_states[n] = states[t];
      n++;
      t_check=0;
      shutter_data_t[n] = time[t];
      shutter_data_x[n] = 0;
      shutter_data_y[n] = 0;
      shutter_data_z[n] = 0;
      
    }
    else
    {
      t_check++;
       shutter_data_x[n] += x[t];
       shutter_data_y[n] += y[t];
       shutter_data_z[n] += z[t];
    }
    
  }
  
  // Save the data_file
  FILE *open_file;
  open_file = fopen(CAM_FILE, "w");
  //fprintf(open_file, "t(s)\tstate\tx(nm)\ty(nm)\tz(nm)\n"); //Header
  for(int i = 0; i < n-1; i++)
     fprintf(open_file, "\%f\t%f\t%f\t%f\t%f\n", shutter_data_t[i], shutter_states[i],shutter_data_x[i],shutter_data_y[i],shutter_data_z[i]);
  
  fclose(open_file);
  
  // Free memory
  free(time);
  free(states);
  free(x);
  free(y);
  free(z);
  free(shutter_data_t);
  free(shutter_data_x);
  free(shutter_data_y);
  free(shutter_data_z);
  return 0;
}


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

