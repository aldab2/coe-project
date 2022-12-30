#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>

#define MAX_ITER 1000
#define MAX 100 //maximum value of the matrix element
#define TOL 0
#define N 3
#define THREAD_COUNT 4
// Generate a random float number with the maximum value of max
float rand_float(int max){
  return ((float)rand()/(float)(RAND_MAX)) * max;
}


// Allocate 2D matrix
void allocate_init_2Dmatrix(float ***mat, int n, int m){
  int i, j;
  *mat = (float **) malloc(n * sizeof(float *));
  for(i = 0; i < n; i++) {
    (*mat)[i] = (float *)malloc(m * sizeof(float));
    for (j = 0; j < m; j++)
      (*mat)[i][j] = rand_float(MAX);
  }
} 

// solver
void solver(float ***mat, int n, int m){
  float diff = 0, temp;
  int done = 0, cnt_iter = 0, i, j;

  while (!done && (cnt_iter < MAX_ITER)){
    diff = 0;
   

#pragma omp parallel for num_threads(THREAD_COUNT) private(temp,j) reduction(+:diff)
    for(int i =0;i<n;i++){
        for(int j = i%2==0?0:1;j<(n/2)+2;j+=2){
          if(i==0 || j==0 || i==n-1 || j==n-1)
          continue;
          	temp = (*mat)[i][j];
	(*mat)[i][j] = 0.2 * ((*mat)[i][j] + (*mat)[i][j - 1] + (*mat)[i - 1][j] + (*mat)[i][j + 1] + (*mat)[i + 1][j]);
	diff += fabsf((*mat)[i][j] - temp);
        }
    }
#pragma omp parallel for num_threads(THREAD_COUNT) private(temp,j) reduction(+:diff)
    for(int i =0;i<n;i++){
        for(int j = i%2==1?0:0;j<(n/2)+2;j+=2){
          if(i==0 || j==0 || i==n-1 || j==n-1)
          continue;
           	temp = (*mat)[i][j];
	(*mat)[i][j] = 0.2 * ((*mat)[i][j] + (*mat)[i][j - 1] + (*mat)[i - 1][j] + (*mat)[i][j + 1] + (*mat)[i + 1][j]);
	diff += fabsf((*mat)[i][j] - temp);
        }
    }
    if (diff/n/n < TOL)
      done = 1; 
    cnt_iter ++;
  }

  if (done)
    printf("Solver converged after %d iterations\n", cnt_iter);
  else
    printf("Solver not converged after %d iterations\n", cnt_iter);
  
}

int main(int argc, char *argv[]) {
  int n;
  float **a;


 

  n = N;
  printf("Matrix size = %d \n", n);

  allocate_init_2Dmatrix(&a, n, n);

  // Initial operation time
  clock_t i_exec_t = clock();

  solver(&a, n, n);

  printf("Results:\n");
  for(int i=0 ; i<N;i++){
    for(int j=0;j<N;j++){
        printf("%.2lf\t",a[i][j]);
    }
    printf("\n");
  }

  // Final operation time
  clock_t f_exec_t = clock();
  float exec_time = (float)(f_exec_t - i_exec_t) / CLOCKS_PER_SEC;
  printf("Operations time: %f\n", exec_time);

  return 0;
}