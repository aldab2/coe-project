#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <omp.h>

int MAX = 50;
int N= 4;

double rand_double(int max)
{
    return ((double)rand() / (double)(RAND_MAX)) * max;
}
void allocate_init_DD_2Dmatrix(double ***mat, double **b, int n, int m)
{
    int i, j;
    *mat = (double **)malloc(n * sizeof(double *));
    *b = (double *)malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
    {
        (*mat)[i] = (double *)malloc(m * sizeof(double));
        for (j = 0; j < m; j++)
        {
            if (i == j)
            {
                (*mat)[i][j] = (N)*MAX; // i + 1
            }
            else
                (*mat)[i][j] = rand_double(MAX);
        }
        (*b)[i] = i + 1;
    }

      printf("Initial Matrix:\n");
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            printf("%.2lf\t", (*mat)[i][j]);
        }
            printf("\n", (*b)[i]);
       // printf("|%.2lf\n", (*b)[i]);
    }

}


int main(){
double(**a); // = malloc(sizeof(double[N][N]));
    double(*err) = malloc(sizeof(double[N])), (*p) = malloc(sizeof(double[N]));
    double(*b);
        allocate_init_DD_2Dmatrix(&a, &b, N, N);

    int WIDTH= N;
    int HEIGHT = N;

    printf("Printing Daigs\n");
    for( int k = 1 ; k <= WIDTH + HEIGHT -3 ; k++ ) {
        for( int j = 0 ; j <= k ; j++ ) {
            int i = k - j;
            if( i < HEIGHT && j < WIDTH ) {
                printf(  "%.2lf\t",a[i][j] );
            }
        }
        printf("\n");
    }


    printf("Printing red elements\n");
    for(int i =0;i<N;i++){
        for(int j = i%2==0?0:1;j<(N/2)+2;j+=2){
           printf("%.2lf\t",a[i][j]);
        }
        printf ("\n");
    }

     printf("Printing black elements\n");
    for(int i =0;i<N;i++){
        for(int j = i%2==1?0:0;j<(N/2)+2;j+=2){
           printf("%.2lf\t",a[i][j]);
        }
        printf ("\n");
    }

}