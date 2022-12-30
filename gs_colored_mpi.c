#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <mpi.h>

#define MAX_ITER 1000
#define MAX 100 // maximum value of the matrix element
#define TOL 0.001
#define N 8
#define THREAD_COUNT 4
// Generate a random float number with the maximum value of max
float rand_float(int max)
{
    return ((float)rand() / (float)(RAND_MAX)) * max;
}

// Allocate 2D matrix
void allocate_init_2Dmatrix(float ***mat1, int n, int m)
{
    int i, j;
    *mat1 = (float **)malloc(n * sizeof(float *));
    for (i = 0; i < n; i++)
    {
        (*mat1)[i] = (float *)malloc(m * sizeof(float));
        for (j = 0; j < m; j++)
        {
            (*mat1)[i][j] = rand_float(MAX);
        }
    }
}

// solver
void mpi_colored_solver(float ***mat, int n, int m, int p_cols,int rank,int size)
{
    float diff = 0, temp;
    int done = 0, cnt_iter = 0, i, j;
    MPI_Status status;

    while (!done && (cnt_iter < MAX_ITER))
    {

        diff = 0;

        for (int i = rank*p_cols; i < p_cols * (rank + 1); i++)
        {
            for (int j = i % 2 == 0 ? 2 : 1; j < n - 1; j += 2)
            {
                if (i == 0 || j == 0 || i == n - 1 || j == n - 1)
                    continue;
                temp = (*mat)[i][j];
                (*mat)[i][j] = 0.2 * ((*mat)[i][j] + (*mat)[i][j - 1] + (*mat)[i - 1][j] + (*mat)[i][j + 1] + (*mat)[i + 1][j]);
                diff += fabs((*mat)[i][j] - temp);
            }
            if(rank == 0 ){
                MPI_Send(&(mat[p_cols-1][0]),n,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
                MPI_Recv(&(mat[p_cols][0]),n,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&status);
            }
            if(rank == size -1){
                MPI_Recv(&(mat[p_cols][0]),n,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&status);
            }
        }

        for (int i = 1; i < n - 1; i++)
        {
            for (int j = i % 2 == 0 ? 1 : 2; j < n - 1; j += 2)
            {
                if (i == 0 || j == 0 || i == n - 1 || j == n - 1)
                    continue;
                temp = (*mat)[i][j];
                (*mat)[i][j] = 0.2 * ((*mat)[i][j] + (*mat)[i][j - 1] + (*mat)[i - 1][j] + (*mat)[i][j + 1] + (*mat)[i + 1][j]);
                diff += fabs((*mat)[i][j] - temp);
            }
        }
        if (diff / n / n < TOL)
            done = 1;
        cnt_iter++;
    }

    if (done)
        printf("Solver converged after %d iterations\n", cnt_iter);
    else
        printf("Solver not converged after %d iterations\n", cnt_iter);
}

int main(int argc, char *argv[])
{
    int n;
    // colored a
    float **ac;
    // sequential a
    float **a;
    n = N;
    printf("Matrix size = %d \n", n);
    allocate_init_2Dmatrix(&ac, n, n);
    int size, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int partition = N/size;

    // Colored Initial operation time
    clock_t colored_i_exec_t = clock();

    mpi_colored_solver(&ac, n, n, partition);

    printf("Colored Results:\n");
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%.2lf\t", ac[i][j]);
        }
        printf("\n");
    }
    // Final operation time
    clock_t colored_f_exec_t = clock();
    float colored_exec_time = (float)(colored_f_exec_t - colored_i_exec_t) / CLOCKS_PER_SEC;
    printf("Colored Operations time: %f\n", colored_exec_time);

    return 0;
}