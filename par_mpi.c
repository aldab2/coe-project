#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
// #include <omp.h>
#include <mpi.h>
/* Arranging given system of linear
   equations in diagonally dominant
   form:
   20x + y - 2z = 17
   3x + 20y -z = -18
   2x - 3y + 20z = 25
*/
// #define N 960
#define MAXITER 1000
#define MAX 50
double accepted_err =0;// 0.1E-5;
int N = 300;
void term1(double(*x), double **b);
void compute_right(double(*x), double(*p), double ***a);
void compute_left(double(*x), double ***a);
double relative_residual(double(*x), double **b, double ***a);
int malloc2ddouble(double ***array, int n, int m);
// Generate a random double number with the maximum value of max
double rand_double(int max)
{
    return ((double)rand() / (double)(RAND_MAX)) * max;
}

void allocate_init_DD_2Dmatrix(double ***mat, double **b, int n, int m);

int main(int argc, char *argv[])
{
    int size, rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double(*x) = malloc(sizeof(double[N]));
    double **a; // = malloc(sizeof(double[N][N]));
    double(*err) = malloc(sizeof(double[N])), (*p) = malloc(sizeof(double[N]));
    double(*b) = malloc(sizeof(double[N]));
    int i, j, k, iter;
    double start_time;
    double avgErr;
    double oldErr;
    double res;
    double oldres;
    MPI_Status status;

    for (j = 0; j < N; j++)
    {

        x[j] = 0;
        p[j] = 0;
        err[j] = 100;
    }
    malloc2ddouble(&a, N, N);

    if (rank == 0)
    {
        printf("Starting...\n");
        // init_DD_2Dmatrix(&a, &b, N, N);
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                if (i == j)
                {
                    a[i][j] = N * MAX;
                }
                else
                    a[i][j] = rand_double(MAX);
            }
            b[i] = i+1;
        }
        printf("Initial Matrix:\n");
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
               // printf("%.2lf\t", (a)[i][j]);
            }

           // printf("|%.2lf\n", (b)[i]);
        }

        
        // MPI_Bcast(&(a[0][0]), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (int i = 0; i < N; i++)
        {
            for (int r = 1; r < size; r++)
            {
                MPI_Send(&b[i], 1, MPI_DOUBLE, r, 0, MPI_COMM_WORLD);
            }
            for (int j = 0; j < N; j++)
            {
                for (int r = 1; r < size; r++)
                {
                    MPI_Send(&a[i][j], 1, MPI_DOUBLE, r, 0, MPI_COMM_WORLD);
                }
            }
        }
        printf("Broadcasting Done\n");
    }
    if (rank != 0)
        for (int i = 0; i < N; i++)
        {

            MPI_Recv(&b[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

            for (int j = 0; j < N; j++)
            {
                MPI_Recv(&a[i][j], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
            }
        }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 3)
    {
        printf("Initial Matrix from Proc %d:\n",rank);
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
               // printf("%.2lf\t", a[i][j]);
            }

           // printf("|%.2lf\n", (b)[i]);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank ==0){
        start_time = MPI_Wtime();
    }

    for (iter = 1; iter <= MAXITER; iter++)
    {
            MPI_Barrier(MPI_COMM_WORLD);



        // compute right
        for (int i = rank; i < N; i += size)
        {  
            x[i] = b[i];
            for (int j = i + 1; j < N; j++)
            {
                x[i] = x[i] - p[j] * a[i][j];
            }
            if (rank != 0)
            {
                //printf("iter %d p %d send x[%d]\n", iter, rank, i);
                MPI_Send(&x[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
            else
            {
                for (k = 1; k < size-1; k++)
                {
                    //printf("iter %d p %d receiving x[%d] from %d\n", iter, rank, (i*size)+k, k);

                    MPI_Recv(&x[(i*size)+k], 1, MPI_DOUBLE,k, 0, MPI_COMM_WORLD, &status);
                }
            }
        }

       MPI_Barrier(MPI_COMM_WORLD);

        if (rank == 0)
        {
            compute_left(x, &a);
            int errGreaterThanMax = 1;
            oldErr = avgErr;
            avgErr = 0;
            for (i = 0; i < N; i++)
            {
                err[i] = fabs(x[i] - p[i]);
                p[i] = x[i];
                avgErr += err[i];

                for (int r = 1; r < size; r++)
            {
                MPI_Send(&p[i], 1, MPI_DOUBLE, r, 0, MPI_COMM_WORLD);
            }
            }
             
            avgErr /= N;
            oldres = res;
            res = relative_residual(x, &b, &a);
            // printf("Iter %d residual %.8lf", iter, res);
            if (fabs(oldres - res) < accepted_err)
            {
                //printf("iter(%d): %.12lf vs %.12lf\n", iter, fabs(oldres - res), accepted_err);

               // printf("Converged after %d iterations ... Stopping\n", iter);
                break;
            }
            if (iter > MAXITER - 5)
                printf("iter(%d): %.12lf vs %.12lf\n", iter, fabs(oldres - res), accepted_err);
        }
         if(rank != 0){
            for (int i = 0; i < N; i++)
        {

            MPI_Recv(&p[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        }

        }
    }

    // printf("Iteration:%d\t%0.4f\t%0.4f\t%0.4f\n", iter, x[0], x[1], x[2]);

    if (rank == 0)
    {
        double run_time = MPI_Wtime() - start_time;
        printf("time %f\n", run_time);
        {
            printf("%.10lf vs %.10lf (%d)\n", res, accepted_err, avgErr >= accepted_err);
            for (i = 0; i < N; i++)
            {
                //printf("\nSolution: x[%d]=%0.3f\n", i, x[i]);
            }
        }
    }

    free(a);
    free(x);
    free(b);
    free(p);
    MPI_Finalize();
    return 0;
}

void term1(double(*x), double **b)
{
    for (int i = 0; i < N; i += 1)
    {
        x[i] = (*b)[i];
    }
}

void compute_right(double(*x), double(*p), double ***a)
{
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            x[j] = x[j] - p[j] * (*a)[i][j];
        }
    }
}

void compute_left(double(*x), double ***a)
{
    for (int i = 0; i < N; i++)
    {
#pragma omp parallel for
        for (int j = 0; j < i; j++)
        {
            // #pragma omp critical
            {
                x[i] = x[i] - x[j] * (*a)[i][j];
            }
        }
        x[i] = x[i] / (*a)[i][i];
    }
}

void init_DD_2Dmatrix(double ***mat, double **b, int n, int m)
{
    int i, j;
    //*mat = (double **)malloc(n * sizeof(double *));
    // b = (double *)malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
    {
        (*mat)[i] = (double *)malloc(m * sizeof(double));
        for (j = 0; j < m; j++)
        {
            if (i == j)
            {
                (*mat)[i][j] = (i + 1) * MAX;
            }
            else
                (*mat)[i][j] = rand_double(MAX);
        }
        //(*b)[i] = i + 1;
    }

    printf("Initial Matrix:\n");
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            printf("%.2lf\t", (*mat)[i][j]);
        }

        // printf("|%.2lf\n", (*b)[i]);
    }
}
double relative_residual(double(*x), double **b, double ***a)
{
    double(*ax) = malloc(sizeof(double[N]));
    double(*ax_minus_b) = malloc(sizeof(double[N]));

    double norm_ax_minus_b = 0, norm_b = 0;

    for (int i = 0; i < N; i++)
    {
        ax[i] = 0;
        for (int j = 0; j < N; j++)
        {
            ax[i] += (*a)[i][j] * x[j];
        }
    }
    for (int i = 0; i < N; i++)
    {
        ax_minus_b[i] = ax[i] - (*b)[i];
        norm_ax_minus_b += fabs(ax_minus_b[i]);
        norm_b += fabs((*b)[i]);
    }

    return norm_ax_minus_b / norm_b;
}

int malloc2ddouble(double ***array, int n, int m)
{

    /* allocate the n*m contiguous items */
    double *p = (double *)malloc(n * m * sizeof(double));
    if (!p)
        return -1;

    /* allocate the row pointers into the memory */
    (*array) = (double **)malloc(n * sizeof(double *));
    if (!(*array))
    {
        free(p);
        return -1;
    }

    /* set up the pointers into the contiguous memory */
    for (int i = 0; i < n; i++)
        (*array)[i] = &(p[i * m]);

    return 0;
}
