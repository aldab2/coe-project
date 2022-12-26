#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <omp.h>
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
double accepted_err = 1E-7;
int N = 50000;
void term1(double(*x), double **b);
void compute_right(double(*x), double(*p), double ***a);
void compute_left(double(*x), double ***a);
// Generate a random double number with the maximum value of max
double rand_double(int max)
{
    return ((double)rand() / (double)(RAND_MAX)) * max;
}

void allocate_init_DD_2Dmatrix(double ***mat, double **b, int n, int m);

int main(int argc, char *argv[])
{
    double(*x) = malloc(sizeof(double[N]));
    double(**a); // = malloc(sizeof(double[N][N]));
    double(*err) = malloc(sizeof(double[N])), (*p) = malloc(sizeof(double[N]));
    double(*b); //= malloc(sizeof(double[N]));
    int i, j, k, iter;
    double start_time;
    double avgErr;
    double oldErr;
    printf("Starting...\n");
    allocate_init_DD_2Dmatrix(&a, &b, N, N);
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
        {

            x[i] = 0;
            p[i] = 0;
            err[i] = 100;
        }


    start_time = omp_get_wtime();

    // Outer Iterater
    for (iter = 1; iter <= MAXITER; iter++)
    {
        term1(x, &b);
        compute_right(x, p, &a);
        compute_left(x, &a);
        int errGreaterThanMax = 1;
        oldErr = avgErr;
        avgErr = 0;
        for (i = 0; i < N; i++)
        {
            err[i] = fabs(x[i] - p[i]);
            p[i] = x[i];
            avgErr += err[i];
        }
        avgErr /= N;

        if (fabs(oldErr-avgErr) < accepted_err)
        {
            printf("Converged after %d iterations ... Stopping\n", iter);
            break;
        }
        if (iter> 995)
        {
            printf(" change %.12lf\n",fabs(oldErr-avgErr));
        }

        // printf("Iteration:%d\t%0.4f\t%0.4f\t%0.4f\n", iter, x[0], x[1], x[2]);
    }
    double run_time = omp_get_wtime() - start_time;
    printf("time %f\n", run_time);
    {
        printf("%.10lf vs %.10lf (%d)\n", avgErr, accepted_err, avgErr >= accepted_err);
    }

    for (i = 0; i < N; i++)
    {
        // printf("\nSolution: x[%d]=%0.3f\n", i, x[i]);
    }
    free(a);
    free(x);
    free(b);
    free(p);
    return 0;
}

void term1(double(*x), double **b)
{
    for (int i = 0; i < N; i++)
    {
        x[i] = (*b)[i];
    }
}

void compute_right(double(*x), double(*p), double ***a)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            x[i] = x[i] - p[j] * (*a)[i][j];
        }
    }
}

void compute_left(double(*x), double ***a)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            {
                x[i] = x[i] - x[j] * (*a)[i][j];
            }
        }
        x[i] = x[i] / (*a)[i][i];
    }
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
                (*mat)[i][j] = n * MAX;
            }
            else
                (*mat)[i][j] = rand_double(MAX);
        }
        (*b)[i] = i + 1;
    }

    // printf("Initial Matrix:\n");
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            // printf("%.2lf\t",(*mat)[i][j] );
        }

        // printf("|%.2lf\n",(*b)[i]);
    }
}
