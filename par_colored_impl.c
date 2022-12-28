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
double accepted_err = 0.1E-5;
int N = 300;
void term1(double(*x), double **b);
void compute_right(double(*x), double(*p), double ***a);
void compute_left(double(*x), double ***a);
double relative_residual(double(*x), double **b, double ***a);
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
    double res;
    double oldres;
    printf("Starting...\n");
    allocate_init_DD_2Dmatrix(&a, &b, N, N);
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
        {

            x[j] = 0;
            p[j] = 0;
            err[j] = 100;
        }

    start_time = omp_get_wtime();

    // Outer Iterater
    for (iter = 1; iter <= MAXITER; iter++)
    {
        // Setup (Calculate x0):
        // term1
        x[0] = b[0];
        // other terms
        for (j = 1; j < N; j++)
        {
            x[0] = x[0] - p[j] * (a)[0][j];
        }

        //devide by coofecient
         x[0] = x[0] / (a)[0][0];
        // update previous
        p[0] = x[0];

        for (i = 1; i < N; i++)
        {
            // term1
            x[i] = b[i];

// Black iteration in parallel
#pragma omp parallel for //reduction(+:x[i])
            for (j = i+1; j < N; j++)
            {
                
                x[i] = x[i] - p[j] * (a)[i][j];
            }
// Blue Iteration in parallel
#pragma omp parallel for// reduction(+:x[i])
            for (j = i - 1; j >= 0; j--)
            {

                x[i] = x[i] - x[j] * (a)[i][j];
            }
            //devide by coofecient
         x[i] = x[i] / (a)[i][i];
            // update previous
            p[i] = x[i];
        }
        oldres = res;
        res = relative_residual(x, &b, &a);
        // printf("Iter %d residual %.8lf", iter, res);
        if (fabs(oldres - res) < accepted_err)
        {
            printf("iter(%d): %.12lf vs %.12lf\n", iter, fabs(oldres - res), accepted_err);

            printf("Converged after %d iterations ... Stopping\n", iter);
            break;
        }
        if (iter > MAXITER - 5)
            printf("iter(%d): %.12lf vs %.12lf\n", iter, fabs(oldres - res), accepted_err);

        // printf("Iteration:%d\t%0.4f\t%0.4f\t%0.4f\n", iter, x[0], x[1], x[2]);
    }
    double run_time = omp_get_wtime() - start_time;
    printf("time %f\n", run_time);
    {
        printf("%.10lf vs %.10lf (%d)\n", res, accepted_err, avgErr >= accepted_err);
    }

     for (i = 0; i < 3; i++)
    {
        printf("\nSolution: x[%d]=%0.3f\n", i, x[i]);
    }
    free(a);
    free(x);
    free(b);
    free(p);
    return 0;
}

void term1(double(*x), double **b)
{
#pragma omp parallel for
    for (int i = 0; i < N; i++)
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
                (*mat)[i][j] = (i + 1) * MAX;
            }
            else
                (*mat)[i][j] = rand_double(MAX);
        }
        (*b)[i] = i + 1;
    }

    //printf("Initial Matrix:\n");
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            //printf("%.2lf\t", (*mat)[i][j]);
        }

        //printf("|%.2lf\n", (*b)[i]);
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
