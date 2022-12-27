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
#define MAXITER 20000
#define MAX 50
double accepted_err = 0.1E-50;
int N = 3;
void omp_term1(double(*x), double **b);
void omp_compute_right(double(*x), double(*p), double ***a);
void omp_compute_left(double(*x), double ***a);

void seq_term1(double(*x), double **b);
void seq_compute_right(double(*x), double(*p), double ***a);
void seq_compute_left(double(*x), double ***a);

double relative_residual(double(*x), double **b, double ***a);
// Generate a random double number with the maximum value of max
double rand_double(int max)
{
    return ((double)arc4random() / (double)(RAND_MAX)) * max;
}

void allocate_init_DD_2Dmatrix(double ***mat, double **b, int n, int m);

int main(int argc, char *argv[])
{
    double(*xp) = malloc(sizeof(double[N]));
    double(*xs) = malloc(sizeof(double[N]));
    double(**a); // = malloc(sizeof(double[N][N]));
    double(*err) = malloc(sizeof(double[N])), (*pp) = malloc(sizeof(double[N])), (*ps) = malloc(sizeof(double[N]));
    double(*b); //= malloc(sizeof(double[N]));
    int i, j, k, iter;
    double start_time, run_time;
    double avgErr;
    double oldErr;
    double res;
    double oldres;
    printf("Starting...\n");
    allocate_init_DD_2Dmatrix(&a, &b, N, N);
    for (j = 0; j < N; j++)
    {
        xp[j] = 0;
        xs[j] = 0;
        pp[j] = 0;
        ps[j] = 0;
        err[j] = 100;
    }

    start_time = omp_get_wtime();

    res = 100;
    oldres = 200;
    for (iter = 1; iter <= MAXITER; iter++)
    {
        seq_term1(xs, &b);
        seq_compute_right(xs, ps, &a);
        seq_compute_left(xs, &a);
        for (i = 0; i < N; i++)
        {
            ps[i] = xs[i];
        }
        oldres = res;
        res = relative_residual(xs, &b, &a);
        // printf("Iter %d residual %.8lf", iter, res);
        if (fabs(oldres - res) < accepted_err)
        {
            printf("iter(%d): %.12lf vs %.12lf\n", iter, fabs(oldres - res), accepted_err);

            printf("Sequential Converged after %d iterations ... Stopping\n", iter);
            break;
        }
        if (iter > MAXITER - 5)
            printf("iter(%d): %.12lf vs %.12lf\n", iter, fabs(oldres - res), accepted_err);

        // printf("Iteration:%d\t%0.4f\t%0.4f\t%0.4f\n", iter, x[0], x[1], x[2]);
    }
    run_time = omp_get_wtime() - start_time;
    printf("sequenatial time %f\n", run_time);
    {
        printf("%.10lf vs %.10lf \n", fabs(oldres - res), accepted_err);
    }

    start_time = omp_get_wtime();

    // Outer Iterater
    res = 100;
    oldres = 200;
    for (iter = 1; iter <= MAXITER; iter++)
    {
        omp_term1(xp, &b);
        omp_compute_right(xp, pp, &a);
        omp_compute_left(xp, &a);
        for (i = 0; i < N; i++)
        {
            err[i] = fabs(xp[i] - pp[i]);
            pp[i] = xp[i];
            avgErr += err[i];
        }
        oldres = res;
        res = relative_residual(xp, &b, &a);
        // printf("Iter %d residual %.8lf", iter, res);
        if (fabs(oldres - res) < accepted_err)
        {
            printf("iter(%d): %.12lf vs %.12lf\n", iter, fabs(oldres - res), accepted_err);

            printf("Parallel Converged after %d iterations ... Stopping\n", iter);
            break;
        }
        if (iter > MAXITER - 5)
            printf("iter(%d): %.12lf vs %.12lf\n", iter, fabs(oldres - res), accepted_err);

        // printf("Iteration:%d\t%0.4f\t%0.4f\t%0.4f\n", iter, x[0], x[1], x[2]);
    }
    run_time = omp_get_wtime() - start_time;
    printf("parallel time %f\n", run_time);
    {
        printf("%.10lf vs %.10lf \n", fabs(oldres - res), accepted_err);
    }

    for (i = 0; i < N; i++)
    {
        printf("\nSolution: x[%d]=%0.3f\n", i, xs[i]);
    }
    for (i = 0; i < N; i++)
    {
        printf("\nSolution: x[%d]=%0.3f\n", i, xp[i]);
    }
    free(a);
    free(xp);
    free(xs);
    free(b);
    free(pp);
    free(ps);
    return 0;
}

void omp_term1(double(*x), double **b)
{
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        x[i] = (*b)[i];
    }
}

void omp_compute_right(double(*x), double(*p), double ***a)
{
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            x[i] = x[i] - p[j] * (*a)[i][j];
        }
    }
}

void omp_compute_left(double(*x), double ***a)
{
    for (int i = 0; i < N; i++)
    {
        // #pragma omp parallel for
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

void seq_term1(double(*x), double **b)
{
    for (int i = 0; i < N; i++)
    {
        x[i] = (*b)[i];
    }
}

void seq_compute_right(double(*x), double(*p), double ***a)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            x[i] = x[i] - p[j] * (*a)[i][j];
        }
    }
}

void seq_compute_left(double(*x), double ***a)
{
    for (int i = 0; i < N; i++)
    {
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

        printf("|%.2lf\n", (*b)[i]);
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
