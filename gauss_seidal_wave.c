#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#define MAX_ITER 1000
#define MAX 100 // maximum value of the matrix element
#define TOL 0.000000

// Generate a random double number with the maximum value of max
double rand_double(int max)
{
    return ((double)rand() / (double)(RAND_MAX)) * max;
}

// Allocate 2D matrix
void allocate_init_2Dmatrix(double ***mat, int n, int m)
{
    int i, j;
    *mat = (double **)malloc(n * sizeof(double *));
    for (i = 0; i < n; i++)
    {
        (*mat)[i] = (double *)malloc(m * sizeof(double));
        for (j = 0; j < m; j++)
            (*mat)[i][j] = rand_double(MAX);
    }
}

// solver
void solver(double ***mat, int n, int m)
{
    double diff = 0, temp;
    int done = 0, cnt_iter = 0, i, j, diag;
    int diagnolas = (2 * n) - 1;
    while (!done && (cnt_iter < MAX_ITER))
    {
        diff = 0;
        // loop over diagnolas sequantioaly
        for (diag = 1; diag < diagnolas - 1; diag++)
        {
// printf("\nN diag - (n - 1): %d, (diagnolas)-diag: %d \n", (diag - (n - 1)), (diagnolas)-diag);
// can be paralleized, reduce diff
//((diagnolas - 1) - diag) + diag - (n - 1): move to the base of start and then add
            #pragma omp parallel for num_threads(4) private(temp, i, j) reduction(+:diff)
                for (int d = diag > (n - 1) ? (diag - (n - 1)) : 0; d <= ((diag > (n - 1)) ? ((diagnolas - 1) - diag) + diag - (n - 1) : diag); d++) // d=0,d=1.  --  d=0,d=1,d=2  -- 0,1
                {
                    i = d;
                    j = abs(diag - d);
                    // printf("i: %d, j: %d - ", i, j);
                    if (i == 0 || j == 0 || i == n - 1 || j == n - 1)
                    {
                        continue;
                    }
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
        printf("\nSolver converged after %d iterations\n", cnt_iter);
    else
        printf("\nSolver not converged after %d iterations\n", cnt_iter);
}

int main(int argc, char *argv[])
{
    int n;
    double **a;

    if (argc < 2)
    {
        printf("Call this program with two parameters: matrix_size communication \n");
        printf("\t matrix_size: Add 2 to a power of 2 (e.g. : 18, 1026)\n");

        exit(1);
    }

    n = atoi(argv[1]);
    printf("Matrix size = %d \n", n);

    allocate_init_2Dmatrix(&a, n, n);

    // Initial operation time
    clock_t i_exec_t = clock();

    solver(&a, n, n);

    // Final operation time
    clock_t f_exec_t = clock();
    double exec_time = (double)(f_exec_t - i_exec_t) / CLOCKS_PER_SEC;
    printf("Operations time: %f\n", exec_time);

    return 0;
}