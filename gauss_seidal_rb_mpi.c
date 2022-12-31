#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// #include <omp.h>
#include <math.h>
#include <mpi.h>

#define MAX_ITER 20
#define MAX 100 // maximum value of the matrix element
#define TOL 0.000001
#define THREAD_COUNT 16

// Generate a random double number with the maximum value of max
double rand_double(int max)
{
    return ((double)rand() / (double)(RAND_MAX)) * max;
}

// Allocate 2D matrix
void allocate_init_2Dmatrix(double ***mat1, int n, int m)
{
    int i, j;
    *mat1 = (double **)malloc(n * sizeof(double *));
    for (i = 0; i < n; i++)
    {
        (*mat1)[i] = (double *)malloc(m * sizeof(double));
        for (j = 0; j < m; j++)
        {
            (*mat1)[i][j] = rand_double(MAX);
        }
    }
}

// solver
void mpi_colored_solver(double ***mat, int n, int m, int p_rows, int rank, int size)
{
    double diff = 0, reduced_diff = 0, temp;
    int done = 0, cnt_iter = 0, i, j;
    MPI_Status status;

    while (!done && (cnt_iter < MAX_ITER))
    {

        diff = 0;
        // printf("rank %d Iter %d\n for int i=%d; i<%d i++\n", rank, cnt_iter, rank == 0 ? (rank * p_rows) + 1 : (rank * p_rows), rank == size - 1 ? (rank * p_rows) + p_rows - 1 : (rank * p_rows) + p_rows);
        for (int i = ((rank == 0) ? 1 : (rank * p_rows)); i < ((rank == size - 1) ? ((rank * p_rows) + p_rows - 1) : ((rank * p_rows) + p_rows)); i++)
        {
            for (int j = i % 2 == 0 ? 2 : 1; j < n - 1; j += 2)
            {
                // printf("rank:%d (%d,%d)\n", rank, i, j);
                if (i == 0 || j == 0 || i == n - 1 || j == n - 1)
                    continue;
                temp = (*mat)[i][j];
                (*mat)[i][j] = 0.2 * ((*mat)[i][j] + (*mat)[i][j - 1] + (*mat)[i - 1][j] + (*mat)[i][j + 1] + (*mat)[i + 1][j]);
                diff += fabs((*mat)[i][j] - temp);
            }
        }
        if (rank < size - 1)
        {
            if (cnt_iter == 0)
            {
                // printf("p:%d sending %.2lf to %d\n", rank, (*mat)[(rank * p_rows) + p_rows - 1][1], rank + 1);
            }
            MPI_Send((*mat)[(rank * p_rows) + p_rows - 1], n, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }
        if (rank > 0)
        {
            MPI_Recv((*mat)[((rank - 1) * p_rows) + p_rows - 1], n, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
            if (cnt_iter == 0)
            {
                // printf("p:%d receiving %.2lf from %d\n", rank, (*mat)[((rank - 1) * p_rows) + p_rows - 1][1], rank - 1);
            }
        }
        if (rank > 0)
        {
            if (cnt_iter == 0)
            {
                // printf("p:%d sending %.2lf to %d\n", rank, (*mat)[rank * p_rows][1], rank - 1);
            }
            MPI_Send(((*mat)[rank * p_rows]), n, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
        }
        if (rank < size - 1)
        {
            MPI_Recv(((*mat)[((rank + 1) * p_rows)]), n, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status);
            if (cnt_iter == 0)
            {
                // printf("p:%d receiving %.2lf from %d\n", rank, (*mat)[((rank + 1) * p_rows)][1], rank + 1);
            }
        }

        for (int i = ((rank == 0) ? 1 : (rank * p_rows)); i < ((rank == size - 1) ? ((rank * p_rows) + p_rows - 1) : ((rank * p_rows) + p_rows)); i++)
        {
            for (int j = i % 2 == 0 ? 1 : 2; j < n - 1; j += 2)
            {
                // printf("rank:%d (%d,%d)\n", rank, i, j);
                if (i == 0 || j == 0 || i == n - 1 || j == n - 1)
                    continue;
                temp = (*mat)[i][j];
                (*mat)[i][j] = 0.2 * ((*mat)[i][j] + (*mat)[i][j - 1] + (*mat)[i - 1][j] + (*mat)[i][j + 1] + (*mat)[i + 1][j]);
                diff += fabs((*mat)[i][j] - temp);
            }
        }
        if (rank < size - 1)
        {
            MPI_Send((*mat)[(rank * p_rows) + p_rows - 1], n, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }
        if (rank > 0)
        {
            MPI_Recv((*mat)[((rank - 1) * p_rows) + p_rows - 1], n, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
        }
        if (rank > 0)
        {
            MPI_Send((*mat)[rank * p_rows], n, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
        }
        if (rank < size - 1)
        {
            MPI_Recv((*mat)[((rank + 1) * p_rows)], n, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status);
        }

        MPI_Reduce(&diff, &reduced_diff, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            if (reduced_diff / n / n < TOL)
                done = 1;
        }
        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD); // *x = 4, 5, 6,7,8,
        cnt_iter++;
    }

    if (rank == 0)
    {
        if (done)
            printf("Solver converged after %d iterations\n", cnt_iter);
        else
            printf("Solver not converged after %d iterations\n", cnt_iter);
    }
    // collect to a
    if (rank > 0)
    {
        for (int i = rank * p_rows; i < (rank * p_rows) + p_rows; i++)
        {
            MPI_Send((*mat)[i], n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
    if (rank == 0)
    {
        for (int i = p_rows; i < n; i++)
        {
            int src = (i / p_rows) % size;
            MPI_Recv((*mat)[i], n, MPI_DOUBLE, src, 0, MPI_COMM_WORLD, &status);
        }
    }
}

int main(int argc, char *argv[])
{
    int n;
    double **a;
	double start, end;

    n = atoi(argv[1]);
    int print = atoi(argv[2]);

    printf("Matrix size = %d \n", n);
    allocate_init_2Dmatrix(&a, n, n);
    int size, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int p_rows = n / size;
    if (rank == 0)
    {
        printf("p_rows: %d\n", p_rows);
    }

    if (rank == 0)
	{
		start = MPI_Wtime();
	}
    mpi_colored_solver(&a, n, n, p_rows, rank, size);

    if (rank == 0)
    {
        end = MPI_Wtime();
        printf("Elapsed time is %f\n", end - start);
        clock_t colored_f_exec_t = clock();
        printf("Colored Results:\n");
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (print)
                    printf("%f, ", a[i][j]);
            }
            if (print)
                printf("\n");
        }
        // double colored_exec_time = (double)(colored_f_exec_t - colored_i_exec_t) / CLOCKS_PER_SEC;
        // printf("Colored Operations time: %f\n", colored_exec_time);
    }

    MPI_Finalize();
    return 0;
}