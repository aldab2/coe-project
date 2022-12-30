#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <mpi.h>

// #define N 960
#define MAXITER 100

int N = 10;
int malloc2ddouble(double ***array, int n, int m);
int free2ddouble(double ***array);
int main(int argc, char *argv[])
{

    if (argc > 1)
        N = atoi(argv[1]);

    double **x, **a, **b, **x_total, **b_total;
    int i, j, k, iter;
    double dtime;
    int ssec, esec, susec, eusec;
    struct timeval tv;

    MPI_Init(NULL, NULL);
    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    malloc2ddouble(&x, N, N);
    malloc2ddouble(&a, N, N);
    malloc2ddouble(&b, N, N);

    malloc2ddouble(&x_total, N, N);
    malloc2ddouble(&b_total, N, N);

    // RANK0
    if (rank == 0)
    {

        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++)
            {
                x[i][j] = rand(); //(i+1) * 1.0; //(double)((rand() % 100) * 1.0);
                a[i][j] = rand(); //(i+1) * 1.0; //(double)((rand() % 100) * 1.0);
                b[i][j] = rand(); //(i+1) * 1.0; //(double)((rand() % 100) * 1.0);
            }
    }

    //  BROADCAST TO ALL RANKS
    MPI_Bcast(&(x[0][0]), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(a[0][0]), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(b[0][0]), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // calculate the partition for each process
    int partitition_size = N / size;

    if (rank == 0)
    {

        gettimeofday(&tv, NULL);
        ssec = tv.tv_sec;
        susec = tv.tv_usec;
    }

    //***************ADI code: main loop start****************************
    // Outer Iterater
    for (iter = 1; iter <= MAXITER; iter++)
    {
        //////ADI forward & backword sweep along rows//////
        // L1: Starts
        for (i = rank * partitition_size; i < partitition_size * (rank + 1); i++)
        {
            for (j = 1; j < N; j++)
            {
                x[i][j] = x[i][j] - x[i][j - 1] * a[i][j] / b[i][j - 1];
                b[i][j] = b[i][j] - a[i][j] * a[i][j] / b[i][j - 1];
            }
            x[i][N - 1] = x[i][N - 1] / b[i][N - 1];
        }
        // L2: Starts
        for (i = rank * partitition_size; i < partitition_size * (rank + 1); i++)
        {
            for (j = N - 2; j > 1; j--)
            {
                x[i][j] = (x[i][j] - a[i][j + 1] * x[i][j + 1]) / b[i][j];
            }
        }

        // Receive data from all processes into process 0
        MPI_Gather(&(x[rank * partitition_size][0]), partitition_size * N, MPI_DOUBLE, &(x_total[0][0]), partitition_size * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&(b[rank * partitition_size][0]), partitition_size * N, MPI_DOUBLE, &(b_total[0][0]), partitition_size * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        ////// ADI forward & backward sweep along columns//////
        // L3: Starts
        if (rank == 0)
        {

            for (j = 0; j < N; j++)
            {
                for (i = 1; i < N; i++)
                {
                    x_total[i][j] = x_total[i][j] - x_total[i - 1][j] * a[i][j] / b_total[i - 1][j];
                    x_total[i][j] = b_total[i][j] - a[i][j] * a[i][j] / b_total[i - 1][j];
                }
                x_total[N - 1][j] = x_total[N - 1][j] / b[N - 1][j];
            }
            // L4: Starts
            for (j = 0; j < N; j++)
            {
                for (i = N - 2; i > 1; i--)
                {
                    x_total[i][j] = (x_total[i][j] - a[i + 1][j] * x_total[i + 1][j]) / b_total[i][j];
                }
            }
        }
        //***************ADI code: main loop ends****************************
        MPI_Scatter(&(x_total[0][0]), partitition_size * N, MPI_DOUBLE, &(x[rank * partitition_size][0]), partitition_size * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(&(b_total[0][0]), partitition_size * N, MPI_DOUBLE, &(b[rank * partitition_size][0]), partitition_size * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    //***************ADI code: main loop end****************************
    if (rank == 0)
    {

        gettimeofday(&tv, NULL);
        esec = tv.tv_sec;
        eusec = tv.tv_usec;

        dtime = ((esec * 1.0) + ((eusec * 1.0) / 1000000.0)) - ((ssec * 1.0) + ((susec * 1.0) / 1000000.0));

        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                printf("%g ", x[i][j]);
            }
            printf("\n");
        }
        printf("%.3f\n", dtime);
    }

    free2ddouble(&x);
    free2ddouble(&a);
    free2ddouble(&b);
    free2ddouble(&x_total);
    free2ddouble(&b_total);

    MPI_Finalize();

    return 0;
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

int free2ddouble(double ***array)
{
    /* free the memory - the first element of the array is at the start */
    free(&((*array)[0][0]));

    /* free the pointers into the memory */
    free(*array);

    return 0;
}
