#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
// #include <omp.h>
/* Arranging given system of linear
   equations in diagonally dominant
   form:
   20x + y - 2z = 17
   3x + 20y -z = -18
   2x - 3y + 20z = 25
*/
//#define N 960
#define MAXITER 20
double max_err = 0.000000000000000000000000001;
int N = 3;
void save_to_file(double (*x)[N], char filename[]);
int main(int argc, char *argv[])
{
    double(*x) = malloc(sizeof(double[N])), (*a)[N] = malloc(sizeof(double[N][N]));
    double(*err) = malloc(sizeof(double[N])), (*c) = malloc(sizeof(double[N]));
    double(*b) = malloc(sizeof(double[N]));
    int i, j, k, iter;
    double dtime;
    int ssec, esec, susec, eusec;
    struct timeval tv;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
        {

            x[i] = 0;
            c[i] = 0;
            err[i] = 0;
        }
    // init a
    //    20x + y - 2z = 17
    //    3x + 20y -z = -18
    //    2x - 3y + 20z = 25
    a[0][0] = 20;
    a[0][1] = 1;
    a[0][2] = -2;
    a[1][0] = 3;
    a[1][1] = 20;
    a[1][2] = -1;
    a[2][0] = 2;
    a[2][1] = -3;
    a[2][2] = 20;
    // init b
    b[0] = 17;
    b[1] = -18;
    b[2] = 25;

    gettimeofday(&tv, NULL);
    ssec = tv.tv_sec;
    susec = tv.tv_usec;

    // Outer Iterater
    for (iter = 1; iter <= MAXITER; iter++)
    {
        for (i = 0; i < N; i++)
        {
            // x[i] = b
            //x[0] = b[0] - x[1]* a[0][1] - x[2]*a[0][2]
            //x[1] = b[1] - x[0]*a[1][0] - x[2]*a[1][2]
            //x[2] = b[1] - x[0]*a[1][0] - x[2]*a[1][2]
            x[i] = b[i];
            for (j = 0; j < N; j++)
            {
                if(i != j){
                   x[i] = x[i] - (x[j] * a[i][j]); 
                }
            }
            x[i] = x[i] / a[i][i];
            err[i] = fabs(x[i] - c[i]);
            c[i] = x[i];
        }
        int errGreaterThanMax = 0;

        for (int k = 0; k < N; k++)
        {

            if (err[k] >= max_err)
            {
                errGreaterThanMax = 1;
            }
        }

        if (errGreaterThanMax == 0)
        {
            printf("Reached max error... Stopping\n");
            break;
        }

        printf("Iteration:%d\t%0.4f\t%0.4f\t%0.4f\n", iter, x[0], x[1], x[2]);
    }
    gettimeofday(&tv, NULL);
    esec = tv.tv_sec;
    eusec = tv.tv_usec;
    dtime = ((esec * 1.0) + ((eusec * 1.0))) - ((ssec * 1.0) + ((susec * 1.0)));
    printf("time %f\n", dtime);

    for (i = 0; i < N; i++)
    {
        printf("\nSolution: x[%d]=%0.3f\n", i, x[i]);
    }
    free(a);
    free(x);
    free(b);
    return 0;
}
