#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

/* BLAS Fortran routines */
extern double dznrm2_(int *n, double _Complex *x, int *incx);
extern void zdotc_(double _Complex *res, int *n,
                   double _Complex *x, int *incx,
                   double _Complex *y, int *incy);
extern void zaxpy_(int *n, double _Complex *alpha,
                   double _Complex *x, int *incx,
                   double _Complex *y, int *incy);

#define N 10
#define DIM 200

int main() {
    int n = DIM;
    int one = 1;

    double _Complex v[N][DIM];
    double _Complex u[N][DIM];

    /* Initialize random vectors */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < DIM; j++) {
            v[i][j] = (double)rand() / RAND_MAX +
                      (double)rand() / RAND_MAX * I;
        }



    }

    /* Modified Gram-Schmidt */
    

       for (int i = 0; i < N; i++) {

        double norm = dznrm2_(&n, v[i], &one);

        for (int k = 0; k < DIM; k++)
            u[i][k] = v[i][k] / norm;

        for (int j = i + 1; j < N; j++) {
            double _Complex beta;
            zdotc_(&beta, &n, u[i], &one, v[j], &one);

            double _Complex minus_beta = -beta;
            zaxpy_(&n, &minus_beta, u[i], &one, v[j], &one);
        }
    }

    /* Check norms */
    printf("Norms of orthonormal vectors:\n");
    for (int i = 0; i < N; i++) {
        double norm = dznrm2_(&n, u[i], &one);
        printf("||u[%d]|| = %.12f\n", i, norm);
    }

    /* Check orthogonality */
    printf("\nInner products:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double _Complex dot;
            zdotc_(&dot, &n, u[i], &one, u[j], &one);
            printf("<u[%d],u[%d]> = %.2e + %.2ei\n",
                   i, j, creal(dot), cimag(dot));
        }
        printf("\n");
    }

    return 0;
}
