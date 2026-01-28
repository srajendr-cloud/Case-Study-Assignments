#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

/* Fortran LAPACK routines */
extern void zgetrf_(int*, int*, double complex*, int*, int*, int*);
extern void zgetrs_(char*, int*, int*, double complex*, int*, int*, double complex*, int*, int*);

/* BLAS */
extern void zgemv_(char*, int*, int*, double complex*, double complex*, int*,
                   double complex*, int*, double complex*, double complex*, int*);

#define MAX_ITER 50

double norm2(int n, double complex *x) {
    double sum = 0.0;
    for (int i = 0; i < n; i++)
        sum += creal(x[i])*creal(x[i]) + cimag(x[i])*cimag(x[i]);
    return sqrt(sum);
}

int main() {
    FILE *fp = fopen("iterations_vs_n.dat", "w");

    for (int n = 2; n <= 500; n++) {

        double complex *A = malloc(n*n*sizeof(double complex));
        double complex *A0 = malloc(n*n*sizeof(double complex));
        double complex *x_true = malloc(n*sizeof(double complex));
        double complex *b = malloc(n*sizeof(double complex));
        double complex *x = malloc(n*sizeof(double complex));
        double complex *r = malloc(n*sizeof(double complex));
        double complex *delta = malloc(n*sizeof(double complex));
        int *ipiv = malloc(n*sizeof(int));
        int info;

        for (int i = 0; i < n*n; i++)
            A[i] = drand48() + I*drand48();

        for (int i = 0; i < n; i++)
            x_true[i] = drand48() + I*drand48();

        char trans = 'N';
        double complex one = 1.0, zero = 0.0;
        int inc = 1;

        zgemv_(&trans, &n, &n, &one, A, &n, x_true, &inc, &zero, b, &inc);

        for (int i = 0; i < n*n; i++) A0[i] = A[i];

        zgetrf_(&n, &n, A, &n, ipiv, &info);

        for (int i = 0; i < n; i++) x[i] = b[i];
        zgetrs_(&trans, &n, &inc, A, &n, ipiv, x, &n, &info);

        int iter;
        double prev = 1e100;

        for (iter = 0; iter < MAX_ITER; iter++) {

            zgemv_(&trans, &n, &n, &one, A0, &n, x, &inc, &zero, r, &inc);
            for (int i = 0; i < n; i++) r[i] -= b[i];

            for (int i = 0; i < n; i++) delta[i] = r[i];
            zgetrs_(&trans, &n, &inc, A, &n, ipiv, delta, &n, &info);

            for (int i = 0; i < n; i++) x[i] -= delta[i];

            double err = norm2(n, delta);
            if (err >= prev) break;
            prev = err;
        }

        fprintf(fp, "%d %d\n", n, iter);

        free(A); free(A0); free(x_true); free(b);
        free(x); free(r); free(delta); free(ipiv);
    }

    fclose(fp);
    return 0;
}
