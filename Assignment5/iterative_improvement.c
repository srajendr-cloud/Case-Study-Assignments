#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

/* Fortran LAPACK routines from MKL */
extern void zgetrf_(int*, int*, double _Complex*, int*, int*, int*);
extern void zgetrs_(char*, int*, int*, double _Complex*, int*, int*, double _Complex*, int*, int*);

#define MAX_ITERS 50

double norm2(int n, double _Complex *v)
{
    double s = 0.0;
    for (int i = 0; i < n; i++)
        s += creal(v[i])*creal(v[i]) + cimag(v[i])*cimag(v[i]);
    return sqrt(s);
}

int main()
{
    FILE *fp = fopen("iterations_vs_n.dat", "w");

    for (int n = 2; n <= 500; n++)
    {
        double _Complex *A  = malloc(n*n*sizeof(*A));
        double _Complex *A0 = malloc(n*n*sizeof(*A0));
        double _Complex *x  = malloc(n*sizeof(*x));
        double _Complex *x0 = malloc(n*sizeof(*x0));
        double _Complex *b  = malloc(n*sizeof(*b));
        double _Complex *r  = malloc(n*sizeof(*r));
        double _Complex *d  = malloc(n*sizeof(*d));
        int *ipiv = malloc(n*sizeof(int));

        /* Build A (column-major) */
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
            {
                A[i + j*n]  = ((double)rand()/RAND_MAX) + I*((double)rand()/RAND_MAX);
                A0[i + j*n] = A[i + j*n];
            }

        for (int i = 0; i < n; i++)
            x0[i] = ((double)rand()/RAND_MAX) + I*((double)rand()/RAND_MAX);

        /* b = A x_true */
        for (int i = 0; i < n; i++)
        {
            b[i] = 0.0;
            for (int j = 0; j < n; j++)
                b[i] += A0[i + j*n] * x0[j];
        }

        int info;
        zgetrf_(&n, &n, A, &n, ipiv, &info);

        for (int i = 0; i < n; i++)
            x[i] = b[i];

        char trans = 'N';
        int nrhs = 1;
        zgetrs_(&trans, &n, &nrhs, A, &n, ipiv, x, &n, &info);

        double prev_err = 1e100;
        int iter;

        for (iter = 0; iter < MAX_ITERS; iter++)
        {
            for (int i = 0; i < n; i++)
            {
                r[i] = -b[i];
                for (int j = 0; j < n; j++)
                    r[i] += A0[i + j*n] * x[j];
            }

            for (int i = 0; i < n; i++)
                d[i] = r[i];

            zgetrs_(&trans, &n, &nrhs, A, &n, ipiv, d, &n, &info);

            for (int i = 0; i < n; i++)
                x[i] -= d[i];

            for (int i = 0; i < n; i++)
                d[i] = x[i] - x0[i];

            double err = norm2(n, d) / norm2(n, x0);
            if (err >= prev_err) break;
            prev_err = err;
        }

        fprintf(fp, "%d %d\n", n, iter);

        free(A); free(A0); free(x); free(x0);
        free(b); free(r); free(d); free(ipiv);
    }

    fclose(fp);
    return 0;
}

