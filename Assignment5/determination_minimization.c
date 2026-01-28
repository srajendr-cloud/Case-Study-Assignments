#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Fortran LAPACK routines */
extern void dgetrf_(int*, int*, double*, int*, int*, int*);
extern void dgetri_(int*, double*, int*, int*, double*, int*, int*);

/* Build A(t) in column-major form */
void build_A(double *A, int n, double t)
{
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
        {
            if (i == j)
                A[i + j*n] = 1.0;
            else if (i < j)
                A[i + j*n] = (j - i + 1);
            else
                A[i + j*n] = t;
        }
}

/* Build dA/dt */
void build_dA(double *B, int n)
{
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            B[i + j*n] = (i > j) ? 1.0 : 0.0;
}

/* Compute d/dt det(A(t)) */
double det_derivative(int n, double t)
{
    double *A = malloc(n*n*sizeof(double));
    double *Ainv = malloc(n*n*sizeof(double));
    double *dA = malloc(n*n*sizeof(double));
    int *ipiv = malloc(n*sizeof(int));
    double *work = malloc(n*sizeof(double));
    int info;

    build_A(A, n, t);
    build_dA(dA, n);

    /* LU factorization */
    dgetrf_(&n, &n, A, &n, ipiv, &info);

    double det = 1.0;
    for (int i = 0; i < n; i++)
        det *= fabs(A[i + i*n]);

    /* Invert A */
    for (int i = 0; i < n*n; i++)
        Ainv[i] = A[i];

    dgetri_(&n, Ainv, &n, ipiv, work, &n, &info);

    /* Trace(A^{-1} dA) */
    double trace = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            trace += Ainv[i + j*n] * dA[j + i*n];

    free(A); free(Ainv); free(dA); free(ipiv); free(work);
    return det * trace;
}

int main()
{
    FILE *fp = fopen("determinant_results.dat", "w");

    for (int n = 3; n <= 9; n += 2)
    {
        double a = 0.0, b = 2.0, m;

        while (fabs(b - a) > 1e-8)
        {
            m = 0.5 * (a + b);
            if (det_derivative(n, a) * det_derivative(n, m) < 0)
                b = m;
            else
                a = m;
        }

        double t_star = 0.5 * (a + b);
        fprintf(fp, "%d %.8f\n", n, t_star);
    }

    fclose(fp);
    return 0;
}

