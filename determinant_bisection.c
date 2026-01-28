#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern void dgetrf_(int*, int*, double*, int*, int*, int*);

double determinant(int n, double *A) {
    int *ipiv = malloc(n*sizeof(int));
    int info;
    dgetrf_(&n, &n, A, &n, ipiv, &info);

    double det = 1.0;
    int swaps = 0;
    for (int i = 0; i < n; i++) {
        det *= A[i*n+i];
        if (ipiv[i] != i+1) swaps++;
    }
    if (swaps % 2) det = -det;

    free(ipiv);
    return fabs(det);
}

int main() {
    FILE *fp = fopen("determinant_results.dat", "w");

    for (int n = 3; n <= 9; n += 2) {
        double t = 1.0;

        double *A = malloc(n*n*sizeof(double));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                A[i*n+j] = (j >= i) ? (j - i + 1) : t;

        double det = determinant(n, A);
        fprintf(fp, "n=%d t*=%.4f det=%.6e\n", n, t, det);
        free(A);
    }

    fclose(fp);
    return 0;
}

