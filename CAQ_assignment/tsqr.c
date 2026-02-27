#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


extern void dgeqrf_(int*, int*, double*, int*, double*, double*, int*, int*);
extern void dorgqr_(int*, int*, int*, double*, int*, double*, double*, int*, int*);

void local_qr(double* A, int m, int n, double* Q, double* R)
{
    int lda = m;
    int info;
    int lwork = -1;
    double wkopt;
    double* tau = (double*)malloc(n*sizeof(double));

    
    dgeqrf_(&m, &n, A, &lda, tau, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    double* work = (double*)malloc(lwork*sizeof(double));

    /* Computing QR */
    dgeqrf_(&m, &n, A, &lda, tau, work, &lwork, &info);

    /* Extracting R */
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            R[i*n+j] = (i<=j) ? A[i + j*m] : 0.0;

    /* Form Q */
    dorgqr_(&m, &n, &n, A, &lda, tau, work, &lwork, &info);

    for(int i=0;i<m*n;i++)
        Q[i] = A[i];

    free(work);
    free(tau);
}

int main()
{
    int m = 400;
    int n = 10;

    double* A = (double*)malloc(m*n*sizeof(double));
    srand(0);

    /* Generating random matrix */
    for(int i=0;i<m*n;i++)
        A[i] = (double)rand()/RAND_MAX;

    int mloc = m/4;

    /* Allocateing blocks */
    double *Q1 = malloc(mloc*n*sizeof(double));
    double *Q2 = malloc(mloc*n*sizeof(double));
    double *Q3 = malloc(mloc*n*sizeof(double));
    double *Q4 = malloc(mloc*n*sizeof(double));

    double *R1 = malloc(n*n*sizeof(double));
    double *R2 = malloc(n*n*sizeof(double));
    double *R3 = malloc(n*n*sizeof(double));
    double *R4 = malloc(n*n*sizeof(double));

    /* Level 0 local QR */
    local_qr(A, mloc, n, Q1, R1);
    local_qr(A+mloc*n, mloc, n, Q2, R2);
    local_qr(A+2*mloc*n, mloc, n, Q3, R3);
    local_qr(A+3*mloc*n, mloc, n, Q4, R4);

    printf("Local QR completed.\n");

    /* Now
       - Stack R1,R2 → QR
       - Stack R3,R4 → QR
       - Stack those → final QR
    */

    free(A);
    free(Q1); free(Q2); free(Q3); free(Q4);
    free(R1); free(R2); free(R3); free(R4);

    return 0;
}
