#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* LAPACK routine declarations */
extern void dgeqrf_(int*, int*, double*, int*, double*, double*, int*, int*);
extern void dorgqr_(int*, int*, int*, double*, int*, double*, double*, int*, int*);

void local_qr(double* A, int m, int n, double* Q, double* R)
{
    int lda = m;
    int info;
    int lwork = -1;
    double wkopt;
    double* tau = (double*)malloc(n*sizeof(double));

    /* Workspace query */
    dgeqrf_(&m, &n, A, &lda, tau, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    double* work = (double*)malloc(lwork*sizeof(double));

    /* Compute QR */
    dgeqrf_(&m, &n, A, &lda, tau, work, &lwork, &info);

    /* Extract R */
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

    for(int i=0;i<m*n;i++)
        A[i] = (double)rand()/RAND_MAX;

    int mloc = m/4;

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



    double *R12 = malloc(2*n*n*sizeof(double));
    double *R34 = malloc(2*n*n*sizeof(double));

    for(int j=0;j<n;j++)
        for(int i=0;i<n;i++){
            R12[i + j*(2*n)] = R1[i + j*n];
            R12[i+n + j*(2*n)] = R2[i + j*n];
        }

    for(int j=0;j<n;j++)
        for(int i=0;i<n;i++){
            R34[i + j*(2*n)] = R3[i + j*n];
            R34[i+n + j*(2*n)] = R4[i + j*n];
        }

    double *Q12 = malloc(2*n*n*sizeof(double));
    double *Q34 = malloc(2*n*n*sizeof(double));
    double *R12_new = malloc(n*n*sizeof(double));
    double *R34_new = malloc(n*n*sizeof(double));

    local_qr(R12, 2*n, n, Q12, R12_new);
    local_qr(R34, 2*n, n, Q34, R34_new);

    printf("Level 1 reduction done.\n");


    double *Rroot = malloc(2*n*n*sizeof(double));

    for(int j=0;j<n;j++)
        for(int i=0;i<n;i++){
            Rroot[i + j*(2*n)] = R12_new[i + j*n];
            Rroot[i+n + j*(2*n)] = R34_new[i + j*n];
        }

    double *Qroot = malloc(2*n*n*sizeof(double));
    double *Rfinal = malloc(n*n*sizeof(double));

    local_qr(Rroot, 2*n, n, Qroot, Rfinal);

    printf("Final reduction done.\n");


    free(R12); free(R34);
    free(Q12); free(Q34);
    free(R12_new); free(R34_new);
    free(Rroot); free(Qroot); free(Rfinal);

    free(A);
    free(Q1); free(Q2); free(Q3); free(Q4);
    free(R1); free(R2); free(R3); free(R4);

    return 0;
}
