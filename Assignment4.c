#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <openblas/cblas.h>

#define IDX(i,j,ld) ((i) + (j)*(ld))

double complex ONE  = 1.0 + 0.0*I;
double complex ZERO = 0.0 + 0.0*I;


double complex rand_complex() {
    return ((double)rand() / RAND_MAX)
         + ((double)rand() / RAND_MAX) * I;
}


void arnoldi_iteration() {
    int m = 100;
    int N = 10;

    double complex *A = malloc(m * m * sizeof(double complex));
    double complex *Q = calloc(m * (N + 1), sizeof(double complex));
    double complex *H = calloc((N + 1) * N, sizeof(double complex));
    double complex *v = malloc(m * sizeof(double complex));

    srand(0);

    
    for (int i = 0; i < m * m; i++)
        A[i] = rand_complex();

    
    for (int i = 0; i < m; i++)
        Q[IDX(i,0,m)] = rand_complex();

    
    double norm = cblas_dznrm2(m, Q, 1);
    cblas_zdscal(m, 1.0 / norm, Q, 1);

    
    for (int n = 0; n < N; n++) {
        
        cblas_zgemv(CblasColMajor, CblasNoTrans,
                    m, m,
                    &ONE,
                    A, m,
                    &Q[IDX(0,n,m)], 1,
                    &ZERO,
                    v, 1);

        for (int j = 0; j <= n; j++) {
            
            cblas_zdotc_sub(m,
                            &Q[IDX(0,j,m)], 1,
                            v, 1,
                            &H[IDX(j,n,N+1)]);

            
            double complex alpha = -H[IDX(j,n,N+1)];
            cblas_zaxpy(m, &alpha,
                        &Q[IDX(0,j,m)], 1,
                        v, 1);
        }

        
        H[IDX(n+1,n,N+1)] = cblas_dznrm2(m, v, 1);

        
        double scale = 1.0 / creal(H[IDX(n+1,n,N+1)]);
        cblas_zdscal(m, scale, v, 1);

        for (int i = 0; i < m; i++)
            Q[IDX(i,n+1,m)] = v[i];
    }

    
    double complex *AQ = malloc(m * N * sizeof(double complex));
    double complex *QH = malloc(m * N * sizeof(double complex));

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                m, N, m,
                &ONE,
                A, m,
                Q, m,
                &ZERO,
                AQ, m);

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                m, N, N+1,
                &ONE,
                Q, m,
                H, N+1,
                &ZERO,
                QH, m);

    double err = 0.0;
    for (int i = 0; i < m * N; i++)
        err += pow(cabs(AQ[i] - QH[i]), 2);

    printf("Arnoldi relation error ||AQ - QH|| = %.3e\n\n", sqrt(err));

    
    double complex *tmp = malloc((N+1) * sizeof(double complex));
    for (int i = 0; i < N+1; i++) {
        cblas_zgemv(CblasColMajor, CblasConjTrans,
                    m, N+1,
                    &ONE,
                    Q, m,
                    &Q[IDX(0,i,m)], 1,
                    &ZERO,
                    tmp, 1);

        printf("Q^H q_%d = ", i+1);
        for (int j = 0; j < N+1; j++)
            printf("(%+.2e) ", cabs(tmp[j]));
        printf("\n");
    }

    free(A); free(Q); free(H); free(v);
    free(AQ); free(QH); free(tmp);
}



void nilpotent_matrices() {
    printf("\nNilpotent matrix test:\n");

    for (int m = 2; m <= 50; m++) {
        double complex *A = calloc(m * m, sizeof(double complex));
        double complex *B = calloc(m * m, sizeof(double complex));

        
        for (int j = 0; j < m; j++)
            for (int i = 0; i < j; i++)
                A[IDX(i,j,m)] = rand_complex();

        for (int i = 0; i < m * m; i++)
            B[i] = A[i];

        int N;
        for (N = 1; N <= m; N++) {
            double norm = 0.0;
            for (int i = 0; i < m * m; i++)
                norm += cabs(B[i]);

            if (norm < 1e-12)
                break;

            cblas_ztrmm(CblasColMajor, CblasRight,
                        CblasUpper, CblasNoTrans, CblasNonUnit,
                        m, m,
                        &ONE,
                        A, m,
                        B, m);
        }

        printf("m = %d, nilpotency index N = %d\n", m, N);

        free(A); free(B);
    }
}



void alternative_orthogonalization() {
    int m = 5;
    int max_iter = 200;
    double h = 1e-3;
    double delta = 1e-5;

    double complex *C = malloc(m * m * sizeof(double complex));
    double complex *M = malloc(m * m * sizeof(double complex));
    double complex *R = malloc(m * m * sizeof(double complex));

    srand(1);
    for (int i = 0; i < m * m; i++)
        C[i] = rand_complex();

    for (int k = 0; k < max_iter; k++) {
        
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                    m, m, m,
                    &ONE,
                    C, m,
                    C, m,
                    &ZERO,
                    M, m);

        
        double trace = 0.0;
        for (int i = 0; i < m; i++)
            trace += creal(M[IDX(i,i,m)]);

        double err = fabs(trace - m);
        printf("Iteration %d: |Tr(CC^H) - m| = %.3e\n", k, err);

        if (err <= m * delta)
            break;

        
        for (int i = 0; i < m * m; i++)
            R[i] = -M[i];
        for (int i = 0; i < m; i++)
            R[IDX(i,i,m)] += 1.0;

        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    m, m, m,
                    &ONE,
                    R, m,
                    C, m,
                    &ZERO,
                    R, m);

        
        for (int i = 0; i < m * m; i++)
            C[i] += h * R[i];
    }

    free(C); free(M); free(R);
}


int main() {
    printf("===== PART 1: Arnoldi Iteration =====\n");
    arnoldi_iteration();

    printf("\n===== PART 2: Nilpotent Matrices =====\n");
    nilpotent_matrices();

    printf("\n===== PART 3: Alternative Orthogonalization =====\n");
    alternative_orthogonalization();

    return 0;
}

