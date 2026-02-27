#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    int rank, size, n;
    double *vec = NULL, *local;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        FILE *f = fopen(argv[1], "r");
        fscanf(f, "%d", &n);
        vec = malloc(n * sizeof(double));
        for (int i = 0; i < n; i++) fscanf(f, "%lf", &vec[i]);
        fclose(f);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int local_n = n / size;
    local = malloc(local_n * sizeof(double));

    MPI_Scatter(vec, local_n, MPI_DOUBLE,
                local, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < local_n; i++) local[i] += 1.0;

    MPI_Gather(local, local_n, MPI_DOUBLE,
               vec, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        FILE *out = fopen("q1_output_collective.txt", "w");
        fprintf(out, "%d\n", n);
        for (int i = 0; i < n; i++) fprintf(out, "%lf ", vec[i]);
        fclose(out);
        free(vec);
    }

    free(local);
    MPI_Finalize();
    return 0;
}

