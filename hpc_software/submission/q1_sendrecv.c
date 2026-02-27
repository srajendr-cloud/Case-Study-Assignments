#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    int rank, size, n;
    double *vec = NULL;
    double *local = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* Rank 0 reads the vector */
    if (rank == 0) {
        FILE *f = fopen(argv[1], "r");
        if (!f) {
            printf("Cannot open file\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fscanf(f, "%d", &n);
        vec = malloc(n * sizeof(double));
        for (int i = 0; i < n; i++)
            fscanf(f, "%lf", &vec[i]);
        fclose(f);
    }

    /* Broadcast vector length */
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int local_n = n / size;
    local = malloc(local_n * sizeof(double));

    /* Distribute using Send/Recv */
    if (rank == 0) {
        /* Copy own chunk */
        for (int i = 0; i < local_n; i++)
            local[i] = vec[i];

        /* Send chunks to other ranks */
        for (int p = 1; p < size; p++) {
            MPI_Send(vec + p * local_n, local_n,
                     MPI_DOUBLE, p, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(local, local_n, MPI_DOUBLE,
                 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    /* Add 1.0 locally */
    for (int i = 0; i < local_n; i++)
        local[i] += 1.0;

    /* Gather back using Send/Recv */
    if (rank == 0) {
        for (int i = 0; i < local_n; i++)
            vec[i] = local[i];

        for (int p = 1; p < size; p++) {
            MPI_Recv(vec + p * local_n, local_n,
                     MPI_DOUBLE, p, 1, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        }

        /* Write output */
        FILE *out = fopen("q1_output_sendrecv.txt", "w");
        fprintf(out, "%d\n", n);
        for (int i = 0; i < n; i++)
            fprintf(out, "%lf ", vec[i]);
        fprintf(out, "\n");
        fclose(out);

        free(vec);
    } else {
        MPI_Send(local, local_n, MPI_DOUBLE,
                 0, 1, MPI_COMM_WORLD);
    }

    free(local);
    MPI_Finalize();
    return 0;
}

