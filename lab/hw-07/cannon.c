#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <hdf5.h>


#ifndef BLOCKSIZE
#define BLOCKSIZE 100
#endif

#define DATASET_NAME "DATASET"

void usage(FILE *out, const char *argv0) {
    fprintf(out, "Usage: %s <matrixA.h5> <matrixB.h5> <matrixC.h5>\n", argv0);
}

void matrix_multiply_block(double *a, double *b, double *c,
                           size_t m, size_t n, size_t p,
                           size_t block_size) {
    memset(c, 0, sizeof(double) * m * p);
    for(size_t ib = 0; ib < m; ib += block_size)
        for(size_t kb = 0; kb < n; kb += block_size)
            for(size_t jb = 0; jb < p; jb += block_size)
                for(size_t i = ib; (i < ib + block_size) && (i < m); i++)
                    for(size_t k = kb; (k < kb + block_size) && (k < n); k++)
                        for(size_t j = jb; (j < jb + block_size) && (j < p); j++)
                            c[i * p + j] += a[i * n + k] * b[k * p + j];
}

void matrix_add(const double *a, const double *b, double *c, size_t m, size_t n) {
    for(size_t i = 0; i < (m * n); i++)
        c[i] = a[i] + b[i];
}

void print_matrix(FILE *out, const double *A, size_t m, size_t n) {
    for(size_t i = 0; i < m; i++) {
        for(size_t j = 0; j < n; j++)
            fprintf(out, "%f ", A[i * n + j]);
        fprintf(out, "\n");
    }
}

int main(int argc, char **argv) {
    int world_rank, world_nprocs;
    herr_t status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_nprocs);

    if(argc != 4) {
        if(world_rank == 0)
            usage(stderr, argv[0]);

        MPI_Finalize();
        exit(-1);
    }

    // Load matrix A
    hid_t matrix_a_file = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t matrix_a_dataset = H5Dopen(matrix_a_file, DATASET_NAME, H5P_DEFAULT);
    hid_t matrix_a_space = H5Dget_space(matrix_a_dataset);

    hsize_t matrix_a_dims[2];
    hsize_t matrix_a_mdims[2];
    status = H5Sget_simple_extent_dims(matrix_a_space, matrix_a_dims, matrix_a_mdims);

    int a_rows = matrix_a_dims[0];
    int a_cols = matrix_a_dims[1];

    double matrix_a[a_rows][a_cols];
    status = H5Dread(matrix_a_dataset, H5T_NATIVE_DOUBLE,
                     H5S_ALL, H5S_ALL, H5P_DEFAULT, &matrix_a[0][0]);

    // Load matrix B
    hid_t matrix_b_file = H5Fopen(argv[2], H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t matrix_b_dataset = H5Dopen(matrix_b_file, DATASET_NAME, H5P_DEFAULT);
    hid_t matrix_b_space = H5Dget_space(matrix_b_dataset);

    hsize_t matrix_b_dims[2];
    hsize_t matrix_b_mdims[2];
    status = H5Sget_simple_extent_dims(matrix_b_space, matrix_b_dims, matrix_b_mdims);

    int b_rows = matrix_b_dims[0];
    int b_cols = matrix_b_dims[1];

    double matrix_b[b_rows][b_cols];
    status = H5Dread(matrix_b_dataset, H5T_NATIVE_DOUBLE,
                     H5S_ALL, H5S_ALL, H5P_DEFAULT, &matrix_b[0][0]);

    int n = a_rows;
    if((a_cols != n) || (b_rows != n) || (b_cols != n)) {
        if(world_rank == 0)
            fprintf(stderr, "Error: both matrices must be N x N\n");

        H5Sclose(matrix_b_space);
        H5Dclose(matrix_b_dataset);
        H5Fclose(matrix_b_file);

        H5Sclose(matrix_a_space);
        H5Dclose(matrix_a_dataset);
        H5Fclose(matrix_a_file);

        MPI_Finalize();
        exit(-1);
    }

    int procs_per_dim = (int) sqrt(world_nprocs);
    int local_n = n / procs_per_dim;

    // Set up our cartesian topology
    MPI_Comm comm_cart;
    int cart_rank, cart_nprocs;

    int cart_dims[2]    = {procs_per_dim, procs_per_dim};
    int cart_periods[2] = {1, 1}; // Wrap both directions

    MPI_Cart_create(MPI_COMM_WORLD, 2, cart_dims, cart_periods, 1, &comm_cart);
    MPI_Comm_rank(comm_cart, &cart_rank);
    MPI_Comm_size(comm_cart, &cart_nprocs);

    int cart_coords[2];
    MPI_Cart_coords(comm_cart, cart_rank, 2, cart_coords);

    int cart_row = cart_coords[0];
    int cart_col = cart_coords[1];

    double local_a[local_n][local_n];
    double local_b[local_n][local_n];

    // Initialize our matrices (FIXME)
    for(int i = 0; i < local_n; i++)
        for(int j = 0; j < local_n; j++)
            local_a[i][j] = matrix_a[local_n * cart_row + i][local_n * cart_col + j];

    for(int i = 0; i < local_n; i++)
        for(int j = 0; j < local_n; j++)
            local_b[i][j] = matrix_b[local_n * cart_row + i][local_n * cart_col + j];

    // Run the iteration
    for(int k = 0; k < procs_per_dim; k++) {
    }

    H5Sclose(matrix_b_space);
    H5Dclose(matrix_b_dataset);
    H5Fclose(matrix_b_file);

    H5Sclose(matrix_a_space);
    H5Dclose(matrix_a_dataset);
    H5Fclose(matrix_a_file);

    MPI_Finalize();
}
