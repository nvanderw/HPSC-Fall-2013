#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/resource.h>

#define RANK_MASTER 0
#define TOLERANCE 1e-5

#define ABS(x) (x < 0) ? -x : x

struct jacobi_parameters {
    // Values for each of the boundaries
    double top;
    double bottom;
    double left;
    double right;

    unsigned int max_k; // Number of iterations
    unsigned int n;     // Dimensions of Jacobi matrix, incl. boundary
};

void usage(FILE *out, char *argv0) {
    fprintf(out, "Usage: %s <Top> <Bottom> <Left> <Right> <Iters> <N>\n", argv0);
}

void update_boundary_values(double *A, // Matrix, including boundary (ghost) layer
                            size_t n,  // A is n x n, so this *includes* the ghost layer
                            const int *coords, // Array of length 2 in (row, col)
                            const int *threads, // Array of length 2
                            MPI_Comm cart_comm, // Cartesian communicator
                            const struct jacobi_parameters *params) {

    int row = coords[0];
    int col = coords[1];

    if(col == 0) { // We're on the left boundary
        double left = params->left;
        for(size_t i = 0; i < n; i++)
            A[i * n] = left;
    }

    else { // Not on left boundary; swap boundary info with left neighbor
        int us, them; // And after all we're only ordinary men
        MPI_Cart_shift(cart_comm, 1, -1, &us, &them);

        MPI_Request handle;
        // Start asynchronous receive of our neighbor's right boundary
        double their_right_boundary[n];
        double our_left_boundary[n];

        for(size_t i = 0; i < n; i++)
            our_left_boundary[i] = A[i * n + 1];

        MPI_Irecv(&their_right_boundary[0], n, MPI_DOUBLE,
                  them, MPI_ANY_TAG, cart_comm, &handle);
        MPI_Send(&our_left_boundary[0], n, MPI_DOUBLE, them,
                 0, cart_comm);
        MPI_Wait(&handle, MPI_STATUS_IGNORE);

        for(size_t i = 0; i < n; i++)
            A[i * n] = their_right_boundary[i];
    }

    if(col == threads[0] - 1) { // We're on the right
        double right = params->right;
        for(size_t i = 0; i < n; i++)
            A[i * n + n - 1] = right;
    }

    else { // Not on right boundary; swap info with right neighbor
        int us, them; // Haven't you heard it's a battle of words?
        MPI_Cart_shift(cart_comm, 1, 1, &us, &them);

        MPI_Request handle;
        double their_left_boundary[n];
        double our_right_boundary[n];

        for(size_t i = 0; i < n; i++)
            our_right_boundary[i] = A[i * n + n - 2];

        MPI_Irecv(&their_left_boundary[0], n, MPI_DOUBLE,
                  them, MPI_ANY_TAG, cart_comm, &handle);
        MPI_Send(&our_right_boundary[0], n, MPI_DOUBLE, them,
                 0, cart_comm);
        MPI_Wait(&handle, MPI_STATUS_IGNORE);

        for(size_t i = 0; i < n; i++)
            A[i * n + n - 1] = their_left_boundary[i];
    }

    if(row == 0) { // We're on the top boundary
        double top = params->top;
        for(size_t i = 0; i < n; i++)
            A[i] = top;
    }

    else { // Receive value from above
        int us, them; // With, without, and who'll deny it's what the
                      // fighting's all about?
        MPI_Cart_shift(cart_comm, 0, -1, &us, &them);

        MPI_Request handle;
        double their_bottom_boundary[n];
        double our_top_boundary[n];

        for(size_t i = 0; i < n; i++)
            our_top_boundary[i] = A[n + i];

        MPI_Irecv(&their_bottom_boundary[0], n, MPI_DOUBLE,
                  them, MPI_ANY_TAG, cart_comm, &handle);
        MPI_Send(&our_top_boundary[0], n, MPI_DOUBLE, them,
                 0, cart_comm);
        MPI_Wait(&handle, MPI_STATUS_IGNORE);

        for(size_t i = 0; i < n; i++)
            A[i] = their_bottom_boundary[i];
    }

    if(row == threads[1] - 1) { // We're on the bottom
        double bottom = params->bottom;
        for(size_t i = 0; i < n; i++)
            A[n * (n - 1) + i] = bottom;
    }

    else { // Receive value from below
        int us, them; // ...
        MPI_Cart_shift(cart_comm, 0, 1, &us, &them);

        MPI_Request handle;
        double their_top_boundary[n];
        double our_bottom_boundary[n];

        for(size_t i = 0; i < n; i++)
            our_bottom_boundary[i] = A[n * (n - 2) + i];

        MPI_Irecv(&their_top_boundary[0], n, MPI_DOUBLE,
                  them, MPI_ANY_TAG, cart_comm, &handle);
        MPI_Send(&our_bottom_boundary[0], n, MPI_DOUBLE, them,
                0, cart_comm);
        MPI_Wait(&handle, MPI_STATUS_IGNORE);

        for(size_t i = 0; i < n; i++)
            A[n * (n - 1) + i] = their_top_boundary[i];
    }
}

int update_interior_values(double *out, const double *A, size_t n, double eps) {
    int converged = 1;
    for(size_t i = 1; i < (n - 1); i++)
        for(size_t j = 1; j < (n - 1); j++) {
            out[i * n + j] = (A[(i - 1) * n + j] + A[(i + 1) * n + j]
                           + A[i * n + (j - 1)] + A[i * n + (j + 1)])/4;
            if(ABS(out[i * n + j] - A[i * n + j]) > eps)
                converged = 0;
        }
    
    return converged;
}

void print_matrix(const double *A, size_t m, size_t n) {
    for(size_t i = 0; i < m; i++) {
        for(size_t j = 0; j < n - 1; j++) {
            printf("%f, ", A[i * n + j]);
        }
        printf("%f\n", A[i * n + n - 1]);
    }
}

int main(int argc, char **argv) {
    // Increase the stack size so we can fit some big matrices
    struct rlimit rlim;
    getrlimit(RLIMIT_STACK, &rlim);
    rlim.rlim_cur = 1024 * 1024 * 1024;
    setrlimit(RLIMIT_STACK, &rlim);

    // Set up OpenMPI
    int rank, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    struct jacobi_parameters params;

    if(argc < 7) {
        if(rank == RANK_MASTER)
            usage(stderr, argv[0]);
        MPI_Finalize();
        return 1;
    }

    if(rank == RANK_MASTER) {
        // Parse command line arguments
        params.top    = atof(argv[1]);
        params.bottom = atof(argv[2]);
        params.left   = atof(argv[3]);
        params.right  = atof(argv[4]);
        params.max_k  = atoi(argv[5]);
        params.n      = atoi(argv[6]);
    }

    // Give everyone the parameters
    MPI_Bcast((void*) &params, sizeof(params), MPI_BYTE, RANK_MASTER, MPI_COMM_WORLD);

    // Assuming square number of processors, calculate number of processors in
    // each dimension
    unsigned int procs_per_dim = (unsigned int) sqrt((double) nprocs);
    unsigned int window_size = (params.n - 2) / procs_per_dim;
    if(rank == RANK_MASTER) {
        printf("Procs per dimension: %d\n", procs_per_dim);
        printf("Window size: %d\n", window_size);
    }

    int dims[2] = {procs_per_dim, procs_per_dim};
    int periods[2] = {0, 0}; // Don't wrap either dimension
    
    MPI_Comm cart_comm; // Our Cartesian communicator
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);
    int cart_rank;
    MPI_Comm_rank(cart_comm, &cart_rank);

    // Figure out our own coordinates
    int coords[2]; // (row, col)
    MPI_Cart_coords(cart_comm, cart_rank, 2, coords);

    // Allocate enough space for our window plus the boundary values on
    // either side
    int n = window_size + 2;
    double A[n][n];
    double B[n][n];

    memset(&A[0][0], 0, sizeof(A));

    for(int k = 0; k < params.max_k; k++) {
        update_boundary_values(&A[0][0], n, coords, dims, cart_comm, &params);
        memcpy(&B[0][0], &A[0][0], sizeof(A));

        int converged = update_interior_values(&A[0][0], &B[0][0], n, TOLERANCE);
        // Check if all processes have converged
        int all_converged;

        MPI_Allreduce(&converged, &all_converged, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        // if(all_converged) {
        //     if(rank == RANK_MASTER)
        //         printf("Terminated early!\n");
        //     break;
        // }
    }

    // for(int i = 0; i < procs_per_dim; i++) {
    //     for(int j = 0; j < procs_per_dim; j++) {
    //         if((coords[0] == i) && (coords[1] == j)) {
    //             printf("(row, col) = (%d, %d):\n", coords[0], coords[1]);
    //             
    //             print_matrix(&A[0][0], window_size + 2, window_size + 2);
    //         }
    //         MPI_Barrier(cart_comm);
    //     }
    // }


    // Now the iterations are done and we need to gather the interior values
    // from all of the nodes
    double interiors[n - 2][n - 2];

    for(int i = 1; i < (n - 1); i++)
        for(int j = 1; j < (n - 1); j++)
            interiors[i - 1][j - 1] = A[i][j];

    //               _________________
    //              /                 \
    //             | Whoo whoop whoop! |
    //              \ ________________/
    //              |/
    //       (V) (;,,;) (V)

    
    // Have the root gather all interior values
    int num_interiors = (n - 2) * (n - 2); // Number of interior values for each process
    double all_interiors[nprocs * num_interiors];

    MPI_Gather(&interiors[0][0], num_interiors, MPI_DOUBLE,
                &all_interiors[0], num_interiors, MPI_DOUBLE,
                RANK_MASTER, cart_comm);

    if(cart_rank == RANK_MASTER) {
        double result_matrix[params.n][params.n];

        for(int i = 0; i < params.n; i++) {
            result_matrix[0][i] = params.top;
            result_matrix[params.n - 1][i] = params.bottom;
            result_matrix[i][0] = params.left;
            result_matrix[i][params.n - 1] = params.right;
        }



        for(int M = 0; M < procs_per_dim; M++) { // Rows of processes
            for(int N = 0; N < procs_per_dim; N++) { // Columns of processes
                int their_coords[2] = {M, N};
                int rank;
                MPI_Cart_rank(cart_comm, &their_coords[0], &rank);
   
                double *their_interior = &all_interiors[rank * num_interiors];
                for(int i = 0; i < num_interiors; i++) {
                    int row = i / (n - 2);
                    int col = i % (n - 2);
                    result_matrix[1 + (n - 2) * M + row][1 + (n - 2) * N + col] = their_interior[i];
                }
            }
        }

        print_matrix(&result_matrix[0][0], params.n, params.n);
    }


    MPI_Finalize();
}

// There is no dark side of the moon really.
// Matter of fact it's all dark.
