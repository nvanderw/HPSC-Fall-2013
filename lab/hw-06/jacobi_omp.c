#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>

#ifndef TOLERANCE
#define TOLERANCE 1e-5
#endif

#define ABS(x) (x < 0) ? -x : x

/*
 * Performs a Jacobi iteration on the given matrix.
 * Arguments:
 *      A: An n x n matrix, which will be modified in place
 *      n: Dimension of the matrix
 *      max_k: Maximum number of iterations
 *      eps: Upper bound on error
 *
 *  Returns: nonzero on convergence, zero otherwise
 */
int jacobi(double *A, size_t n, int max_k, double eps) {
    double B[n * n];
    memcpy(&B[0], A, sizeof(double) * n * n);

    double *new = A;
    double *old = &B[0];
    double *tmp;

    int converged;
    
    double delta[n][n];

    for(int k = 0; k < max_k; k++) {
        converged = 1;
        // Update the interior values
        #pragma omp parallel for
        for(size_t i = 1; i < n - 1; i++)
            for(size_t j = 1; j < n - 1; j++) {
                new[i * n + j] = (old[(i - 1) * n + j] + old[(i + 1) * n + j]
                                + old[i * n + (j - 1)] + old[i * n + (j + 1)])/4;

                if(ABS(new[i * n + j] - old[i * n + j]) > eps)
                    converged = 0;
            }

        if(converged) {
            if(new != A)
                memcpy(A, new, sizeof(double) * n * n);
            return 1;
        }

        // Swap the pointers
        tmp = old;
        old = new;
        new = tmp;
    }

    if(new != A)
        memcpy(A, new, sizeof(double) * n * n);
    return 0;
}

void init_jacobi(double *A, size_t n, double left, double right, double top, double bottom) {
    memset(A, 0, sizeof(double) * n * n);
    // Initialize the top
    for(size_t i = 0; i < n; i++) {
        A[i] = top;
    }

    for(size_t i = n * n - n; i < n * n; i++) {
        A[i] = bottom;
    }

    for(size_t i = 0; i < n * n; i += n) {
        A[i] = left;
    }
    
    for(size_t i = n - 1; i < n * n; i += n) {
        A[i] = right;
    }
}

void print_matrix(FILE *out, double *A, size_t m, size_t n) {
    for(size_t i = 0; i < m; i++) {
        for(size_t j = 0; j < n - 1; j++) {
            fprintf(out, "%f, ", A[i * n + j]);
        }
        fprintf(out, "%f\n", A[i * n + n - 1]);
    }
}

void usage(FILE *out, char *argv0) {
    fprintf(out, "Usage: %s <Top> <Bottom> <Left> <Right> <Iters> <N> <outfile.csv>\n", argv0);
}

int main(int argc, char **argv) {
    struct rlimit rlim;
    getrlimit(RLIMIT_STACK, &rlim);
    rlim.rlim_cur = 1024 * 1024 * 1024;
    setrlimit(RLIMIT_STACK, &rlim);

    if(argc != 8) {
        usage(stderr, argv[0]);
        exit(1);
    }

    double top, bottom, left, right;
    int max_k, n;

    top = atof(argv[1]);
    bottom = atof(argv[2]);
    left = atof(argv[3]);
    right = atof(argv[4]);
    max_k = atoi(argv[5]);
    n = atoi(argv[6]);
    const char *outpath = argv[7];


    double A[n][n];
    init_jacobi(&A[0][0], n, left, right, top, bottom);

    int convergence = jacobi(&A[0][0], n, max_k, TOLERANCE);
    // if(convergence)
    //     printf("Converged\n");
    // else
    //     printf("Did not converge\n");
    //
    
    FILE *outfile = fopen(outpath, "w");
    print_matrix(outfile, &A[0][0], n, n);
    fclose(outfile);
}
