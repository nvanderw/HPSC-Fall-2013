#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>

/*
 * Performs a Jacobi iteration on the given matrix.
 * Arguments:
 *      A: An n x n matrix, which will be modified in place
 *      n: Dimension of the matrix
 *      max_k: Maximum number of iterations
 *
 */
void jacobi(double *A, size_t n, int max_k) {
    double B[n * n];
    memcpy(&B[0], A, sizeof(double) * n * n);

    double *new = A;
    double *old = &B[0];
    double *tmp;

    for(int k = 0; k < max_k; k++) {
        // Update the interior values
        #pragma omp parallel for
        for(size_t i = 1; i < n - 1; i++)
            for(size_t j = 1; j < n - 1; j++)
                new[i * n + j] = (old[(i - 1) * n + j] + old[(i + 1) * n + j]
                                + old[i * n + (j - 1)] + old[i * n + (j + 1)])/4;

        // Swap the pointers
        tmp = old;
        old = new;
        new = tmp;
    }

    if(new != A)
        memcpy(A, new, sizeof(double) * n * n);
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

void print_matrix(double *A, size_t m, size_t n) {
    for(size_t i = 0; i < m; i++) {
        for(size_t j = 0; j < n; j++) {
            printf("%f ", A[i * n + j]);
        }
        printf("\n");
    }
}

void usage(FILE *out, char *argv0) {
    fprintf(out, "Usage: %s <Top> <Bottom> <Left> <Right> <Iters> <N>\n", argv0);
}

int main(int argc, char **argv) {
    struct rlimit rlim;
    getrlimit(RLIMIT_STACK, &rlim);
    rlim.rlim_cur = 1024 * 1024 * 1024;
    setrlimit(RLIMIT_STACK, &rlim);

    if(argc != 7) {
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

    double A[n][n];
    init_jacobi(&A[0][0], n, left, right, top, bottom);

    jacobi(&A[0][0], n, max_k);
    print_matrix(&A[0][0], n, n);
}
