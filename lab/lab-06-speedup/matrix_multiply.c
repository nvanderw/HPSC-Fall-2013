#include <stdlib.h>
#include <string.h> // for memset
#include <stdio.h>
#include <errno.h>

#include <sys/resource.h>

#define ASSERT(x) if(!x) return 0;
#define TOLERANCE 0.00001
// #define TEST

/*
 * Multiplies a matrix a of size m x n and a matrix b of size n x p
 * into an allocation c of size m x p.
 */
void matrix_multiply(float *a, float *b, float *c,
                     size_t m, size_t n, size_t p) {
    memset(c, 0, sizeof(float) * m * p);
    #pragma omp parallel
    {
        #pragma omp for
        for(size_t i = 0; i < m; i++) // Rows of c
            for(size_t k = 0; k < n; k++)
                for(size_t j = 0; j < p; j++) // Columns of c
                    c[i * p + j] += a[i * n + k] * b[k * p + j];
    }
}

void matrix_print(FILE *out, float *matrix, size_t m, size_t n) {
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            fprintf(out, "%f ", matrix[i * n + j]);
        }
        fprintf(out, "\n");
    }
}

#ifdef TEST
/*
 * Compares two floating-point numbers within a specified tolerance. Returns
 * a number less than, equal to, or greater than zero.
 */
int float_cmp(float a, float b, float tol) {
    float delta = a - b;
    if(delta < (-tol)) return -1;
    if(delta > tol) return 1;
    else return 0;
}

/**
 * Uses float_cmp to compare an array element-by-element using a specified tolerance.
 * Returns a number less than, equal to, or greater than zero.
 */
int float_cmp_array(float *a, float *b, float tol, size_t nelem) {
    for(size_t i = 0; i < nelem; i++) {
        int res = float_cmp(a[i], b[i], tol);
        if(res != 0)
            return res;
    }
    return 0;
}

// A few unit tests

int test_unary() {
    size_t m = 1;
    size_t n = 1;
    size_t p = 1;

    float a = 1.0f;
    float b = 1.0f;
    float c;

    matrix_multiply(&a, &b, &c, m, n, p);
    ASSERT(c == 1.0f);
    return 1;
}

int test_131() {
    size_t m = 1;
    size_t n = 3;
    size_t p = 1;

    float a[3] = {1, 2, 3};
    float b[3] = {1, 2, 3};
    float c[1];

    matrix_multiply(a, b, c, m, n, p);
    ASSERT(!float_cmp(c[0], 14.0f, TOLERANCE));
    return 1;
}

int test_313() {
    size_t m = 3;
    size_t n = 1;
    size_t p = 3;

    float a[3] = {1, 2, 3};
    float b[3] = {1, 2, 3};
    float c[9];

    float reference[9] = {1.0f, 2.0f, 3.0f,
                          2.0f, 4.0f, 6.0f,
                          3.0f, 6.0f, 9.0f};

    matrix_multiply(a, b, c, m, n, p);

    ASSERT(!float_cmp_array(c, reference, TOLERANCE, 9));
    return 1;
}

/**
 * Check left and right multiplication by the identity
 */
int test_identity() {
    float i[4] = {1, 0,
                  0, 1};
    float a[4] = {1, 2,
                  3, 4};
    float b[4];

    matrix_multiply(i, a, b, 2, 2, 2);
    ASSERT(!float_cmp_array(b, a, TOLERANCE, 4));

    matrix_multiply(a, i, b, 2, 2, 2);
    ASSERT(!float_cmp_array(b, a, TOLERANCE, 4));

    return 1;
}

#endif

/**
 * Reads up to nelem elements from the file pointer input into the array of
 * floats pointed to by outarray.
 *
 * Returns the number of elements actually read, which may be less than nelem.
 */
size_t load_array(FILE *input, float *outarray, size_t nelem) {
    size_t count;
    for(count = 0; count < nelem; count++) {
        int result = fscanf(input, "%f ", &outarray[count]);
        if(result < 1) break;
    }

    return count;
}

void usage(char *name, FILE *out) {
    fprintf(out, "Usage: %s <matrix> <matrix>\n", name);
}

#ifndef TEST
int main(int argc, char **argv) {
    int status;

    if(argc < 3) {
        usage(argv[0], stderr);
        return 1;
    }

    // First, bump up our stack size.
    struct rlimit rlim;

    status = getrlimit(RLIMIT_STACK, &rlim);
    if(status == -1) {
        int err = errno;
        fprintf(stderr, "Could not get system resource limits. %s\n", strerror(err));
        return 1;
    }

    rlim.rlim_cur = 1024 * 1024 * 1024; // 1 GB
    status = setrlimit(RLIMIT_STACK, &rlim);
    if(status == -1) {
        int err = errno;
        fprintf(stderr, "Could not set system stack size. %s\n", strerror(err));
        return 1;
    }

    // Read in the file headers
    FILE *a_file, *b_file;
    int m_a, n_a, m_b, n_b;

    a_file = fopen(argv[1], "r");
    if(a_file == NULL) {
        int err = errno;
        fprintf(stderr, "Could not open file %s: %s\n", argv[1], strerror(err));
        return 1;
    }
    b_file = fopen(argv[2], "r");
    if(b_file == NULL) {
        int err = errno;
        fprintf(stderr, "Could not open file %s: %s\n", argv[2], strerror(err));
        return 1;
    }

    status = fscanf(a_file, "%d %d", &m_a, &n_a);
    if(status < 2) {
        fprintf(stderr, "Could not parse file header for %s\n", argv[1]);
        return 1;
    }

    status = fscanf(b_file, "%d %d", &m_b, &n_b);
    if(status < 2) {
        fprintf(stderr, "Could not parse file header for %s\n", argv[2]);
        return 1;
    }

    printf("%d %d\n", m_a, n_b);
    // Allocate the buffers
    float a_matrix[m_a][n_a];
    float b_matrix[m_b][n_b];

    status = load_array(a_file, &a_matrix[0][0], m_a * n_a);
    if(status < m_a * n_a) {
        fprintf(stderr, "Error: could not parse all entries of %s\n", argv[1]);
        return 1;
    }

    status = load_array(b_file, &b_matrix[0][0], m_b * n_b);
    if(status < m_b * n_b) {
        fprintf(stderr, "Error: could not parse all entries of %s\n", argv[2]);
        return 1;
    }

    float c_matrix[m_a][n_b];
    matrix_multiply(&a_matrix[0][0],
                    &b_matrix[0][0],
                    &c_matrix[0][0],
                    m_a, m_b, n_b);

    matrix_print(stdout, &c_matrix[0][0], m_a, n_b);

    return 0;
}
#else
int main() {
    if(test_unary() && test_131() && test_313() && test_identity()) {
        printf("All tests pass.\n");
    }
    else {
        printf("At least one test didn't pass.\n");
    }
}
#endif
