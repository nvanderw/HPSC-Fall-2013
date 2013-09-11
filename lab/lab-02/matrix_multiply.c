#include <stdlib.h>
#include <string.h> // for memset
#include <stdio.h>
#include <errno.h>

#define ASSERT(x) if(!x) return 0;

/*
 * Multiplies a matrix a of size m x n and a matrix b of size n x p
 * into an allocation c of size m x p.
 */
void matrix_multiply(float *a, float *b, float *c,
                     size_t m, size_t n, size_t p) {
    // Zero out the output matrix
    memset((void *) c, 0, sizeof(float) * m * p);

    for(size_t j = 0; j < m; j++) // Rows of c
        for(size_t i = 0; i < p; i++) // Columns of c
            for(size_t k = 0; k < n; k++)
                c[j * p + i] += a[i * n + k] * b[k * p + j];
}

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
    ASSERT(!float_cmp(c[0], 14.0f, 0.0001));
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

    ASSERT(!float_cmp_array(c, reference, 0.0001, 9));
    return 1;
}

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

int main(int argc, char **argv) {
    if(argc < 3) {
        usage(argv[0], stderr);
        return 1;
    }

    FILE *a_file, *b_file;
    int m_a, n_a, m_b, n_b;

    a_file = fopen(argv[1], "r");
    if(a_file == NULL) {
        int err = errno;
        fprintf(stderr, "Could not open file %s: %s\n", argv[1], strerror(errno));
        return 1;
    }
    b_file = fopen(argv[2], "r");
    if(b_file == NULL) {
        int err = errno;
        fprintf(stderr, "Could not open file %s: %s\n", argv[2], strerror(errno));
        return 1;
    }

    int status;

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

    if(n_a != m_b) {
        fprintf(stderr, "Matrix dimensions do not match.\n");
        return 1;
    }

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

    for(size_t i = 0; i < m_a; i++) {
        for(size_t j = 0; j < n_b; j++) {
            printf("%f ", c_matrix[i][j]);
        }
        printf("\n");
    }

    return 0;
}
