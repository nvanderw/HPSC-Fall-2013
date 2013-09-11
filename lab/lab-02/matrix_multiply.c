#include <stdlib.h>
#include <string.h> // for memset
#include <stdio.h>

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

int main() {
}
