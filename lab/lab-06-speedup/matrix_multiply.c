#include <stdlib.h>
#include <string.h> // for memset
#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <hdf5.h>

#include <sys/resource.h>

#ifndef BLOCKSIZE
#define BLOCKSIZE 100
#endif

#define ASSERT(x) if(!x) return 0;
#define TOLERANCE 0.00001
// #define TEST


void matrix_multiply_block(double *a, double *b, double *c,
                           size_t m, size_t n, size_t p,
                           size_t block_size) {
    memset(c, 0, sizeof(float) * m * p);
    #pragma omp parallel for
    for(size_t i = 0; i < m; i++) // Rows of c
        for(size_t k = 0; k < n; k++)
            for(size_t j = 0; j < p; j++) // Columns of c
                c[i * p + j] += a[i * n + k] * b[k * p + j];
}

/*
 * Multiplies a matrix a of size m x n and a matrix b of size n x p
 * into an allocation c of size m x p.
 */
void matrix_multiply(double *a, double *b, double *c,
                     size_t m, size_t n, size_t p) {
    matrix_multiply_block(a, b, c, m, n, p, 1);
}

void matrix_print(FILE *out, double *matrix, size_t m, size_t n) {
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            fprintf(out, "%f ", matrix[i * n + j]);
        }
        fprintf(out, "\n");
    }
}

#ifdef TEST
/*
 * Compares two doubleing-point numbers within a specified tolerance. Returns
 * a number less than, equal to, or greater than zero.
 */
int double_cmp(double a, double b, double tol) {
    double delta = a - b;
    if(delta < (-tol)) return -1;
    if(delta > tol) return 1;
    else return 0;
}

/**
 * Uses double_cmp to compare an array element-by-element using a specified tolerance.
 * Returns a number less than, equal to, or greater than zero.
 */
int double_cmp_array(double *a, double *b, double tol, size_t nelem) {
    for(size_t i = 0; i < nelem; i++) {
        int res = double_cmp(a[i], b[i], tol);
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

    double a = 1.0;
    double b = 1.0;
    double c;

    matrix_multiply(&a, &b, &c, m, n, p);
    ASSERT(c == 1.0f);
    return 1;
}

int test_131() {
    size_t m = 1;
    size_t n = 3;
    size_t p = 1;

    double a[3] = {1, 2, 3};
    double b[3] = {1, 2, 3};
    double c[1];

    matrix_multiply(a, b, c, m, n, p);
    ASSERT(!double_cmp(c[0], 14.0, TOLERANCE));
    return 1;
}

int test_313() {
    size_t m = 3;
    size_t n = 1;
    size_t p = 3;

    double a[3] = {1, 2, 3};
    double b[3] = {1, 2, 3};
    double c[9];

    double reference[9] = {1.0, 2.0, 3.0,
                           2.0, 4.0, 6.0,
                           3.0, 6.0, 9.0};

    matrix_multiply(a, b, c, m, n, p);

    ASSERT(!double_cmp_array(c, reference, TOLERANCE, 9));
    return 1;
}

/**
 * Check left and right multiplication by the identity
 */
int test_identity() {
    double i[4] = {1, 0,
                   0, 1};

    double a[4] =  {1, 2,
                    3, 4};
    double b[4];

    matrix_multiply(i, a, b, 2, 2, 2);
    ASSERT(!double_cmp_array(b, a, TOLERANCE, 4));

    matrix_multiply(a, i, b, 2, 2, 2);

    ASSERT(!double_cmp_array(b, a, TOLERANCE, 4));

    return 1;
}

int test_4x4_block_identity() {
    double i[16] = {1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1};

    double a[16] = {1,  2,  3,  4,
                    5,  6,  7,  8,
                    9,  10, 11, 12,
                    13, 14, 15, 16};
    double b[16];

    matrix_multiply_block(i, a, b, 4, 4, 4, 3);
    ASSERT(!double_cmp_array(b, a, TOLERANCE, 16));

    return 1;
}

double random_double() {
    int r = rand();
    return ((double) r) / ((double) RAND_MAX);
}

/* Uses a geometric RV to get a small number >= 1 */
size_t get_small_dimension(double p) {
    int dim = 1;
    while(random_double() > p)
        dim++;
    return dim;
}

size_t get_block_size(size_t m, size_t n) {
    size_t max = (m > n) ? m : n;
    return (int) (random_double() * (max + 1));
}

void randomize(double *a, size_t nelem) {
    for(size_t i = 0; i < nelem; i++)
        a[i] = (double) random_double();
}

int test_block_random() {
    int NUMTESTS = 10;
    for(int test = 0; test < NUMTESTS; test++) {
        // Generate random small m, n, p
        size_t m = get_small_dimension(0.1); // Expected value: 10
        size_t n = get_small_dimension(0.1);
        size_t p = get_small_dimension(0.1);

        double a[m * n];
        double b[n * p];
        double c1[m * p];
        double c2[m * p];

        randomize(a, m * n);
        randomize(b, n * p);

        printf("(m, n, p) = (%d, %d, %d)\n", m, n, p);
        matrix_multiply_block(a, b, c1, m, n, p, get_block_size(m, p));
        matrix_multiply_block(a, b, c2, m, n, p, get_block_size(m, p));
        ASSERT(!double_cmp_array(c1, c2, TOLERANCE, m * p));
    }

    return 1;
}

#endif

/**
 * Reads up to nelem elements from the file pointer input into the array of
 * doubles pointed to by outarray.
 *
 * Returns the number of elements actually read, which may be less than nelem.
 */
size_t load_array(FILE *input, double *outarray, size_t nelem) {
    size_t count;
    for(count = 0; count < nelem; count++) {
        int result = fscanf(input, "%f ", &outarray[count]);
        if(result < 1) break;
    }

    return count;
}

#ifndef TEST
void usage(FILE *out, char *argv0) {
    fprintf(out, "Usage: %s <A.h5> <B.h5> <C.h5>\n", argv0);
}

int main(int argc, char **argv) {
    struct rlimit rlim;
    getrlimit(RLIMIT_STACK, &rlim);
    rlim.rlim_cur = 1024 * 1024 * 1024;
    setrlimit(RLIMIT_STACK, &rlim);

    hid_t file_id, dataset_id, dataspace_id, status, property_id;
    hsize_t dims[2];

    if(argc < 4) {
        usage(stderr, argv[0]);
        exit(1);
    }

    file_id = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset_id = H5Dopen2(file_id, "x", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);

    status = H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

    hsize_t a_rows = dims[0], a_cols = dims[1];
    double a[a_rows][a_cols];
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, a);

    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);

    file_id = H5Fopen(argv[2], H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset_id = H5Dopen2(file_id, "x", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);

    status = H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

    hsize_t b_rows = dims[0], b_cols = dims[1];
    double b[b_rows][b_cols];
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, b);

    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);


    if(a_cols != b_rows) {
        fprintf(stderr, "Error: matrix dimension mismatch.\n");
        exit(1);
    }

    size_t m = a_rows;
    size_t n = a_cols;
    size_t p = b_cols;

    double c[m][p];
    matrix_multiply_block(&a[0][0], &b[0][0], &c[0][0], m, n, p, BLOCKSIZE);

    // Write out the output
    dims[0] = m;
    dims[1] = p;

    file_id = H5Fcreate(argv[3], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Screate_simple(2, dims, NULL);
    property_id = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_layout(property_id, H5D_CONTIGUOUS);

    dataset_id = H5Dcreate(file_id, "x", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, property_id, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &c[0][0]);

    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);
    status = H5Pclose(property_id);
}
#else
int main() {
    srand(time(NULL));
    if(test_block_random()) printf("All random tests passed.\n");
    else printf("A test failed.\n");
}
#endif
