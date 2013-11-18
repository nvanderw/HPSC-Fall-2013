#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include "functions.h"

double integrate(double (*f)(double), double x0, double h, int iters) {
    double sum = 0;

    for(int i = 0; i < iters; i++) {
        sum += f(x0) + 4*f(x0 + 0.5 * h) + f(x0 + h);
        x0 += h;
    }
    return h/6 * sum;
}

void usage(FILE *out, char *argv0) {
    fprintf(out, "Usage: %s <lower bound> <upper bound> <a|b|c|d> <iters>\n", argv0);
}

int main(int argc, char **argv) {
    int rank, nprocs;
    double sum;
    double (*f)(double); // Pointer to function to integrate over

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(argc < 5) {
        if(rank == 0) /* Master */
            usage(stderr, argv[0]);
        MPI_Finalize();
        return 1;
    }

    /* Read arguments */
    double lower = atof(argv[1]);
    double upper = atof(argv[2]);
    char   *fun  =      argv[3];
    int    iters = atoi(argv[4]);

    if(!strcmp(fun, "a"))
        f = func_a;
    else if(!strcmp(fun, "b"))
        f = func_b;
    else if(!strcmp(fun, "c"))
        f = func_c;
    else if(!strcmp(fun, "d"))
        f = func_d;
    else {
        if(rank == 0)
            printf("Invalid function name. Must match (a|b|c|d).\n");
        MPI_Finalize();
        return 1;
    }


    /* Figure out how big of a contiguous region of the function's domain
     * each process is responsible for. */
    double region_size = (upper - lower) / nprocs;

    /* What's the lower bound on my region? */
    double my_lower = region_size * rank + lower;

    /* Integrate over that region */
    double result = integrate(f, my_lower, region_size / iters, iters);

    MPI_Reduce((void*) &result, (void*) &sum, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);

    if(rank == 0)
        printf("area = %f\n", sum);

    MPI_Finalize();
}
