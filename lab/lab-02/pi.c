#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>

#define SLICES 1e9
#define CLOCK_FREQ 2800 // MHz
#define F(X) 4/(1.0 + x * x)

double current_time() {
    struct timeval time;
    gettimeofday(&time, NULL);
    return time.tv_sec + 1e-6 * time.tv_usec;
}

int main(int argc, char **argv) {
    double sum = 0;
    double x;
    double dx = 1.0 / SLICES;

    double start_time = current_time();

    for(int i = 0; i < SLICES; i++) {
        sum += F(x);
        x   += dx;
    }

    double end_time = current_time();

    struct timeval stop_time;
    gettimeofday(&stop_time, NULL);

    printf("Pi: %f\n", sum * dx);
    double time = end_time - start_time;
    printf("time = %f\n", time);
    // Assuming 1 FLOP per loop iteration
    double mflops = SLICES / time / 1e6;
    printf("MFLOPS = %f\n", mflops);
    double latency = CLOCK_FREQ / mflops;
    printf("latency = %f\n", latency);
}
