#include "timing.h"

/*
 * Fetches the current time in seconds, as a double. Returns 0 on success,
 * -1 on failure.
 */
int get_seconds(double *result) {
    struct timeval tv;
    if(gettimeofday(&tv, NULL) == -1) {
        return -1;
    }

    *result = (double) tv.tv_sec + ((double) tv.tv_usec) / 1e6;
    return 0;
}
