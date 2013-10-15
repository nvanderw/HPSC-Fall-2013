#!/usr/bin/env bash

for nthreads in `seq 12`; do
    export OMP_NUM_THREADS="$nthreads"
    /usr/bin/time -f "$nthreads,%e" ./jacobi 2 1 2 1 20000 1000 > /dev/null
done
