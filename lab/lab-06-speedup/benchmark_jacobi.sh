#!/usr/bin/env bash

for nthreads in `seq 12`; do
    export OMP_NUM_THREADS="$nthreads"
    for i in `seq 5`; do
        /usr/bin/time -f "$nthreads,$i,%e" ./jacobi 2 1 2 1 20000 1000 > /dev/null
    done
done
