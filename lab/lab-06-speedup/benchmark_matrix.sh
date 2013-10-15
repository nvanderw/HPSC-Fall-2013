#!/usr/bin/env bash

for nthreads in `seq 12`; do
    export OMP_NUM_THREADS="$nthreads"
    for i in `seq 5`; do
        /usr/bin/time -f "$nthreads,$i,%e" ./matrix_multiply matrix.hdf5 matrix.hdf5 output.hdf5
    done
done
