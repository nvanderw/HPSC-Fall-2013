#!/usr/bin/env bash

sizes=`seq 1000 1000 10000`
for n in $sizes; do
    python generate_matrix.py -m $n -n $n -f matrix."$n"x"$n".txt -t txt
    python generate_matrix.py -m $n -n $n -f matrix."$n"x"$n".h5 -t h5

    ./io_text matrix.*.txt
    ./io_hdf matrix.*.h5
done
