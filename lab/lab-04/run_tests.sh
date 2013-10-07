#!/bin/bash

echo "N, your_time, dot"
ns="10 100 200 400 800 1600 3200"
for n in $ns; do
    python ../lab-03/generate_matrix.py -m $n -n $n -t h5 -f matrixA.h5
    python ../lab-03/generate_matrix.py -m $n -n $n -t h5 -f matrixB.h5
    results="`python compare.py -a matrixA.h5 -b matrixB.h5`"
    dot=`echo "$results" | grep dot | awk '{print $2}'`
    your_time=`echo "$results" | grep your | awk '{print $3}'`
    echo "$n, $your_time, $dot"
done
