#!/usr/bin/env python
import sys
import argparse
import numpy as np

def get_args(argv):
    parser = argparse.ArgumentParser()  
    parser.add_argument('-a', help='e.g. -a matrix.A')
    parser.add_argument('-b', help='e.g. -n matrix.B')
    parser.add_argument('-c', help='e.g. -f matrix.C')
    return parser.parse_args(argv)

if __name__ == '__main__':

    args = get_args(sys.argv[1:])

    if not args.a or not args.b or not args.c:
        print 'please specify values -a, -b, and -c'
        sys.exit()

    # Read in values
    A = np.matrix(np.loadtxt(args.a, delimiter=' ', skiprows=1)) 
    print 'A', np.shape(A)
    B = np.matrix(np.loadtxt(args.b, delimiter=' ', skiprows=1)) 
    print 'B', np.shape(B)

    M = np.shape(A)[0]
    N = np.shape(B)[1]
    K = np.shape(A)[1]
    # Allocate solution
    C = np.zeros((M,N))
    print 'C', np.shape(C)

    # The right way...
    solution = np.dot(A,B)

    # Matrix Multiply... The homework way
    for m in xrange(M):
        for n in xrange(N):
            for k in xrange(K):
                C[m,n] += A[m,k]*B[k,n]

    assert (C == solution).all()

    #Write to file
    np.savetxt(args.c, C, fmt='%.5f', header="{0} {1}".format(M, N), comments='')
