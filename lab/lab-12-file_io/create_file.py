#!/usr/bin/env python
import sys
import argparse
import numpy as np
import tables as tb

def get_args(argv):
    parser = argparse.ArgumentParser()  
    parser.add_argument('-m', help='e.g. -m 100')
    parser.add_argument('-n', help='e.g. -n 10')
    parser.add_argument('-f', help='e.g. -f matrix.10x10' )
    return parser.parse_args(argv)

def read_h5(filename):
    with tb.openFile(filename, mode = "r") as h5:
    	X = h5.root.DATASET.read()
    return X

def write_h5(filename, X):
	with tb.openFile(filename, mode = "w") as h5:
		h5.createArray(h5.root,'DATASET',X)

if __name__ == '__main__':

	args = get_args(sys.argv[1:])

	if not args.m or not args.n or not args.f:
		print 'please specify values -m, -n, -f'
		sys.exit()

	M = int(args.m)
	N = int(args.n)

	a = np.random.rand(M,N)

	write_h5(args.f, a)

	b = read_h5(args.f)

	diff = np.square(a-b).sum()
	assert diff < 1e-6

