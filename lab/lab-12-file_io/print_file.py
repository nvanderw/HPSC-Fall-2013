#!/usr/bin/env python
import sys
import tables as tb

def read_h5(filename):
    with tb.openFile(filename, mode = "r") as h5:
    	X = h5.root.DATASET.read()
    return X

if __name__ == '__main__':


	b = read_h5(sys.argv[1])

	print b.shape

	print b