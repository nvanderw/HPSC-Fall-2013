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
    parser.add_argument('-t', help='e.g. -t h5 or -t txt' )
    return parser.parse_args(argv)

if __name__ == '__main__':

    args = get_args(sys.argv[1:])

    if not args.m or not args.n or not args.f and not args.t:
        print 'please specify values -m, -n, -f, and -t'
        sys.exit()

    M = int(args.m)
    N = int(args.n)

    a = np.random.rand(M,N)

    # Save as txt
    if args.t == 'txt':
        filename = args.f
        with open(filename,'w') as outfile:
            outfile.write(str(M)+' '+str(N)+'\n')
            np.savetxt(outfile, a, fmt='%.5f')

    # Save as HDF5
    if args.t == 'h5':
        filename = args.f
        # h5file = tb.openFile(filename, mode = "w")
        # h5file.createArray(h5file.root,'x',a)
        # h5file.close()
        with tb.openFile(filename, mode = "w") as h5:
            h5.createArray(h5.root,'x',a)
















