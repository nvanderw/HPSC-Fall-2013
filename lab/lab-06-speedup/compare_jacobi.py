#!/usr/bin/env python
import numpy as np
import argparse
import sys
import subprocess
import time

def get_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-t')
    parser.add_argument('-b')
    parser.add_argument('-r')
    parser.add_argument('-l')
    parser.add_argument('-i')
    parser.add_argument('-n')
    return parser.parse_args(argv)

def local_jacobi(args):
    nx = int(args.n)
    ny = int(args.n)
    A = np.zeros((nx, ny))

    A[0,:] = int(args.t)   # Top
    A[nx-1,:] = int(args.b)  # Bottom
    A[:,0] = int(args.l)     # Left
    A[:,ny-1] = int(args.r)    # Right

    A_new = A.copy()
    error = A_new - A

    k = 1 
    epsilon = 1e-5
    max_iterations = int(args.i)
    while k <= max_iterations:
        for i in np.arange(1., nx-1):
            for j in np.arange(1.,  ny-1):
                A_new[i,j] = (A[i,j+1] + A[i,j-1] + A[i-1,j] + A[i+1,j])*0.25
                error[i,j] = np.abs(A_new[i,j] - A[i,j])
        
        A = A_new.copy()
        k += 1
        errmax = np.max(np.max(error))

        if errmax < epsilon:
            break
    return A

if __name__ == '__main__':

    args = get_args(sys.argv[1:])
    A = local_jacobi(args)


    cmd = 'make run'
    cmd += ' T=' + str(args.t)
    cmd += ' B=' + str(args.b)
    cmd += ' R=' + str(args.r)
    cmd += ' L=' + str(args.l)
    cmd += ' I=' + str(args.i)
    cmd += ' N=' + str(args.n)
    print cmd
    tic = time.time()
    pid = subprocess.Popen(cmd, shell=True, 
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE)
    pid.wait()
    processtime = time.time() - tic
    out, err = pid.communicate()
    
    B = np.zeros(A.shape)
    i = 0
    for line in out.split('\n'):
        #print line
        if line.startswith(args.l):
            for j, tmp in enumerate(line.split(',')):
                #print i,j, tmp
                B[i, j] = float(tmp)
            i += 1

    N = A.shape[0]     
    print 'difference: ', np.abs(np.sum(A[1:N-1,1:N-1]-B[1:N-1,1:N-1]))
    print 'time:       ', processtime


    