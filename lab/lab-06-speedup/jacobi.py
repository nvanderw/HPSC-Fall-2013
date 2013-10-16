#!/usr/bin/env python
import numpy as np
import sys

T,B,R,L,I,N = sys.argv[1:]
print sys.argv
nx = int(N)
ny = int(N)
A = np.zeros((nx, ny))

A[0,:] = int(T)   # Top
A[nx-1,:] = int(B)  # Bottom
A[:,0] = int(L)     # Left
A[:,ny-1] = int(R)    # Right

A_new = A.copy()
error = A_new - A

k = 1 
epsilon = 1e-5
max_iterations = int(I)
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

for i in np.arange(0,nx):
	for j in np.arange(0,ny-1):
		print str(A[i,j]) + ',',
	print str(A[i,ny-1])


