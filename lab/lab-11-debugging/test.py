#!/usr/bin/env python
from mpi4py import MPI
import parallel_jacobi_working as pj

if __name__ == '__main__':

    jac = pj.ParallelJacobi()

    # jac.iterate()
    # jac.communicate()
    # jac.printstate()

    jac.test_gather()
    jac.test_shift()