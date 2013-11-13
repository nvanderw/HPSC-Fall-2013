#Homework 7: The Last One

This assignment is due **Friday, December 6th**.

In this assignment you will implement [Cannon's algorithm](http://en.wikipedia.org/wiki/Cannon's_algorithm) for parallel matrix multiply using [parallel IO with HDF5](http://www.nersc.gov/users/training/online-tutorials/introduction-to-scientific-i-o/?show_all=1#toc-anchor-4). 

##Interface

Please have your *Makefile* accept 4 values:

- **A**: name of matrix A in HDF5 format.
- **B**: name of matrix B in HDF5 format.
- **C**: name of the output file.  This should be in the same format as the previous assignments.
- **P**: number of processors to use.

If we type

	make A=matrix_a.hdf5 B=matrix_b.hdf5 C=output_test1.hdf5 P=10

your *Makefile* should call `mpirun -np $(P)` and pass the executable with the values of A, B, and C. 

##File IO

We will be running code in the lustre scratch directory for this assignment.  You will need to submit and check your code from here:

	/lustre/janus_scratch/hpcs/assignments/assignment-07/$USER

This directory will not be setup until Friday, November 22nd.

##Grading

You will be graded on three characteristics:

- **threshold**: there will be a minimum accetable runtime for this assignemnt (very easy to make if you have a working algorithm) and you must run under this threshold to get any points.
- **accuracy**: you must get the correct answer to get any points.
- **runtime**: you will recieve additional points if you are below a performance runtime.
- **performance**: **10%** of your grade will be based on your rank in the class.  I'll tell you your rank when you run the check script.

I will have the check scripts ready by Friday, November 22nd.





