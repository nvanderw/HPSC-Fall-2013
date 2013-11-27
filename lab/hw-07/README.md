#Homework 7: The Last One

This assignment is due **Friday, December 6th**.

There will **NOT** be a late period for this assignment.

In this assignment you will implement [Cannon's algorithm](http://en.wikipedia.org/wiki/Cannon's_algorithm) for parallel matrix multiply using [parallel IO with HDF5](http://www.nersc.gov/users/training/online-tutorials/introduction-to-scientific-i-o/?show_all=1#toc-anchor-4). 

##Interface

Please have your *Makefile* accept 4 values:

- **A**: name of matrix A in HDF5 format.
- **B**: name of matrix B in HDF5 format.
- **C**: name of the output file.  This should be in the same format as the previous assignments.
- **P**: number of processors to use.

If we type

	make run A=matrix_a.hdf5 B=matrix_b.hdf5 C=output_test1.hdf5 P=10

your *Makefile* should call `mpirun -np $(P)` and pass the executable with the values of A, B, and C. 

##File IO

We will be running code in the lustre scratch directory for this assignment.  You will need to submit and check your code from here:

	/lustre/janus_scratch/hpcs/assignments/assignment-07/$USER

This directory is available.

##Grading

You will be graded on three characteristics:

- **threshold**: there will be a minimum accetable runtime for this assignemnt (very easy to make if you have a working algorithm) and you must run under this threshold to get any points.
- **accuracy**: you must get the correct answer to get any points.
- **runtime**: you will recieve additional points if you are below a performance runtime.
- **performance**: **10%** of your grade will be based on your rank in the class.

There has been at least one serious complaint aobut the check scripts.  As a result, I am not planning on providing a check script for this assignment.  We will discuss this after the break.  You should test your implementation carefully to make sure that:

- You can run problems using `P={4,9,16,25,36,64}` processors.  Note I will only test your code using these P values.
- The HDF5 files I will be using are similar to the previous assignemnts, specifically using `double` values and a dataset called `DATASET`. Please see [this example](https://github.com/ResearchComputing/HPSC-Fall-2013/blob/master/lab/lab-12-file_io/create_file.py) for details on how to create this file.
- I will be testing your code with any size, yet square, matrix.
- I will be using your `Makefile` to test your code.  Please make sure it works. 
	- 'make' must build your code.
	- The 'make run' command must execute the trail.
	- We will not be fixing Makefiles or code to make this process work.
	- There will be no partial credit for failed builds or trials.
- This assignment will have 10 tests, each worth 15 points:
	- 5 points for minimum threshold
	- 5 points for accuracy
	- 5 points for performance threshold
- The final 15 points will be based on your rank.

If you have any questions, please ask!

Also, please start early as this  is a difficult assignment too.

##Advice

- If you don't have a working program, you will not get any points.
- With that in mind, start simple and make sure you have the parallel algorithm down first.
	- e.g. You can put off the parallel IO until the end if necessary.
- The gather method that you wrote for the previous assignment should drop right in.  This is a great debugging tool.



