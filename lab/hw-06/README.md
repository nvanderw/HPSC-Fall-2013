#Homework 6: Jacobi OpenMP, MPI, and Hybrid

This assignment is due on **Friday, November 15th, 2013**.

There are four parts to assignment 6.

- [Simpson's rule](https://github.com/ResearchComputing/HPSC-Fall-2013/tree/master/lab/lab-08-integrate): Please submit your code with this assignment instead of assignment-08.

- [Jacobi with OpenMP](https://github.com/ResearchComputing/HPSC-Fall-2013/tree/master/lab/lab-06-speedup): We are going to grade your code here instead of the previous submission in assignment-06.

- **MPI**: Write an MPI Jacobi implementation.

- **Hybrid OpenMP-MPI**: Write an implementation that uses both OpenMP and MPI.

**NOTE**: The `N` that we pass includes the boundary.

## Interface

We are changing the interface slightly to accommodate all the different programs.

### Simpson's

We will test your Simpson's code by typing:

	make simpsons L=5 U=10 F=a N=1000 P=24

**NOTE** this has changed from **make run** to **make simpsons**.

Also, if you would like to resubmit your report you may do so.

### OpenMP

We will test your OpenMP Jacobi implementation by typing:

	export OMP_NUM_THREADS=12; make openmp N=10 L=2 T=2 R=1 B=1  I=10000 O=output_2.csv

**NOTE** the following changes:

- Command changed from **make run** to **make openmp**.
- There is an additional argument O for an output file that you should write your results to.  This should avoid some of the issues with the screen dump.  

### MPI

We will test your MPI implementation using the following command:

	export OMP_NUM_THREADS=1; make mpi N=11 L=2 T=2 R=1 B=1  I=20000 O=output_2.csv P=9

This is the same arguments used for the OpenMP version with one exception: we are now passing the variable `P`, which is the number of processors.  As with the Simpson's rule, you need to call **mpirun** in your make file.

	mpirun -np $(P)$ your_mpi_implementation
	
You are required to use the **grid** implementation described in class.  We will look at code to verify.  No credit given for other implementations.


### Hybrid

We will test your MPI implementation using the following command (example):

	export OMP_NUM_THREADS=3; make hybrid N=1002 L=2 T=2 R=1 B=1  I=20000 O=output_2.csv P=4

This is the same arguments used for the MPI versions with one important Makefile implementation detail: your hybrid command must call mpirun in the following way.

	mpirun -np $(P)$ -bynode -bind-to-core your_mpi_implementation

Again, use the **grid** implementation from above.

## Grading

Grading instructions will be posted on *Wednesday, November 30th*.

## Submission

Please use the new submission instructions. Your directory is located here:

	/curc/tools/grading/assignments/assignment-06/$USER

At the minimum, you should have:

- Makefile: make, make clean, make run
- ALL your source code
- **A report in PDF format for the OpenMP version**.
- Please no tar, zip, or bzip files. From your test directory, just type:

		cp * /curc/tools/grading/assignments/assignment-06/$USER

NOTES:

make clean should remove all executables.
