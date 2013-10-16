#Lab 06: Speedup, Efficiency, and Karp-Flatt

__DUE: October 15th, 2013__

You will compare the speedup, efficiency, and the Karp-flatt metric for two OpenMP-enabled codes:

* Matrix multiply
* Jacobi iterative solution to the heat equation.

## Jacobi Iteration

Please create a `Makefile` that accepts a `run` target and 6 parameters: 

* T: top 
* B: bottom
* L: left
* R: right
* I: max iterations
* N: matrix size

If we wanted to test your code based on the example from lecture, we will typing the following:

	make run T=2 B=1 R=1 L=2 I=10 N=4

Your output should be a matrix that represent the solution.  In the previous example, this would look like the following:

	2.0,  2.0,  2.0,  1.0 
	2.0,  1.75,  1.5,  1.0
	2.0,  1.5,  1.25,  1.0
	2.0,  1.0,  1.0,  1.0 
	
Note that the `N=4` includes the entire matrix and that the boundary conditions are included in the output.

You should use a value of `epsilon=1e-5` as a stopping criteria.

Use the following module for building your code:

        module load HPSC_CLASS/lab06

## Report

Perform a scaling study for each code.  Specifically run and time each code based on the following parameters:

Matrix Multiply

* 1000 x 1000 problem size
* Threads: 1,2,4,8,12

Jacobi

* make run T=2 B=1 R=2 L=1 I=20000 N=1000
* Threads: 1,2,4,8,12

Create three plots:

* Speedup: Plot the speedup of both algorithms on a single plot. 
* Efficiency: Plot the efficiency of both algorithms on a single plot.
* Karp-Flatt: Both algorithms on one plot.

Please clearly label your graphs and provide a brief explanation of what you observe.


## Grading

Your report will be worth 10 points per plot plus 5 for the description.  This is a total of 45 points. 
We run 3 tests on your jacobi code to verify it works.  The is worth 30 points.  This will be based on accuracy
and speedup.

Total: 75 point

I have added a `Makefile` and a `compare_jacobi.py`.  We will test your code by running the following command several times:

	python compare_jacobi.py -t=2 -b=1 -r=1 -l=2 -i=1000 -n=10
	
We will change the input parameters.  It is necessary that your `Makefile` is working as the `compare_jacobi.py` calls

	make run 

Here are a few examples from my Makefile.

With `N=10`

	login01 >make check N=10 B=1 L=1 I=1000
	python compare_jacobi.py -t=2 -b=1 -r=1 -l=1 -i=1000 -n=10
	make run T=2 B=1 R=1 L=1 I=1000 N=10
	difference:  4.35576019697e-11
	time:        0.498589992523

With `N=20`

	login01 >make check N=20 B=1 L=1 I=1000
	python compare_jacobi.py -t=2 -b=1 -r=1 -l=1 -i=1000 -n=20
	make run T=2 B=1 R=1 L=1 I=1000 N=20
	difference:  4.12878620182e-11
	time:        4.25668883324


## What to turn in?

Please use the new submission instructions. Your directory is located here:

    /curc/tools/grading/assignments/assignment-06/$USER

At the minimum, you should have:

- **Makefile**.  Include `make` to build your code, `make clean` to remove object files and executables, and `make run` as described above.
- **All the source code** (at least one source file) for the heat equation.
- **A report in PDF format**.

Please do not check in a tar file.  Place all file at the root of your directory `assignment-06/$USER`.





