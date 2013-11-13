#Lab 08: Simpson`s Rule

__DUE: October 25th, 2013__

Simpson`s rule is a widely-used technique for integrating an arbitrary function over specific bounds.  For more details, please see the [wikipedia article](http://en.wikipedia.org/wiki/Simpson's_rule).

Here's the [notebook](http://nbviewer.ipython.org/urls/raw.github.com/ResearchComputing/HPSC-Fall-2013/master/lab/lab-08-integrate/simpsons.ipynb).


You will implement an `MPI` version of Simpson's rule capable of running on any number of **even** processors.  Your program will use the functions defined in the `functions.h` file.  We will pass you the following parameters via your `make run` command.

- `L`: lower bound
- `U`: upper bound
- `F`: [a,b,c,d]
	- a: `func_a`
	- b: `func_b`
	- c: `func_c`
	- d: `func_d`
- `N`: number of iterations
- `P`: number of processors, for example see `Makefile`

		mpirun -np $(P$) your_exe

## Grading

Please use this module:

	module load HPSC_CLASS/lab08 

	login04 >which mpicxx
	/curc/tools/x_86_64/rh6/openmpi/1.6.4/intel/13.0.0/torque/4.1.4/ib/bin/mpicxx


We will call your `make run` command from `python` (as usual):

	make run L=5 U=10 F=a N=1000 P=24

Your ouput shoule be parsable and look like the following:

	area = 44.222211332

## What to turn in?

Please see [lab 10]() for new submission instructions.

