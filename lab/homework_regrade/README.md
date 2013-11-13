# Homework Regrade

We have opened up the following directories:

	/curc/tools/grading/assignments/assignment-01/$USER
	/curc/tools/grading/assignments/assignment-02/$USER
	/curc/tools/grading/assignments/assignment-04/$USER
	/curc/tools/grading/assignments/assignment-05/$USER

If please go into your directory and run the grading scripts.  You may make changes to *any part of your code* you'd like.  On **November 15th, 2013** we will close the directories and regrade the assignments.  No changes will be accepted after that point.  Sound fair?

## Grading scripts

Navigate to your directory.

	/curc/tools/grading/assignments/assignment-01/$USER

Type

	/curc/tools/grading/scripts/check 1
	/curc/tools/grading/scripts/check 2
	...
	/curc/tools/grading/scripts/check N

where *N* is the number of trials for the homework.  Each assignment has a different number of trials.  The first trial, N=1, always runs **make clean** and **make**, and there are not points given.  Make sure you run this trial in your directory.

## Example: Assignment 5

Go to assignment 5 directory.

	cd /curc/tools/grading/assignments/assignment-05/$USER

Run the first trial.

	/curc/tools/grading/scripts/check 1

As long as your executables are built, it doesn't matter what the output says.  No points either way.


Run trial number 2.

	/curc/tools/grading/scripts/check 2

output

	checking assignment 5: trial  2 of 11

	75 points total

	PASSED +1 point: file exists
	PASSED +1 point: compiles without openmp
	PASSED +1 point: compiles with openmp
	PASSED +1 point: Correct number of threads
	PASSED +1 point: Output looks correct

Run trial number 8.

	/curc/tools/grading/scripts/check 8

output

	checking assignment 5: trial  8 of 11

	75 points total

	accuracy:    1.56723217515e-07
	1 thread:    3.52934908867
	12 threads:  0.433204889297
	speedup:     8.1470666095
	threshold:     7

	PASSED +5 points:  accuracy
	PASSED +5 points:  for speedup 7

## How we grade?

You can simulate the grading by combining all the checks in a single script.

	/curc/tools/grading/scripts/check 1
	/curc/tools/grading/scripts/check 2
	/curc/tools/grading/scripts/check 3
	/curc/tools/grading/scripts/check 4 
	/curc/tools/grading/scripts/check 5
	/curc/tools/grading/scripts/check 6
	/curc/tools/grading/scripts/check 7
	/curc/tools/grading/scripts/check 8
	/curc/tools/grading/scripts/check 9
	/curc/tools/grading/scripts/check 10
	/curc/tools/grading/scripts/check 11
	/curc/tools/grading/scripts/check 12
	/curc/tools/grading/scripts/check 14
	/curc/tools/grading/scripts/check 15

Run this in each directory, send the output to a file, and count up the number of points you have.

- Assignment 01
	* 2: 10 points total 

- Assignment 02
	* **46 points**
	* 2 - 4: 3 points each 
	* 5 - 9: 3 points each
	* 10: 13 points
	* 9 points free 

- Assignment 04
	* **100 points**
	* 2-6: 5 points each (25)
	* 7-11: 15 points each (75)
		- 5 accuracy
		- 5 efficiency
		- 5 better efficiency

- Assignment 05
	* **75 points**
	* 2: hello world 5 points
	* 3-6: 5 points each (20)
	* 7-11: 10 points each (50)
		- 5 accuracy
		- 5 efficiency

- Assignment 06
	* 175 points **
	* 2-4: simpsons 10 points each (40)
		- 5 accuracy
		- 5 efficiency
	* 6-8: openmp 10 each (30)
		- 5 accuracy
		- 5 efficiency
	* 9-11: mpi 10 each (30)
		- 5 accuracy
		- 5 efficiency
	* 12-14: hybrid 10 each (30)
		- 5 accuracy
		- 5 efficiency
	* pdf report: 45 points


