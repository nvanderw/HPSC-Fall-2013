#Lab 05: OpenMP

## Homework to turn in

This assignment *should* be quick.  We will do more with OpenMP in the next few weeks as we talk more about parallelism.

__DUE: October 11, 2013__

## In lab discussion

Please submit the following programs (hello_worlrd.c and the corrected openmp programs) to the dropbox on D2L

### hello_world.c (5 points)

Write a program that prints the following output

Total number of tasks: #   
Hello world from tasks: # 

Can you also make it so that it also works with a non-openmp compiler

### Bug 1-4 (5 points for each bug fix)

- The first program doesn't compile.  Getting it to compile is sufficient.
- The next three either seg-fault or give the wrong answer.  You can figure out the correct results by runnign the code on a single thread.

        export OMP_NUMP_THREADS=1

        
- The second program computes the wrong results
- The third program will segmentation fault
- The fourth program will give you a wrong result (also, you will need to pass `sum` to the function to get the correct result).

Download the code and fix the bugs.

##Matrix Multiply (10 points)

Submitt per instructions below

Take you matrix multiply from the assignment four and make it even faster by adding shared-memory parallelism with OpenMP.  As in [lab 04](https://github.com/ResearchComputing/HPSC-Fall-2013/tree/master/lab/lab-04), please test your code with the `compare.py` file.  We will be using this to test your code also.
**Make it easier on yourself**: don't worry about blocking and simply parallelize the naive with the correct loop order.

Use the following module for building your code:

        module load HPSC_CLASS/lab05

## Grading

We will run your code on 5 different problem sizes using 1,6, and 12 threads.  Any test that gives the incorrect result will get zero points.  You will be graded on your accuracy and relative speed of the multi-thread examples.  Each test is worth 10 points for a total of 50 points.  

## What to turn in?

Please use the new [submission instructions](https://github.com/ResearchComputing/HPSC-Fall-2013/tree/master/lab/submission_instructions).  Your 
directory is located here:
        
        /curc/tools/grading/assignments/assignment-05/$USER

At the minimum, you should have a:

- Makefile
- At least one source file




