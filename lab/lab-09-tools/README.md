# Lab 09: The Tools

## Steps to Optimal Performance

- **Verification**
	- Are you getting the correct answer?
- **Serial Optimization**
	- Have you optimized your serial code?
	- Flags: -O3, -Ofast ...
	- Cache: Are you using the best serial algorithm?
	- Memory: Do you have any memory leaks?
	- Profile: Will parallelization help? Where?
- **Parallel Optimization**
	- Are you still getting the correct answer?
	- Where are the bottlenecks?

## Tool set

- **gprof**: profiling
- **psrun**: cache, pofiling
- **valgrind**: memory leaks
- **HPCToolkit** parallel profiling
- **TAU**: parallel profiling
- **TotalView**: parallel debugging


### [gprof](https://sourceware.org/binutils/docs/gprof/)

- Compile your code with the `-pg` option.  
- Run your code, it will produce a file called **gmon.out**.
	
		./<your_application>

- Interpret your results:

		gprof <your_application> gmon.out

### [psrun](http://perfsuite.ncsa.illinois.edu/psrun/)
	
	module load perfsuite/perfsuite-1.1.2_openmpi-1.6.4_intel-13.0.0_ib
	psrun <your_application>
	psprocess <your_application>.*.xml

### [valgrind](http://valgrind.org/)

	module load valgrind/valgrind-3.8.1
	valgrind --tool=memcheck <your_serial_application>

### [HPCToolkit](http://hpctoolkit.org/)

https://gist.github.com/mlunacek/7118061

### [TAU](http://www.cs.uoregon.edu/research/tau/home.php)

https://gist.github.com/mlunacek/7118492

### [TotalView](http://www.roguewave.com/products/totalview.aspx)

(check back)