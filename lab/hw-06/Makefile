mpi_simpsons: mpi_simpsons.c
	mpicc -o mpi_simpsons -std=gnu99 -O3 mpi_simpsons.c

simpsons: mpi_simpsons
	mpirun -np $(P) mpi_simpsons $(L) $(U) $(F) $(N)
