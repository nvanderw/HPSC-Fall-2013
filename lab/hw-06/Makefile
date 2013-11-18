CFLAGS=-std=gnu99 -O3
CC=gcc
EXECUTABLES=mpi_simpsons jacobi_mpi jacobi_omp jacobi_hybrid

all: $(EXECUTABLES)

mpi_simpsons: mpi_simpsons.c
	mpicc -o mpi_simpsons -std=gnu99 -O3 mpi_simpsons.c

simpsons: mpi_simpsons
	mpirun -np $(P) mpi_simpsons $(L) $(U) $(F) $(N)

jacobi_omp: jacobi_omp.c
	$(CC) -o jacobi_omp $(CFLAGS) -fopenmp jacobi_omp.c

jacobi_mpi: jacobi_mpi.c
	mpicc -o jacobi_mpi $(CFLAGS) jacobi_mpi.c

jacobi_hybrid: jacobi_mpi.c
	mpicc -o jacobi_hybrid $(CFLAGS) -fopenmp jacobi_mpi.c

clean:
	rm -f $(EXECUTABLES)
