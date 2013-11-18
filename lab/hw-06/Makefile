CFLAGS=-std=gnu99 -O3
CC=gcc
EXECUTABLES=mpi_simpsons jacobi_mpi jacobi_omp jacobi_hybrid

all: $(EXECUTABLES)

mpi_simpsons: mpi_simpsons.c
	mpicc -o mpi_simpsons -std=gnu99 -O3 mpi_simpsons.c

jacobi_omp: jacobi_omp.c
	$(CC) -o jacobi_omp $(CFLAGS) -fopenmp jacobi_omp.c

jacobi_mpi: jacobi_mpi.c
	mpicc -o jacobi_mpi $(CFLAGS) jacobi_mpi.c

jacobi_hybrid: jacobi_mpi.c
	mpicc -o jacobi_hybrid $(CFLAGS) -fopenmp jacobi_mpi.c

simpsons: mpi_simpsons
	mpirun -np $(P) mpi_simpsons $(L) $(U) $(F) $(N)

openmp: jacobi_omp
	./jacobi_omp $(T) $(B) $(L) $(R) $(I) $(N) $(O)

mpi: jacobi_mpi
	mpirun -np $(P) jacobi_mpi $(T) $(B) $(L) $(R) $(I) $(N) $(O)

hybrid: jacobi_hybrid
	mpirun -np $(P) -bynode -bind-to-core jacobi_hybrid $(T) $(B) $(L) $(R) $(I) $(N) $(O)

clean:
	rm -f $(EXECUTABLES)
