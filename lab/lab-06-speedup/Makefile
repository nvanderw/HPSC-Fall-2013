CC=gcc
CFLAGS=-std=gnu99 -Wall -Ofast -fopenmp

all: jacobi matrix_multiply

jacobi: jacobi.c
	$(CC) $(CFLAGS) -o $@ $<

run: jacobi
	./jacobi $(T) $(B) $(L) $(R) $(I) $(N)

matrix_multiply: matrix_multiply.c
	$(CC) $(CFLAGS) -o $@ $< -lhdf5

clean:
	rm -f *.o jacobi matrix_multiply
