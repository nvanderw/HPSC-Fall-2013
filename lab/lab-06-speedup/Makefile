T=2
B=1
R=1
L=2
I=1000
N=8

CC=gcc
CFLAGS=-std=gnu99 -Wall -Ofast -fopenmp

all: jacobi matrix_multiply

jacobi: jacobi.c
	$(CC) $(CFLAGS) -o $@ $<

run: jacobi
	./jacobi $(T) $(B) $(L) $(R) $(I) $(N)

matrix_multiply: matrix_multiply.c
	$(CC) $(CFLAGS) -o $@ $< -lhdf5

check: compare_jacobi.py 
	python compare_jacobi.py -t=$(T) -b=$(B) -r=$(R) -l=$(L) -i=$(I) -n=$(N)

clean:
	rm -f *.o jacobi matrix_multiply

report.pdf:
	pdflatex report.tex
