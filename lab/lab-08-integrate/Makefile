
P=2

simpsons: simpsons.c
	mpicc -o simpsons simpsons.c -std=gnu99 -lm
run: simpsons
	mpirun -np $(P) ./simpsons $(L) $(U) $(F) $(N)
clean:
	rm simpsons
