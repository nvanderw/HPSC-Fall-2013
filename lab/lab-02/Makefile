GCC = gcc -std=gnu99
ICC = icc
PI = pi.c
PIEXE = pi

VA = matrixA
VB = matrixB
VC = matrixC

intel-slow: $(PI)
	$(ICC)  -O0 -no-vec $(PI) -o $(PIEXE)
	./$(PIEXE)

intel-novec: $(PI)
	$(ICC) -O2 -no-vec $(PI) -o $(PIEXE)
	./$(PIEXE)

intel-fast: $(PI)
	$(ICC) -Ofast $(PI) -o $(PIEXE)
	./$(PIEXE)

gcc-slow: $(PI)
	$(GCC) -O0 $(PI) -o $(PIEXE)
	./$(PIEXE)

gcc-novec: $(PI)
	$(GCC) -O2 $(PI) -o $(PIEXE)
	./$(PIEXE)	

gcc-fast: $(PI)
	$(GCC) -O3 $(PI) -o $(PIEXE)
	./$(PIEXE)

python-matrix-multiply: matrix_multiply.py
	python matrix_multiply.py -a $(VA) -b $(VB) -c $(VC)

matrix_multiply: matrix_multiply.c
	$(GCC) -O3 -o matrix_multiply -std=gnu99 matrix_multiply.c

matrix-multiply: matrix_multiply
	./matrix_multiply $(VA) $(VB) > $(VC)

cpuinfo:
	cat /proc/cpuinfo | grep 'cpu MHz'

check:
	/curc/tools/grading/scripts/check

clean:
	rm -f $(PIEXE) matrix_multiply
