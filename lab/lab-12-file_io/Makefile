#module load hdf5/hdf5-1.8.11_szip-2.1_zlib-1.2.78_jpeglib-8d_openmpi-1.6.4_gcc-4.7.2_ib

CC=mpicxx
OPT=-O3
CFLAGS=-c

HDF5_INCLUDE=/curc/tools/x_86_64/rh6/hdf5/1.8.11/szip/2.1/zlib/1.2.78/jpeglib/8d/openmpi/1.6.4/gcc/4.7.2/torque/4.1.4/ib/include
HDF5_LIB_DIR=/curc/tools/x_86_64/rh6/hdf5/1.8.11/szip/2.1/zlib/1.2.78/jpeglib/8d/openmpi/1.6.4/gcc/4.7.2/torque/4.1.4/ib/lib
HDF5_LIBS=hdf5

EXES=hdf_read hdf_write hdf_read_partial hdf_write_independent

all: $(EXES) 

%: %.cpp
	$(CC) $(OPT) $< -I $(HDF5_INCLUDE) -L $(HDF5_LIB_DIR) -l $(HDF5_LIBS) -o $@

parallel: all
	mpirun -np 10 ./hdf_write_independent test.h5
	python print_file.py test.h5

serial: all
	./hdf_write test.h5
	python print_file.py test.h5
	./hdf_read test.h5
	./hdf_read_partial test.h5 0 0
	./hdf_read_partial test.h5 4 4	

clean:
	rm -rf $(EXES) *.h5
