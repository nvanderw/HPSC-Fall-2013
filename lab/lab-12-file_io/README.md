# Scientific IO

This is based on the excellent tutorial from [NERSC](http://www.nersc.gov/users/training/online-tutorials/introduction-to-scientific-i-o/?show_all=1).

We use file IO for:

- Writing output to a file
- Reading input for analysis
- Implementing **out-of-core** algorithms

## Methods for File IO

- File-pre-process
- Shared-file
- Collective buffering

<img src="http://www.nersc.gov/assets/Training-Materials/HowisonParallelIOFolder/_resampled/resizedimage600564-io-types.png">

###File-per-process

- Very efficient.
- Metadata for a large collection of files introduces a bottleneck.
- Hard to manage.
- Requires same number of processors each time.

###Shared-file

- Many processors share the same file.
- Write to exclusive regions.
- HDF5: portable and independent of the number of processors used.

###Collective buffering

- Only a subset of processors interact with the file.
- Reduces contention.

##Lustre

- Parallel file system
- Store a single file across several disks and servers to increase parallel performance.

### Basics

- **stripe** is a unit of storage in Lustre.
- Each stripe has a default **size**.
- **striping** divides up a shared file across many servers and disks.

<img src="http://www.nersc.gov/assets/Training-Materials/HowisonParallelIOFolder/_resampled/resizedimage600152-lustre-arch.png">

When you stripe a file you increase bandwidth to that file. The total bandwidth of the shared, striped file is the aggregate of the individual banwidths.

Janus

- 1 **Metadata Server (MDS)** 
- 60 **Object Target Servers (OTS)** (I/O server)

### The `lfs` Utility

List of the OTSs:

	lsf osts

What is the current strip setting?

	lsf getstripe

Set stripe settings:

	lfs setstripe [--size|-s stripe_size] [--count|-c stripe_cnt]  <dirname|filename>

Example:

	bash-janus> mkdir temp_dir

	bash-janus> lfs setstripe -c 3 -s 32m temp_dir

	bash-janus> touch temp_dir/temp_file

	bash-janus> lfs getstripe temp_dir

	temp_dir
	stripe_count:   3 stripe_size:    33554432 stripe_offset:  -1 

	temp_dir/temp_file
	lmm_stripe_count:   3
	lmm_stripe_size:    33554432
	lmm_stripe_offset:  18
		obdidx		 objid		objid		 group
		    18	      12787913	     0xc320c9	             0
		     7	      12863377	     0xc44791	             0
		    23	      12496893	     0xbeaffd	             0

The size is *33554432* **bytes** (32mb) and the count is *3*.  So this file is striped across 3 Object Target Servers with a size of 32mb.

Issues:

- Striping parameters must be set before a file is written to.
- One solution is to touch the file, set the stripe, then write.

##Parallel HDF5

**Hierarchical Data Format v5**

HDF5 provides an abstracted filesystem-within-a-file.
	
- Portable to any system with the HDF5 library.
- Multiple groups and datasets.

<img src="http://www.nersc.gov/assets/Training-Materials/HowisonParallelIOFolder/_resampled/resizedimage600330-hdf5-file.png">

### Simple Read and Write

#### Read

	std::string filename = argv[1];
    herr_t status;

Create a new file using the default properties.

    hid_t file_id = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t dataset_id = H5Dopen(file_id, "DATASET", H5P_DEFAULT);
    hid_t space_id = H5Dget_space(dataset_id);
    
Get the length of the dataset.

    int length = H5Sget_simple_extent_npoints(space_id);

Get the dimensions of the dataset.  Here we know it's 2D.

    hsize_t  dims[2];
    hsize_t  mdims[2];
    status = H5Sget_simple_extent_dims(space_id,dims,mdims);
    
    int m_rows = dims[0];
    int m_cols = dims[1];

Allocated space in `data` and read the dataset.

    double * data = new double[m_rows*m_cols];
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

Close the resources.

    status = H5Sclose(space_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);

#### Write

Create some data.

	int m_rows = 12;
    int m_cols = 8;
    double * data = new double[m_rows*m_cols];

    // Fill with some data
    for(int r=0; r<m_rows; ++r)
        for(int c=0; c<m_cols; ++c)
            data[r + c*m_rows] = r + c*m_rows+1.0;

Create a new file using the default properties.

    hid_t file_id = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

Create dataspace.  Setting maximum size to NULL sets the maximum size to be the current size.
  
    herr_t status;
    hsize_t  dims[2] = {m_rows,m_cols};
    hid_t space_id = H5Screate_simple (2, dims, NULL);

Create the dataset creation property list, set the layout to contiguous.

    hid_t property_id = H5Pcreate (H5P_DATASET_CREATE);
    status = H5Pset_layout (property_id, H5D_CONTIGUOUS);

Create the dataset. 

    hid_t dataset_id = H5Dcreate (file_id, "DATASET", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, property_id, H5P_DEFAULT);
   
Write the data to the dataset.

    status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

Clean up.

    status = H5Sclose(space_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);
    status = H5Pclose(property_id);

### Offset, Count, Stride

- Offset: specifies the starting coordinates of the hyperslab to select.
- Count: specifies how many positions to select from the dataspace in each dimension.
- Stride: specifies which elements along a dimension are selected.

Each process defines:

- the memory hyperslab
- the file hyperslab

Each process executes a partial write or read call.

- collective
- independent

<img src="http://davis.lbl.gov/Manuals/HDF5-1.4.3/Tutor/image/pimg024.gif">

### Examples

#### Partial Read

Create file space. 
    int m_rows = 4;
    int m_cols = 4; 
	  hsize_t offset[2] = {  What goes here? };
   	hsize_t count[2] = { what goes here? };
    status = H5Sselect_hyperslab (space, H5S_SELECT_SET, offset, NULL, count, NULL);

Create memory space.

    hsize_t dim_out[2] = {m_rows,m_cols};
    hsize_t offset_out[2] = {what goes here?;
    hsize_t count_out[2] = {what goes here?};
    hid_t memspace = H5Screate_simple(2,dim_out,NULL); 

    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);

Read the data from the file into the memory.

    status = H5Dread (dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, data);
    
To run this example on JANUS, you can execute:

	./hdf_write test

	./hdf_read_partial test
      20      21      22      23
      28      29      30      31
      36      37      38      39
      44      45      46      47

    ./hdf_read_partial test 0 0
       1       2       3       4
       9      10      11      12
      17      18      19      20
      25      26      27      28

#### Parallel Independent Write

Set up file access property list with parallel I/O access.

    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

Create a new file collectively and release property list identifier.
    
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

Assume we want our file to look like the following when I run with 10 processors.

    python print_file.py test.h5
    (12, 8)
    [[ 0.  0.  0.  0.  0.  0.  0.  0.]
     [ 1.  1.  1.  1.  1.  1.  1.  1.]
     [ 2.  2.  2.  2.  2.  2.  2.  2.]
     [ 3.  3.  3.  3.  3.  3.  3.  3.]
     [ 4.  4.  4.  4.  4.  4.  4.  4.]
     [ 5.  5.  5.  5.  5.  5.  5.  5.]
     [ 6.  6.  6.  6.  6.  6.  6.  6.]
     [ 7.  7.  7.  7.  7.  7.  7.  7.]
     [ 8.  8.  8.  8.  8.  8.  8.  8.]
     [ 9.  9.  9.  9.  9.  9.  9.  9.]
     [ 0.  0.  0.  0.  0.  0.  0.  0.]
     [ 0.  0.  0.  0.  0.  0.  0.  0.]]


How should I setup the count and offset?

	hsize_t count[2] = { What goes here? };   
    hsize_t offset[2] = {What goes here? };

Simple 2D memory space.

    hid_t memspace = H5Screate_simple(2, count, NULL);

Select hyperslab in the file.

    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

Finally, write.  

	plist_id = H5Pcreate(H5P_DATASET_XFER);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);

### Collective

	filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

Exclude certain ranks.

    if (mpi_rank == 3) H5Sselect_none(filespace);

Speficy collective and write.

 	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);














