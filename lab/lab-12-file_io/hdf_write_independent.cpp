/*  
 *   *  This example writes data to the HDF5 file by rows.
 *    *  Number of processes is assumed to be 1 or multiples of 2 (up to 8)
 *     */
#include <iostream>
#include <stdlib.h>   
#include "mpi.h"
#include "hdf5.h"


int main (int argc, char **argv)
{

    // Get the file
    if(argc < 2){
        std::cout << "please specify a file" << std::endl;
        return 0;
    }
    std::string filename = argv[1];

    // Get MPI into this
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  

    int m_rows = 12;
    int m_cols = 8;
    if( argc == 4 )
    {
        m_rows = atoi(argv[2]);
        m_cols = atoi(argv[3]);
    }


    herr_t	status;

    //Set up file access property list with parallel I/O access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    //int use_gpfs_hints=1;
    //H5Pset_fapl_mpiposix(plist_id, MPI_COMM_WORLD, use_gpfs_hints);

    //Create a new file collectively and release property list identifier.
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
   
    //Create the dataspace for the dataset.
    hsize_t  dimsf[2] = {m_rows, m_cols};
    hid_t filespace = H5Screate_simple(2, dimsf, NULL); 

    //Create the dataset with default properties and close filespace.
    hid_t property_id = H5Pcreate (H5P_DATASET_CREATE);
    status = H5Pset_layout (property_id, H5D_CONTIGUOUS);

    hid_t dset_id = H5Dcreate(file_id, "DATASET", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, property_id, H5P_DEFAULT);
     
    H5Sclose(filespace);

    //hyperslab selection parameters
    hsize_t count[2] = {};  
    hsize_t offset[2] = {};

    hid_t memspace = H5Screate_simple(2, count, NULL);

    //Select hyperslab in the file.
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    //Initialize data buffer 
    const int length = count[0]*count[1];
    double * data = new double[length];
    for (int i=0; i < length; i++) {
        data[i] = rank;
    }

    //Create property list for collective dataset write.
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
    
    delete [] data;


    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Fclose(file_id);
 
    MPI_Finalize();

    return 0;
}     
