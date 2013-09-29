#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "hdf5.h"
#include "timing.h"

namespace util
{

void read_hdf(const std::string &filename, double* &data, int &m_rows, int &m_cols){

    hid_t file_id, dataset_id, space_id, property_id; 
    herr_t status;

    //Create a new file using the default properties.
    file_id = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset_id = H5Dopen(file_id, "x", H5P_DEFAULT);
    space_id = H5Dget_space(dataset_id);
    int length = H5Sget_simple_extent_npoints(space_id);
    hsize_t  dims[2];
    hsize_t  mdims[2];
    status = H5Sget_simple_extent_dims(space_id,dims,mdims);
    m_rows = dims[0];
    m_cols = dims[1];
    
    data = new double[length];
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, data);

    status = H5Sclose(space_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);

}

void write_hdf(const std::string &filename, double* &data, int &m_rows, int &m_cols){

    hid_t file_id, dataset_id, space_id, property_id; 
    herr_t status;

    hsize_t  dims[2] = {m_rows,m_cols};
   
    
    //Create a new file using the default properties.
    file_id = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    //Create dataspace.  Setting maximum size to NULL sets the maximum
    //size to be the current size.
    space_id = H5Screate_simple (2, dims, NULL);

    //Create the dataset creation property list, set the layout to compact.
    property_id = H5Pcreate (H5P_DATASET_CREATE);
    status = H5Pset_layout (property_id, H5D_CONTIGUOUS);

    // Create the dataset. 
    dataset_id = H5Dcreate (file_id, "x", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, property_id, H5P_DEFAULT);
   
    //Write the data to the dataset.
    status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    status = H5Sclose(space_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);
    status = H5Pclose(property_id);

}

};

void usage(FILE *out, char *argv0) {
    fprintf(out, "Usage: %s <testfile>\n", argv0);
}

int main(int argc, char ** argv)
{
    if(argc < 2) {
        usage(stderr, argv[0]);
        exit(1);
    }

    for(int arg = 1; arg < argc; arg++) {
        std::string filename = argv[arg];
        std::string prefix("output.");
        double * data = 0;
        int m_rows, m_cols = 0;
        
        double start_time, stop_time;
        
        // Read the file
        get_seconds(&start_time);
        util::read_hdf(filename, data, m_rows, m_cols);
        get_seconds(&stop_time);
        std::cout << filename << "," << stop_time - start_time << std::endl;
        filename = prefix+filename;
        
        // Write the file
        get_seconds(&start_time);
        util::write_hdf(filename, data, m_rows, m_cols);
        get_seconds(&stop_time);
        std::cout << filename << "," << stop_time - start_time << std::endl;
    }


    return 0;
}
