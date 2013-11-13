#include <iostream>
#include <fstream>
#include <iomanip>  
#include "hdf5.h"

int main(int argc, char ** argv)
{

    if(argc < 2){
        std::cout << "please specify a file" << std::endl;
        return 0;
    }

    std::string filename = argv[1];
    
    herr_t status;

    //Create a new file using the default properties.
    hid_t file_id = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t dataset_id = H5Dopen(file_id, "DATASET", H5P_DEFAULT);
    hid_t space_id = H5Dget_space(dataset_id);
    
    int length = H5Sget_simple_extent_npoints(space_id);

    hsize_t  dims[2];
    hsize_t  mdims[2];
    status = H5Sget_simple_extent_dims(space_id,dims,mdims);
    
    int m_rows = dims[0];
    int m_cols = dims[1];
    
    double * data = new double[m_rows*m_cols];
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    status = H5Sclose(space_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);

    // Print to screen
    for(int r=0; r<m_rows; ++r){
        for(int c=0; c<m_cols; ++c){
            std::cout << std::setfill (' ') << std::setw (8);
            std::cout << data[c + r*m_cols];
        }
        std::cout << std::endl;
    }

    delete [] data;

    return 0;
}