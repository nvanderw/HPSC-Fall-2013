#include <iostream>
#include <fstream>
#include <stdlib.h>   
#include <iomanip>  
#include "hdf5.h"

int main(int argc, char ** argv)
{

    if(argc < 2){
        std::cout << "please specify a file" << std::endl;
        return 0;
    }
    std::string filename = argv[1];

    hsize_t offset[2] = {2,3};
    if (argc == 4){
        offset[0] = atoi(argv[2]);
        offset[1] = atoi(argv[3]);
    }

    int m_rows = 4;
    int m_cols = 4; 

    herr_t  status;
    hid_t file = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t dset = H5Dopen (file, "DATASET", H5P_DEFAULT);

    hid_t dcpl = H5Dget_create_plist (dset);
    H5D_layout_t layout = H5Pget_layout (dcpl);

    hid_t space = H5Dget_space (dset);
    hsize_t dims[2], mdims[2];
    status = H5Sget_simple_extent_dims(space,dims,mdims);

    int dim1 = dims[0];
    int dim2 = dims[1];

    const int length = m_rows*m_cols;
    double * data = new double[length];

    // Select from the file

    hsize_t count[2] = {};
    status = H5Sselect_hyperslab (space, H5S_SELECT_SET, offset, NULL, count, NULL);

    // Create memory space
    hsize_t dim_out[2] = {};
    hsize_t offset_out[2] = {};
    hsize_t count_out[2] = {};

    hid_t memspace = H5Screate_simple(2,dim_out,NULL);   
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);

    status = H5Dread (dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, data);
    
    // Print to screen
    for(int r=0; r<m_rows; ++r){
        for(int c=0; c<m_cols; ++c){
            std::cout << std::setfill (' ') << std::setw (8);
            std::cout << data[c + r*m_cols];
        }
        std::cout << std::endl;
    }
             
    delete [] data;

    status = H5Pclose (dcpl);
    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Fclose (file);


    return 0;
}