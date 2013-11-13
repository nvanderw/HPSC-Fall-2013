#include <iostream>
#include <fstream>
#include "stdlib.h"
#include "hdf5.h"

int main(int argc, char ** argv)
{

    if(argc < 2){
        std::cout << "please specify a file" << std::endl;
        return 0;
    }

    std::string filename = argv[1];

    int m_rows = 12;
    int m_cols = 8;
    if( argc == 4 )
    {
        m_rows = atoi(argv[2]);
        m_cols = atoi(argv[3]);
    }

    double * data = new double[m_rows*m_cols];

    // Fill with some data
    for(int r=0; r<m_rows; ++r)
        for(int c=0; c<m_cols; ++c)
            data[r + c*m_rows] = r + c*m_rows+1.0;

    herr_t status;
    hsize_t  dims[2] = {m_rows,m_cols};
   
    //Create a new file using the default properties.
    hid_t file_id = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    //Create dataspace.  Setting maximum size to NULL sets the maximum
    //size to be the current size.
    hid_t space_id = H5Screate_simple (2, dims, NULL);

    //Create the dataset creation property list, set the layout to compact.
    hid_t property_id = H5Pcreate (H5P_DATASET_CREATE);
    status = H5Pset_layout (property_id, H5D_CONTIGUOUS);

    // Create the dataset. 
    hid_t dataset_id = H5Dcreate (file_id, "DATASET", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, property_id, H5P_DEFAULT);
   
    //Write the data to the dataset.
    status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    status = H5Sclose(space_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);
    status = H5Pclose(property_id);

    return 0;
}