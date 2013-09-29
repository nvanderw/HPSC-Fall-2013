#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "timing.h"

namespace util
{
    void read_text(const std::string &filename, 
                   double* &data, 
                   int &m_rows, 
                   int &m_cols){

        std::ifstream input_file(filename.c_str());
        input_file >> m_rows;
        input_file >> m_cols;

        data = new double[m_rows*m_cols];

        for(int r=0; r<m_rows; ++r)
        {
            for(int c=0; c<m_cols; ++c)
            {
                const int index = r*m_cols + c;
                input_file >> data[index];
            }
        }

    }

    void write_text(const std::string &filename, 
                          double* &data, 
                          int &m_rows, 
                          int &m_cols){

        std::ofstream output_file(filename.c_str());
        
        output_file << m_rows << " " << m_cols << std::endl;
        
        for(int r=0; r<m_rows; ++r)
        {
            for(int c=0; c<m_cols; ++c)
            {
                const int index = r*m_cols + c;
                output_file << data[index] << " ";
            }
            output_file << std::endl;
        }
        output_file << std::endl;
        output_file.close();

    }

};

void usage(FILE *out, char *argv0) {
    fprintf(out, "Usage: %s <testfile1> <testfile2> ...\n", argv0);
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
        get_seconds(&start_time);
        
        // Read the file
        util::read_text(filename, data, m_rows, m_cols);
        get_seconds(&stop_time);
        std::cout << filename << "," << stop_time - start_time << std::endl;
        filename = prefix+filename;
        
        // Write the file
        get_seconds(&start_time);
        util::write_text(filename, data, m_rows, m_cols);
        get_seconds(&stop_time);
        std::cout << filename << "," << stop_time - start_time << std::endl;
    }


    return 0;
}
