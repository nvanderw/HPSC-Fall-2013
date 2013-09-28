#include <iostream>
#include <fstream>
#include <stdlib.h>

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
    fprintf(out, "Usage: %s <testfile>\n", argv0);
    exit(1);
}

int main(int argc, char ** argv)
{
    if(argc < 2) {
        usage(stderr, argv[0]);
    }

    std::string filename = argv[1];
    std::string prefix("output.");
    double * data = 0;
    int m_rows, m_cols = 0;
    
    // Read the file
    util::read_text(filename, data, m_rows, m_cols);
    filename = prefix+filename;
    
    // Write the file
    util::write_text(filename, data, m_rows, m_cols);


    return 0;
}
