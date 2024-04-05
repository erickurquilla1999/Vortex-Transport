#include <iostream>

#include "Parameters.H"
#include "Utilities.H"
#include "Meshgeneration.H"

int main(int argc, char* argv[]) {
    
    // read simulation paramaters
    parameters parms = read_input_files(argc, argv);
    
    float mesh = generate_mesh(parms);

    // // Example usage
    // std::string filename = "example.h5";
    // std::string dataset_name = "data";
    // float data[] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
    // int size = sizeof(data) / sizeof(data[0]);
    // write_hdf5_dataset(filename, dataset_name, data, size);

    return 0;
}









