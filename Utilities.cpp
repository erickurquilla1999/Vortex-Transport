#include <iostream>
#include <H5Cpp.h>
#include <fstream>
#include <string>
#include <vector>

#include "Utilities.H"

// write an hdf5 with name filemane, the data set is dataset_name and the data must be a float array, size is the size of the float array
void write_hdf5_dataset(const std::string& filename, const std::string& dataset_name, const float* data, const int& size) {
    // Open or create the HDF5 file
    H5::H5File file(filename, H5F_ACC_TRUNC);

    // Create a dataspace for the dataset
    hsize_t dims[1] = {static_cast<hsize_t>(size)};
    H5::DataSpace dataspace(1, dims);

    // Create a dataset within the file
    H5::DataSet dataset = file.createDataSet(dataset_name, H5::PredType::NATIVE_FLOAT, dataspace);

    // Write the data to the dataset
    dataset.write(data, H5::PredType::NATIVE_FLOAT);

    // Close the dataset and file
    dataset.close();
    file.close();
}

// write an text file with name filemane, the data is given by an string array named line, each element of the array will be a line
void writeToFile(const std::string& filename, const std::vector<std::string>& lines) {
    // Open the file for writing
    std::ofstream outfile(filename);

    // Check if the file is open
    if (!outfile.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Write each line to the file
    for (const std::string& line : lines) {
        outfile << line << std::endl;
    }

    // Close the file
    outfile.close();
}
