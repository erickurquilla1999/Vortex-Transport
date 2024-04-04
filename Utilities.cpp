#include <iostream>
#include <H5Cpp.h>

#include "Utilities.H"

void write_hdf5_dataset(const std::string& filename, const std::string& dataset_name, const float* data, int size) {
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