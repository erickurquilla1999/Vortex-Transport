#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>

// write an hdf5 with name filemane, the data set is dataset_name and the data must be a float array, size is the size of the float array
void write_hdf5_dataset(const std::string& filename, const std::string& dataset_name, const float* data, const int& size);

#endif // UTILITIES_H