#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include <string>

// write an hdf5 with name filemane, the data set is dataset_name and the data must be a float array, size is the size of the float array
void write_hdf5_dataset(const std::string& filename, const std::string& dataset_name, const float* data, const int& size);

// write an text file with name filemane, the data is given by an string array named line, each element of the array will be a line
void writeToFile(const std::string& filename, const std::vector<std::string>& lines);

#endif // UTILITIES_H