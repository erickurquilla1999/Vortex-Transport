#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>

void write_hdf5_dataset(const std::string& filename, const std::string& dataset_name, const float* data, int size);

#endif // UTILITIES_H