#ifndef MESHGENERATION_H
#define MESHGENERATION_H
#include <vector>

#include "Parameters.H"

struct mesh{
    std::vector<int> element_number;
    std::vector<std::vector<std::vector<double>>> element_coordinates;
    std::vector<std::vector<int>> elements_at_boundary;
    std::vector<int> element_type;
    std::vector<std::vector<double>> nodes_reference_space;
};

// generate mesh
mesh generate_mesh(const parameters& parms);

#endif // MESHGENERATION_H