#ifndef MESHGENERATION_H
#define MESHGENERATION_H
#include <vector>

#include "Parameters.H"

struct mesh{
    std::vector<int> element_number;
    std::vector<std::vector<std::vector<double>>> element_coordinates;
    std::vector<std::vector<int>> elements_at_boundary;
    std::vector<int> element_type;
};

// generate mesh
mesh generate_mesh(const parameters& parms);

std::vector<std::vector<double>> generate_nodes_reference_space(const parameters& parms);

#endif // MESHGENERATION_H