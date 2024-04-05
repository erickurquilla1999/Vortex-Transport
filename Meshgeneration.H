#ifndef MESHGENERATION_H
#define MESHGENERATION_H
#include <vector>

#include "Parameters.H"

struct mesh{
    std::vector<int> element_number;
    std::vector<std::vector<double>> grid_points;
    std::vector<std::vector<int>> element_to_node;
    std::vector<std::vector<int>> elements_at_boundary;
    std::vector<int> element_type;

};

// generate mesh
mesh generate_mesh(const parameters& parms);

#endif // MESHGENERATION_H