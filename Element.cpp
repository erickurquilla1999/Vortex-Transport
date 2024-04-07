#include <vector>
#include <string>

#include "Element.H"

Element::Element() {}

Element::Element(const int& ele_num, const mesh& mesh_info,  const std::vector<std::vector<double>>& nods_ref_spa){

    // Initialize element properties    
    time = 0.0; // Simulation initial time
    number = ele_num; // Element number
    type = mesh_info.element_type[ele_num]; // Element type 0 of squere angle is down and 1 if up
    right_element = mesh_info.elements_at_boundary[ele_num][0]; // Elemnt to the right
    left_element = mesh_info.elements_at_boundary[ele_num][1]; // Elemnt to the left
    vertical_element = mesh_info.elements_at_boundary[ele_num][2]; // Elemnt in the vertical direction
    vertices_coords_phys_space = mesh_info.element_coordinates[ele_num]; // Coordinates of the verctices of the element in physical space. First index represent vertices of the element (run between 0 and 2): 0 is the vertice at the square angle, 1 and 2 are the other vertices going counter clock wise. Second index represent x:0 and y:1 position
    nods_coords_refe_space = nods_ref_spa; // Coordinates of the interior nodes of the element in reference space. The first index in the node number, the second inxed runs between 0 and 1. 0: xi and 1: eta.
}

