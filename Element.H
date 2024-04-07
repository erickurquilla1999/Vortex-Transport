#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>
#include <string>

#include "Meshgeneration.H"
#include "Parameters.H"

class Element {
public:

    // Default constructor
    Element();
 
    // Element default constructor with arguments
    Element(const int& ele_num, const mesh& mesh_info, const std::vector<std::vector<double>>& nods_ref_spa, const parameters& parms);

    // Element properties
    double time; // storage time in the dimulation for each element
    int number; // element number
    int type; // element type 0 of squere angle is down and 1 if up
    int right_element; // elemnt to the right
    int left_element; // elemnt to the left
    int vertical_element; // elemnt in the vertical direction
    std::vector<std::vector<double>> vertices_coords_phys_space; // coordinates of the verctices of the element in physical space
    std::vector<std::vector<double>> nods_coords_refe_space; // coordinates of the interior nodes of the element in reference space
    std::vector<std::vector<double>> nods_coords_phys_space; // coordinates of the interior nodes of the element in physical space
    std::vector<std::vector<double>> hidrodynamics_vector_u; // containg the hydrodinamics vecto u for each interior node
    std::vector<std::vector<std::vector<double>>> hidrodynamics_vector_f; // containg the hydrodinamics vecto f for each interior node

    // // Element functions
    // void initialize_hydrodinamics(int id); // initialize the hydronimics quantities u and f
    // void save_data(int id); // save element data

};

#endif // ELEMENT_H