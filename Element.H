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

    std::vector<std::vector<double>> jacobian; // jacobian between transformation from reference space to physical space d vec{x} / d vec{xi} = [ [ x2 - x1 , x3 - x1 ] , [ y2 - y1 , y3 - y1 ] ]
    double determinant_jacobian; // determinant of jacobian
    std::vector<std::vector<double>> inverse_jacobian; // jacobian between transformation from physical space to reference space d vec{xi} / d vec{x} = ( 1 / det( J ) ) * [ [ y3 - y1 , x1 - x3 ] , [ y1 - y2 , x2 - x1 ] ]
    double determinant_inverse_jacobian; // determinant of the inverse of the jacobian

    std::vector<std::vector<double>> inverse_mass_matrix_physical_space; // mass_ij = int in T phi_i phi_j dT . size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2
    std::vector<std::vector<std::vector<double>>> stiffness_matrix_physical_space; // S_ij = integral in T of ( Nabla phi_i ) phi_j dT . form :  hat{e}_x * matrix[ ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ] + hat{e}_y * matrix[ ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ] 
    
    // Element functions
    void initialize_hydrodinamics(const parameters& parms); // initialize the hydronimics quantities u and f
    void write_data(const parameters& parms, const int& step_num); // write element data in output directory

};

#endif // ELEMENT_H