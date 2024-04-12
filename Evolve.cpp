#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "Element.H"
#include "Evolve.H"
#include "Lagrangebasis.H"

Evolve_element::Evolve_element() {}

Evolve_element::Evolve_element(Element* this_elem, Element* right_elem, Element* left_elem, Element* vertical_elem, const std::vector<std::vector<double>>& gau_int_l, const int& p):
    
    // Initialize Evolve_element properties    
    this_element(this_elem),            // this is the element to be evolved in time
    right_element(right_elem),          // this the element to the right to the element to be evolved in time
    left_element(left_elem),            // this the element to the left to the element to be evolved in time
    vertival_element(vertical_elem),    // this the element in the vertical direction to the element to be evolved in time
    gau_integ_line(gau_int_l)           // contains gauss quadrature information for line integral, first index in the gauss quadrature coordinate number, second index runs from 0 to 1. 0 is xi coordinate. 1 is weight.

    {

    //ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
    std::cout << "Element : " << this->this_element->hidrodynamics_vector_u[0][0] << std::endl;
    this->this_element->hidrodynamics_vector_u[0][0] = 9999999.9;
    //ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd

    // get size, that is the number of quadrature point for line integration
    int size = this->gau_integ_line.size();

    // this vector store the values of the lagrange polinomial in evaluated in the quadrature points for side 1, 2 and 3
    // first index runs over interior nodes number, that is, the lagrange polinimial that is one on this node
    // second item runs over the evaluation of the lagrange poliniam in the quadrature points for side 1, 2 and 3
    std::vector<std::vector<double>> phi_in_quadrature_points_side_1( ( p + 1 ) * ( p + 2 ) / 2 , std::vector<double>( size ) );
    std::vector<std::vector<double>> phi_in_quadrature_points_side_2( ( p + 1 ) * ( p + 2 ) / 2 , std::vector<double>( size ) );
    std::vector<std::vector<double>> phi_in_quadrature_points_side_3( ( p + 1 ) * ( p + 2 ) / 2 , std::vector<double>( size ) );

    // initialize counter
    int counter = 0;

    // evaluate the lagrange polinomial in the quadrature points for line integral
    for (int i = 0; i < size; ++i) {

        std::vector<double> xi_eta_gauss_side_1(2);
        xi_eta_gauss_side_1[0] = gau_integ_line[i][0];
        xi_eta_gauss_side_1[1] = 0.0;

        std::vector<double> xi_eta_gauss_side_2(2);
        xi_eta_gauss_side_2[0] = gau_integ_line[i][0];
        xi_eta_gauss_side_2[1] = gau_integ_line[i][0];

        std::vector<double> xi_eta_gauss_side_3(2);
        xi_eta_gauss_side_3[0] = 0.0;
        xi_eta_gauss_side_3[1] = 1.0 - gau_integ_line[i][0];

        std::vector<double> phi_in_xi_eta_gauss_side_1 = lagrange_basis_reference_space( p , xi_eta_gauss_side_1 ); 
        std::vector<double> phi_in_xi_eta_gauss_side_2 = lagrange_basis_reference_space( p , xi_eta_gauss_side_2 ); 
        std::vector<double> phi_in_xi_eta_gauss_side_3 = lagrange_basis_reference_space( p , xi_eta_gauss_side_3 ); 

        // evaluate the lagrange polinomial in the quadrature points
        for (int j = 0; j < ( p + 1 ) * ( p + 2 ) / 2; ++j) {
            phi_in_quadrature_points_side_1[j][counter] = phi_in_xi_eta_gauss_side_1[j];
            phi_in_quadrature_points_side_2[j][counter] = phi_in_xi_eta_gauss_side_2[j];
            phi_in_quadrature_points_side_3[j][counter] = phi_in_xi_eta_gauss_side_3[j];
        }
        counter++;
    }
}