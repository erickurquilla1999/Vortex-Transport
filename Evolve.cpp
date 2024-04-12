#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "Element.H"
#include "Evolve.H"
#include "Lagrangebasis.H"

Evolve_element::Evolve_element() {}

Evolve_element::Evolve_element(Element* this_elem, Element* right_elem, Element* left_elem, Element* vertical_elem, const std::vector<std::vector<double>>& gau_int_l, const int& p_parms):
    
    // Initialize Evolve_element properties
    p(p_parms),                               // lagrange polinomial order    
    this_element(this_elem),            // this is the element to be evolved in time
    right_element(right_elem),          // this the element to the right to the element to be evolved in time
    left_element(left_elem),            // this the element to the left to the element to be evolved in time
    vertival_element(vertical_elem),    // this the element in the vertical direction to the element to be evolved in time
    gau_integ_line(gau_int_l),           // contains gauss quadrature information for line integral, first index in the gauss quadrature coordinate number, second index runs from 0 to 1. 0 is xi coordinate. 1 is weight.

    // this vector store the values of the lagrange polinomial in evaluated in the quadrature points for side 1, 2 and 3
    // first index runs over interior nodes number, that is, the lagrange polinimial that is one on this node
    // second item runs over the evaluation of the lagrange poliniam in the quadrature points for side 1, 2 and 3
    phi_in_quadrature_points_side_1( ( this->p + 1 ) * ( this->p + 2 ) / 2 , std::vector<double>( this->gau_integ_line.size() ) ),
    phi_in_quadrature_points_side_2( ( this->p + 1 ) * ( this->p + 2 ) / 2 , std::vector<double>( this->gau_integ_line.size() ) ),
    phi_in_quadrature_points_side_3( ( this->p + 1 ) * ( this->p + 2 ) / 2 , std::vector<double>( this->gau_integ_line.size() ) ),

    // to following vectors store the interpolation of the hidrodynamic vector u from the interior nodes of the elements to the quadrature points for the line integrals
    // first index runs over the quadrature points number, second index runs over 0 and 3 for the hidrodynamics quantities
    // plus for current element
    u_plus_side_1( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    u_plus_side_2( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    u_plus_side_3( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    // minus for elements in the boundaries
    u_minus_side_1( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    u_minus_side_2( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    u_minus_side_3( this->gau_integ_line.size() , std::vector<double>( 4 ) )

    {

    // get size, that is the number of quadrature point for line integration
    int size = this->gau_integ_line.size();

    // evaluate the lagrange polinomial in the quadrature points for line integral
    for (int i = 0; i < size; ++i) {

        // compute the coordinates in the 2d reference space of the quadrature point for the line integral for side 1
        std::vector<double> xi_eta_gauss_side_1(2);
        xi_eta_gauss_side_1[0] = gau_integ_line[i][0];
        xi_eta_gauss_side_1[1] = 0.0;
        
        // compute the coordinates in the 2d reference space of the quadrature point for the line integral for side 2
        std::vector<double> xi_eta_gauss_side_2(2);
        xi_eta_gauss_side_2[0] = gau_integ_line[i][0];
        xi_eta_gauss_side_2[1] = gau_integ_line[i][0];

        // compute the coordinates in the 2d reference space of the quadrature point for the line integral for side 3
        std::vector<double> xi_eta_gauss_side_3(2);
        xi_eta_gauss_side_3[0] = 0.0;
        xi_eta_gauss_side_3[1] = 1.0 - gau_integ_line[i][0];

        // evaluate the lagrange polinomial in the quadrature points for line integral
        std::vector<double> phi_in_xi_eta_gauss_side_1 = lagrange_basis_reference_space( this->p , xi_eta_gauss_side_1 ); 
        std::vector<double> phi_in_xi_eta_gauss_side_2 = lagrange_basis_reference_space( this->p , xi_eta_gauss_side_2 ); 
        std::vector<double> phi_in_xi_eta_gauss_side_3 = lagrange_basis_reference_space( this->p , xi_eta_gauss_side_3 ); 

        // save the phi(xi,eta) in the class variables
        for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {
            // this vector store the values of the lagrange polinomial in evaluated in the quadrature points
            // first index runs over interior nodes number, that is, the lagrange polinimial that is one on this node
            // second item runs over the evaluation of the lagrange poliniam in the quadrature points
            this->phi_in_quadrature_points_side_1[j][i] = phi_in_xi_eta_gauss_side_1[j];
            this->phi_in_quadrature_points_side_2[j][i] = phi_in_xi_eta_gauss_side_2[j];
            this->phi_in_quadrature_points_side_3[j][i] = phi_in_xi_eta_gauss_side_3[j];
        }
    }
}

// this function compute the numerical flux on the element boundaries, side 1, 2 and 3. 
void Evolve_element::compute_numerical_flux(){

    // get size, that is the number of quadrature point for line integration
    int size = this->gau_integ_line.size();

    // interpolate the u from the nodes to the gauss quadrature for side 1, 2 and 3

    // this loop runs over gauss quadrature points
    for (int i = 0; i < size; ++i) {

        // initialize all the new state vector at the quadrature point in zero to start interpolation
        // this loops runs over hidrodynamic index
        for (int k = 0; k < 3; ++k) {
            // initialize all the new state vector at the quadrature point in zero to start interpolation
            this->u_plus_side_1[i][k] = 0.0;
            this->u_plus_side_2[i][k] = 0.0;
            this->u_plus_side_3[i][k] = 0.0;
            // initialize all the new state vector at the quadrature point in zero to start interpolation
            this->u_minus_side_1[i][k] = 0.0;
            this->u_minus_side_2[i][k] = 0.0;
            this->u_minus_side_3[i][k] = 0.0;
        }

        // this loops runs over hidrodynamic index
        for (int m = 0; m < 3; ++m) {
            // this loop runs over interior nodes
            for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {    
                // do interpolation u = sum phi(x,y) * U 
                // this element
                this->u_plus_side_1[i][m] += this->phi_in_quadrature_points_side_1[j][i] * this->this_element->hidrodynamics_vector_u[j][m];
                this->u_plus_side_2[i][m] += this->phi_in_quadrature_points_side_2[j][i] * this->this_element->hidrodynamics_vector_u[j][m];
                this->u_plus_side_3[i][m] += this->phi_in_quadrature_points_side_3[j][i] * this->this_element->hidrodynamics_vector_u[j][m];
                // boundary elements
                this->u_minus_side_1[i][m] += this->phi_in_quadrature_points_side_1[j][i] * this->vertival_element->hidrodynamics_vector_u[j][m];
                this->u_minus_side_2[i][m] += this->phi_in_quadrature_points_side_2[j][i] * this->vertival_element->hidrodynamics_vector_u[j][m];
                this->u_minus_side_3[i][m] += this->phi_in_quadrature_points_side_3[j][i] * this->vertival_element->hidrodynamics_vector_u[j][m];
            }
        }
    }
}
