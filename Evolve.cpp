#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "Element.H"
#include "Evolve.H"
#include "Lagrangebasis.H"
#include "Numericalflux.H"

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
    // plus represent this element
    plus_phi_in_quadrature_points_side_1( ( this->p + 1 ) * ( this->p + 2 ) / 2 , std::vector<double>( this->gau_integ_line.size() ) ),
    plus_phi_in_quadrature_points_side_2( ( this->p + 1 ) * ( this->p + 2 ) / 2 , std::vector<double>( this->gau_integ_line.size() ) ),
    plus_phi_in_quadrature_points_side_3( ( this->p + 1 ) * ( this->p + 2 ) / 2 , std::vector<double>( this->gau_integ_line.size() ) ),
    // minus represent boundary element
    minus_phi_in_quadrature_points_side_1( ( this->p + 1 ) * ( this->p + 2 ) / 2 , std::vector<double>( this->gau_integ_line.size() ) ),
    minus_phi_in_quadrature_points_side_2( ( this->p + 1 ) * ( this->p + 2 ) / 2 , std::vector<double>( this->gau_integ_line.size() ) ),
    minus_phi_in_quadrature_points_side_3( ( this->p + 1 ) * ( this->p + 2 ) / 2 , std::vector<double>( this->gau_integ_line.size() ) ),

    // to following vectors store the interpolation of the hidrodynamic vector u from the interior nodes of the elements to the quadrature points for the line integrals
    // first index runs over the quadrature points number, second index runs over 0 and 3 for the hidrodynamics quantities
    // plus for current element
    u_plus_side_1( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    u_plus_side_2( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    u_plus_side_3( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    // minus for elements in the boundaries
    u_minus_side_1( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    u_minus_side_2( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    u_minus_side_3( this->gau_integ_line.size() , std::vector<double>( 4 ) ),

    // contains the numerical flux (roe flux) at the quadrature points for sides 1, 2 and 3
    // first index runs over the gauss quadrature points number, second index runs over 0 and 3 for the hidrodynamics quantities
    numerical_flux_side_1( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    numerical_flux_side_2( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    numerical_flux_side_3( this->gau_integ_line.size() , std::vector<double>( 4 ) )

    {

}

// this function compute the numerical flux on the element boundaries, side 1, 2 and 3. 
void Evolve_element::evaluate_basis_in_quadrature_poits(){

    // this is the number of quadrature point for line integration
    int num_gauss_quad_pnts = this->gau_integ_line.size();

    // evaluate the lagrange polinomial in the quadrature points for line integral
    for (int i = 0; i < num_gauss_quad_pnts; ++i) {

        std::vector<double> plus_xi_eta_gauss_side_1(2);
        std::vector<double> minus_xi_eta_gauss_side_1(2);
        plus_xi_eta_gauss_side_1[0]  = gau_integ_line[i][0];
        plus_xi_eta_gauss_side_1[1]  = 0.0;
        minus_xi_eta_gauss_side_1[0] = 1 - gau_integ_line[i][0];
        minus_xi_eta_gauss_side_1[1] = 0.0;

        std::vector<double> plus_xi_eta_gauss_side_2(2);
        std::vector<double> minus_xi_eta_gauss_side_2(2);
        plus_xi_eta_gauss_side_2[0]  = 1 - gau_integ_line[i][0];
        plus_xi_eta_gauss_side_2[1]  = gau_integ_line[i][0];
        minus_xi_eta_gauss_side_2[0] = gau_integ_line[i][0];
        minus_xi_eta_gauss_side_2[1] = 1 - gau_integ_line[i][0];

        std::vector<double> plus_xi_eta_gauss_side_3(2);
        std::vector<double> minus_xi_eta_gauss_side_3(2);
        plus_xi_eta_gauss_side_3[0]  = 0.0;
        plus_xi_eta_gauss_side_3[1]  = 1-gau_integ_line[i][0];
        minus_xi_eta_gauss_side_3[0] = 0.0;
        minus_xi_eta_gauss_side_3[1] = gau_integ_line[i][0];
    
        // evaluate the lagrange polinomial in the quadrature points for line integral
        std::vector<double> plus_phi_in_xi_eta_gauss_side_1 = lagrange_basis_reference_space( this->p , plus_xi_eta_gauss_side_1 ); 
        std::vector<double> plus_phi_in_xi_eta_gauss_side_2 = lagrange_basis_reference_space( this->p , plus_xi_eta_gauss_side_2 ); 
        std::vector<double> plus_phi_in_xi_eta_gauss_side_3 = lagrange_basis_reference_space( this->p , plus_xi_eta_gauss_side_3 ); 

        // evaluate the lagrange polinomial in the quadrature points for line integral
        std::vector<double> minus_phi_in_xi_eta_gauss_side_1 = lagrange_basis_reference_space( this->p , minus_xi_eta_gauss_side_1 ); 
        std::vector<double> minus_phi_in_xi_eta_gauss_side_2 = lagrange_basis_reference_space( this->p , minus_xi_eta_gauss_side_2 ); 
        std::vector<double> minus_phi_in_xi_eta_gauss_side_3 = lagrange_basis_reference_space( this->p , minus_xi_eta_gauss_side_3 ); 

        // save the phi(xi,eta) in the class variables
        for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {
            // this vector store the values of the lagrange polinomial in evaluated in the quadrature points
            // first index runs over interior nodes number, that is, the lagrange polinimial that is one on this node
            // second item runs over the evaluation of the lagrange poliniam in the quadrature points
            this->plus_phi_in_quadrature_points_side_1[j][i] = plus_phi_in_xi_eta_gauss_side_1[j];
            this->plus_phi_in_quadrature_points_side_2[j][i] = plus_phi_in_xi_eta_gauss_side_2[j];
            this->plus_phi_in_quadrature_points_side_3[j][i] = plus_phi_in_xi_eta_gauss_side_3[j];
            this->minus_phi_in_quadrature_points_side_1[j][i] = minus_phi_in_xi_eta_gauss_side_1[j];
            this->minus_phi_in_quadrature_points_side_2[j][i] = minus_phi_in_xi_eta_gauss_side_2[j];
            this->minus_phi_in_quadrature_points_side_3[j][i] = minus_phi_in_xi_eta_gauss_side_3[j];
        }
    }

    //ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
    // test to see if the coordinates in reference space for the boundary elements match the coordinates in physical space
    if (this->this_element->type == 0){
        // side 1 type 0
        std::cout << " " << std::endl;
        for (int i = 0; i < num_gauss_quad_pnts; ++i) {
            std::vector<double> plus_xi_eta_gauss_side_1(2);
            std::vector<double> minus_xi_eta_gauss_side_1(2);
            plus_xi_eta_gauss_side_1[0] = gau_integ_line[i][0];
            plus_xi_eta_gauss_side_1[1] = 0.0;
            minus_xi_eta_gauss_side_1[0] = 1 - gau_integ_line[i][0];
            minus_xi_eta_gauss_side_1[1] = 0.0;
            std::vector<double> r1 = reference_to_physical_space(plus_xi_eta_gauss_side_1, this->this_element->vertices_coords_phys_space);
            std::vector<double> r2 = reference_to_physical_space(minus_xi_eta_gauss_side_1, this->vertival_element->vertices_coords_phys_space);
            std::cout << "Side 1 : This element "<< this->this_element->number << " ( " << r1[0] << " , " << r1[1] << " ). Vertical element " << this->vertival_element->number << " ( " << r2[0] << " , " << r2[1] << " )." << std::endl;
        }
        // side 2 type 0
        std::cout << " " << std::endl;
        for (int i = 0; i < num_gauss_quad_pnts; ++i) {
            std::vector<double> plus_xi_eta_gauss_side_2(2);
            std::vector<double> minus_xi_eta_gauss_side_2(2);
            plus_xi_eta_gauss_side_2[0] = 1 - gau_integ_line[i][0];
            plus_xi_eta_gauss_side_2[1] = gau_integ_line[i][0];
            minus_xi_eta_gauss_side_2[0] = gau_integ_line[i][0];
            minus_xi_eta_gauss_side_2[1] = 1 - gau_integ_line[i][0];
            std::vector<double> r1 = reference_to_physical_space(plus_xi_eta_gauss_side_2, this->this_element->vertices_coords_phys_space);
            std::vector<double> r2 = reference_to_physical_space(minus_xi_eta_gauss_side_2, this->right_element->vertices_coords_phys_space);
            std::cout << "Side 2 : This element "<< this->this_element->number << " ( " << r1[0] << " , " << r1[1] << " ). Right element " << this->right_element->number << " ( " << r2[0] << " , " << r2[1] << " )." << std::endl;
        }
        // side 3 type 0
        std::cout << " " << std::endl;
        for (int i = 0; i < num_gauss_quad_pnts; ++i) {
            std::vector<double> plus_xi_eta_gauss_side_3(2);
            std::vector<double> minus_xi_eta_gauss_side_3(2);
            plus_xi_eta_gauss_side_3[0] = 0.0;
            plus_xi_eta_gauss_side_3[1] = 1-gau_integ_line[i][0];
            minus_xi_eta_gauss_side_3[0] = 0.0;
            minus_xi_eta_gauss_side_3[1] = gau_integ_line[i][0];
            std::vector<double> r1 = reference_to_physical_space(plus_xi_eta_gauss_side_3, this->this_element->vertices_coords_phys_space);
            std::vector<double> r2 = reference_to_physical_space(minus_xi_eta_gauss_side_3, this->left_element->vertices_coords_phys_space);
            std::cout << "Side 3 : This element "<< this->this_element->number << " ( " << r1[0] << " , " << r1[1] << " ). Left element " << this->left_element->number << " ( " << r2[0] << " , " << r2[1] << " )." << std::endl;
        }
    }
    else if (this->this_element->type == 1){
        // side 1 type 1
        std::cout << " " << std::endl;
        for (int i = 0; i < num_gauss_quad_pnts; ++i) {
            std::vector<double> plus_xi_eta_gauss_side_1(2);
            std::vector<double> minus_xi_eta_gauss_side_1(2);
            plus_xi_eta_gauss_side_1[0] = gau_integ_line[i][0];
            plus_xi_eta_gauss_side_1[1] = 0.0;
            minus_xi_eta_gauss_side_1[0] = 1 - gau_integ_line[i][0];
            minus_xi_eta_gauss_side_1[1] = 0.0;
            std::vector<double> r1 = reference_to_physical_space(plus_xi_eta_gauss_side_1, this->this_element->vertices_coords_phys_space);
            std::vector<double> r2 = reference_to_physical_space(minus_xi_eta_gauss_side_1, this->vertival_element->vertices_coords_phys_space);
            std::cout << "Side 1 : This element "<< this->this_element->number << " ( " << r1[0] << " , " << r1[1] << " ). Vertical element " << this->vertival_element->number << " ( " << r2[0] << " , " << r2[1] << " )." << std::endl;
        }
        // side 2 type 1
        std::cout << " " << std::endl;
        for (int i = 0; i < num_gauss_quad_pnts; ++i) {
            std::vector<double> plus_xi_eta_gauss_side_2(2);
            std::vector<double> minus_xi_eta_gauss_side_2(2);
            plus_xi_eta_gauss_side_2[0] = 1 - gau_integ_line[i][0];
            plus_xi_eta_gauss_side_2[1] = gau_integ_line[i][0];
            minus_xi_eta_gauss_side_2[0] = gau_integ_line[i][0];
            minus_xi_eta_gauss_side_2[1] = 1 - gau_integ_line[i][0];
            std::vector<double> r1 = reference_to_physical_space(plus_xi_eta_gauss_side_2, this->this_element->vertices_coords_phys_space);
            std::vector<double> r2 = reference_to_physical_space(minus_xi_eta_gauss_side_2, this->left_element->vertices_coords_phys_space);
            std::cout << "Side 2 : This element "<< this->this_element->number << " ( " << r1[0] << " , " << r1[1] << " ). Left element " << this->left_element->number << " ( " << r2[0] << " , " << r2[1] << " )." << std::endl;
        }
        // side 3 type 1
        std::cout << " " << std::endl;
        for (int i = 0; i < num_gauss_quad_pnts; ++i) {
            std::vector<double> plus_xi_eta_gauss_side_3(2);
            std::vector<double> minus_xi_eta_gauss_side_3(2);
            plus_xi_eta_gauss_side_3[0] = 0.0;
            plus_xi_eta_gauss_side_3[1] = 1-gau_integ_line[i][0];
            minus_xi_eta_gauss_side_3[0] = 0.0;
            minus_xi_eta_gauss_side_3[1] = gau_integ_line[i][0];
            std::vector<double> r1 = reference_to_physical_space(plus_xi_eta_gauss_side_3, this->this_element->vertices_coords_phys_space);
            std::vector<double> r2 = reference_to_physical_space(minus_xi_eta_gauss_side_3, this->right_element->vertices_coords_phys_space);
            std::cout << "Side 3 : This element "<< this->this_element->number << " ( " << r1[0] << " , " << r1[1] << " ). Right element " << this->right_element->number << " ( " << r2[0] << " , " << r2[1] << " )." << std::endl;
        }
    }else{
        std::cerr << " Error: Unsoported element type: it must be 0 or 1! " << std::endl;
        exit(EXIT_FAILURE);
    }
    //ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd

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
                this->u_plus_side_1[i][m] += this->plus_phi_in_quadrature_points_side_1[j][i] * this->this_element->hidrodynamics_vector_u[j][m];
                this->u_plus_side_2[i][m] += this->plus_phi_in_quadrature_points_side_2[j][i] * this->this_element->hidrodynamics_vector_u[j][m];
                this->u_plus_side_3[i][m] += this->plus_phi_in_quadrature_points_side_3[j][i] * this->this_element->hidrodynamics_vector_u[j][m];
                // boundary elements
                if (this->this_element->type == 0){
                    this->u_minus_side_1[i][m] += this->minus_phi_in_quadrature_points_side_1[j][i] * this->vertival_element->hidrodynamics_vector_u[j][m];
                    this->u_minus_side_2[i][m] += this->minus_phi_in_quadrature_points_side_2[j][i] * this->right_element->hidrodynamics_vector_u[j][m];
                    this->u_minus_side_3[i][m] += this->minus_phi_in_quadrature_points_side_3[j][i] * this->left_element->hidrodynamics_vector_u[j][m];
                }
                else if (this->this_element->type == 1){
                    this->u_minus_side_1[i][m] += this->minus_phi_in_quadrature_points_side_1[j][i] * this->vertival_element->hidrodynamics_vector_u[j][m];
                    this->u_minus_side_2[i][m] += this->minus_phi_in_quadrature_points_side_2[j][i] * this->left_element->hidrodynamics_vector_u[j][m];
                    this->u_minus_side_3[i][m] += this->minus_phi_in_quadrature_points_side_3[j][i] * this->right_element->hidrodynamics_vector_u[j][m];
                }else{
                    std::cerr << " Error: Unsoported element type: it must be 0 or 1! " << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    // compute the roe flux
    // this loop runs over gauss quadrature points
    for (int i = 0; i < size; ++i) {
        // initialize the numerical flux at quadrature point in zero
        // this loops runs over hidrodynamic index
        for (int k = 0; k < 3; ++k) {
            this->numerical_flux_side_1[i][k] = 0.0;
            this->numerical_flux_side_2[i][k] = 0.0;
            this->numerical_flux_side_3[i][k] = 0.0;
        }

        for (int a = 0; a < 4; ++a) {
            std::cout << "This element "<< this->this_element->number << " , hidro indx " << a << " : "<<this->u_plus_side_1[i][a] << " , " << this->u_minus_side_1[i][a] << std::endl;
        }

        std::vector<double> num_flux_side_1 = numerical_flux(this->u_plus_side_1[i], this->u_minus_side_1[i], this->this_element->units_vectors_perpendicular_to_element_boundary[0]);
        std::vector<double> num_flux_side_2 = numerical_flux(this->u_plus_side_2[i], this->u_minus_side_2[i], this->this_element->units_vectors_perpendicular_to_element_boundary[1]);
        std::vector<double> num_flux_side_3 = numerical_flux(this->u_plus_side_3[i], this->u_minus_side_3[i], this->this_element->units_vectors_perpendicular_to_element_boundary[2]);
        
        // save numerical flux in evolve_element object
        // this loops runs over hidrodynamic index
        for (int k = 0; k < 3; ++k) {
            this->numerical_flux_side_1[i][k] = num_flux_side_1[k];
            this->numerical_flux_side_2[i][k] = num_flux_side_2[k];
            this->numerical_flux_side_3[i][k] = num_flux_side_3[k];
        }
    }    
}