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
    vertical_element(vertical_elem),    // this the element in the vertical direction to the element to be evolved in time
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

    // the following vectors store the interpolation of the hidrodynamic vector u from the interior nodes of the elements to the quadrature points for the line integrals
    // first index runs over the quadrature points number, second index runs over 0 and 3 for the hidrodynamics quantities
    // plus for current element
    U_plus_side_1( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    U_plus_side_2( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    U_plus_side_3( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    // minus for elements in the boundaries
    U_minus_side_1( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    U_minus_side_2( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    U_minus_side_3( this->gau_integ_line.size() , std::vector<double>( 4 ) ),

    // contains the numerical flux (roe flux) at the quadrature points for sides 1, 2 and 3
    // first index runs over the gauss quadrature points number, second index runs over 0 and 3 for the hidrodynamics quantities
    numerical_flux_side_1( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    numerical_flux_side_2( this->gau_integ_line.size() , std::vector<double>( 4 ) ),
    numerical_flux_side_3( this->gau_integ_line.size() , std::vector<double>( 4 ) ),

    // Discontinuos Galerkin method variables
    DG_numerical_flux_integration( ( this->p + 1 ) * ( this->p + 2 ) / 2 ,  std::vector<double>( 4 ) ), // DG vector that results from the integration of the numerical flux ( integral phi_i hat_{F} dl ). First index runs over interior nodes. Second index runs between 0 and 3 and represend hidrodynamics variables.
    DG_stiffness_vector( ( this->p + 1 ) * ( this->p + 2 ) / 2 ,  std::vector<double>( 4 ) ), // Stiffness vector: DG vector that results from area integral over element of nabla phi_i dot F dOmega. First index runs over interior nodes. Second index runs between 0 and 3 and represend hidrodynamics variables.
    DG_residual_vector( ( this->p + 1 ) * ( this->p + 2 ) / 2 ,  std::vector<double>( 4 ) ), // residual vector: DG vector that results from stiffness vector minus vector result of the numerical flux integration. First index runs over interior nodes. Second index runs between 0 and 3 and represend hidrodynamics variables.
    DG_time_derivative_U( ( this->p + 1 ) * ( this->p + 2 ) / 2 ,  std::vector<double>( 4 ) ) // Time_derivative of state vector U(t): DG vector that results from dot{U}_{i} = M^{-1}_{ij} R_{j} ---> mass matrix inverse times the residual vector. First index runs over interior nodes. Second index runs between 0 and 3 and represend hidrodynamics variables.

    {

}

// this function evaluate lagrande polinomials in the gauss quadrature at the element boundaries
void Evolve_element::evaluate_basis_in_quadrature_poits(){

    // number of quadrature points
    int number_quadrature_points = this->gau_integ_line.size(); 

    std::vector<double> quadrature_point_reference_space(2);
    std::vector<double> plus_phi_side_1( ( this->p + 1 ) * ( this->p + 2 ) / 2 );
    std::vector<double> plus_phi_side_2( ( this->p + 1 ) * ( this->p + 2 ) / 2 );
    std::vector<double> plus_phi_side_3( ( this->p + 1 ) * ( this->p + 2 ) / 2 );
    std::vector<double> minus_phi_side_1( ( this->p + 1 ) * ( this->p + 2 ) / 2 );
    std::vector<double> minus_phi_side_2( ( this->p + 1 ) * ( this->p + 2 ) / 2 );
    std::vector<double> minus_phi_side_3( ( this->p + 1 ) * ( this->p + 2 ) / 2 );

    // loop over quadrature points
    for (int i = 0; i < number_quadrature_points; ++i) {

        // side 1
        //                                 (                    sigma  , 0.0 )
        quadrature_point_reference_space = { this->gau_integ_line[i][0], 0.0 };
        plus_phi_side_1 = lagrange_basis_reference_space(this->p, quadrature_point_reference_space);
        //                                 (                    1.0 - sigma  , 0.0 )
        quadrature_point_reference_space = { 1.0 - this->gau_integ_line[i][0], 0.0 };
        minus_phi_side_1 = lagrange_basis_reference_space(this->p, quadrature_point_reference_space);

        // side 2
        //                                 (                     1.0 - sigma , sigma )
        quadrature_point_reference_space = { 1.0 - this->gau_integ_line[i][0], this->gau_integ_line[i][0] };
        plus_phi_side_2 = lagrange_basis_reference_space(this->p, quadrature_point_reference_space);
        //                                 (                    sigma  , 1.0 - sigma )
        quadrature_point_reference_space = { this->gau_integ_line[i][0], 1.0 - this->gau_integ_line[i][0] };
        minus_phi_side_2 = lagrange_basis_reference_space(this->p, quadrature_point_reference_space);

        // side 3
        //                                 ( 0.0  , 1.0 - sigma )
        quadrature_point_reference_space = { 0.0 , 1.0 - this->gau_integ_line[i][0] };
        plus_phi_side_3 = lagrange_basis_reference_space(this->p, quadrature_point_reference_space);
        //                                 ( 0.0  , sigma )
        quadrature_point_reference_space = { 0.0 , this->gau_integ_line[i][0] };
        minus_phi_side_3 = lagrange_basis_reference_space(this->p, quadrature_point_reference_space);

        for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {
            
            plus_phi_in_quadrature_points_side_1[j][i]  = plus_phi_side_1[j];
            plus_phi_in_quadrature_points_side_2[j][i]  = plus_phi_side_2[j];
            plus_phi_in_quadrature_points_side_3[j][i]  = plus_phi_side_3[j];
            minus_phi_in_quadrature_points_side_1[j][i] = minus_phi_side_1[j];
            minus_phi_in_quadrature_points_side_2[j][i] = minus_phi_side_2[j];
            minus_phi_in_quadrature_points_side_3[j][i] = minus_phi_side_3[j];
        
        }
    }
}

// this function compute the hidrodynamic vector U on the element boundaries, side 1, 2 and 3. 
void Evolve_element::compute_U_plus_minus(){

    // number of quadrature points
    int number_quadrature_points = this->gau_integ_line.size();

    // loop over quadrature points
    for (int i = 0; i < number_quadrature_points; ++i) {

        U_plus_side_1[i]  = { 0.0 , 0.0, 0.0, 0.0 };
        U_plus_side_2[i]  = { 0.0 , 0.0, 0.0, 0.0 };
        U_plus_side_3[i]  = { 0.0 , 0.0, 0.0, 0.0 };
        U_minus_side_1[i] = { 0.0 , 0.0, 0.0, 0.0 };
        U_minus_side_2[i] = { 0.0 , 0.0, 0.0, 0.0 };
        U_minus_side_3[i] = { 0.0 , 0.0, 0.0, 0.0 };

        // loop over hidrodynamics indices
        for (int k = 0; k < 4; ++k) {
            // loop over all interior nodes
            for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {
                U_plus_side_1[i][k] += plus_phi_in_quadrature_points_side_1[j][i] * this->this_element->hidrodynamics_vector_U[j][k];
                U_plus_side_2[i][k] += plus_phi_in_quadrature_points_side_2[j][i] * this->this_element->hidrodynamics_vector_U[j][k];
                U_plus_side_3[i][k] += plus_phi_in_quadrature_points_side_3[j][i] * this->this_element->hidrodynamics_vector_U[j][k];
            }
        }

        if ( this->this_element->type == 0 ) {
            // loop over hidrodynamics indices            
            for (int k = 0; k < 4; ++k) {
                // loop over all interior nodes
                for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {
                    U_minus_side_1[i][k] += minus_phi_in_quadrature_points_side_1[j][i] * this->vertical_element->hidrodynamics_vector_U[j][k];
                    U_minus_side_2[i][k] += minus_phi_in_quadrature_points_side_2[j][i] * this->right_element->hidrodynamics_vector_U[j][k];
                    U_minus_side_3[i][k] += minus_phi_in_quadrature_points_side_3[j][i] * this->left_element->hidrodynamics_vector_U[j][k];

                }
            }
        } 
        
        if ( this->this_element->type == 1 ) {
            // loop over hidrodynamics indices
            for (int k = 0; k < 4; ++k) {
                // loop over all interior nodes
                for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {
                    U_minus_side_1[i][k] += minus_phi_in_quadrature_points_side_1[j][i] * this->vertical_element->hidrodynamics_vector_U[j][k];
                    U_minus_side_2[i][k] += minus_phi_in_quadrature_points_side_2[j][i] * this->left_element->hidrodynamics_vector_U[j][k];
                    U_minus_side_3[i][k] += minus_phi_in_quadrature_points_side_3[j][i] * this->right_element->hidrodynamics_vector_U[j][k];
                }
            }
        }
    }
}

// this function compute the numerical flux on the element boundaries, side 1, 2 and 3. 
void Evolve_element::compute_numerical_flux(){

    // number of quadrature points
    int number_quadrature_points = this->gau_integ_line.size();

    // loop over quadrature points
    for (int i = 0; i < number_quadrature_points; ++i) {
        numerical_flux_side_1[i]  = numerical_flux(U_plus_side_1[i], U_minus_side_1[i], this->this_element->units_vectors_perpendicular_to_element_boundary[0]);
        numerical_flux_side_2[i]  = numerical_flux(U_plus_side_2[i], U_minus_side_2[i], this->this_element->units_vectors_perpendicular_to_element_boundary[1]);
        numerical_flux_side_3[i]  = numerical_flux(U_plus_side_3[i], U_minus_side_3[i], this->this_element->units_vectors_perpendicular_to_element_boundary[2]);
    }
}

// this function create the DG vector that results from the integration of the numerical flux ( integral phi_i hat_{F} dl )
void Evolve_element::integrate_numerical_flux(){

    // number of quadrature points
    int number_quadrature_points = this->gau_integ_line.size();

    // loop over all interior nodes
    for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {

        DG_numerical_flux_integration[i] = { 0.0 , 0.0, 0.0 , 0.0 };
        
        // loop over hidrodynamics indices
        for (int k = 0; k < 4; ++k) {
            // loop over quadrature points
            for (int j = 0; j < number_quadrature_points; ++j) {
                DG_numerical_flux_integration[i][k] += this->this_element->sides_lenght[0] * plus_phi_in_quadrature_points_side_1[i][j] * numerical_flux_side_1[j][k] * this->gau_integ_line[j][1];
                DG_numerical_flux_integration[i][k] += this->this_element->sides_lenght[1] * plus_phi_in_quadrature_points_side_2[i][j] * numerical_flux_side_2[j][k] * this->gau_integ_line[j][1];
                DG_numerical_flux_integration[i][k] += this->this_element->sides_lenght[2] * plus_phi_in_quadrature_points_side_3[i][j] * numerical_flux_side_3[j][k] * this->gau_integ_line[j][1];
            }
        }
    }
}

// compute stiffness vector (area integral over element of nabla phi_i dot F dOmega)
void Evolve_element::compute_stiffness_vector(){
    // loop over all the interior nodes of this element
    for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {
        // loop over hidrodynamics indices
        for (int k = 0; k < 4; ++k) {
            // initialize the DG_stiffness_vector[i][j] values to zero
            this->DG_stiffness_vector[i][k] = 0.0;
            // loop over all the interior nodes of this element    
            for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {
                //             S_i              +=                             S^{x}_{ij}                       *               f^{x}_{j}                             +                             S^{y}_{ij}                       *               f^{y}_{j}
                this->DG_stiffness_vector[i][k] += this->this_element->stiffness_matrix_physical_space[0][i][j] * this->this_element->hidrodynamics_vector_F[j][0][k] + this->this_element->stiffness_matrix_physical_space[1][i][j] * this->this_element->hidrodynamics_vector_F[j][1][k];
                // Stiffness vector: DG vector that results from area integral over element of nabla phi_i dot F dOmega.
                // First index runs over interior nodes. 
                // Second index runs between 0 and 3 and represend hidrodynamics variables.
            }       
        }
    }

    
    int elem_number = 0;
    if ( this->this_element->number == elem_number ){
        for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {
            for (int k = 0; k < 4; ++k) {
                std::cout << " DG_stiffness_vector["<< i <<"]["<< k <<"] = " << this->DG_stiffness_vector[i][k] << std::endl;
            }    
        }
    }

}

// compute  residual vector: DG vector that results from stiffness vector minus vector result of the numerical flux integration
void Evolve_element::compute_residual_vector(){
    // compute DG_residual_vector[i][j]
    // loop over all the interior nodes of this element
    for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {
        // loop over hidrodynamics indices
        for (int k = 0; k < 4; ++k) {
            //           R_i               =               S_i               -            int num F_i
            // this->DG_residual_vector[i][k] =  (this->DG_stiffness_vector[i][k] - this->DG_numerical_flux_integration[i][k] < 1.0e-8 ? 0 : this->DG_stiffness_vector[i][k] - this->DG_numerical_flux_integration[i][k]); 
            this->DG_residual_vector[i][k] = this->DG_stiffness_vector[i][k] - this->DG_numerical_flux_integration[i][k]; 
        }       
    }

    int elem_number = 0;
    if ( this->this_element->number == elem_number ){
        for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {
            for (int k = 0; k < 4; ++k) {
                std::cout << " DG_residual_vector["<< i <<"]["<< k <<"] = " << this->DG_residual_vector[i][k] << std::endl;
            }    
        }
    }
}

// compute residial vector stiffness vector minus vector result of the numerical flux integration
void Evolve_element::compute_time_derivative_U(){
    // compute DG_time_derivative_U[i][j]
    // loop over all the interior nodes of this element
    for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {
        // loop over hidrodynamics indices
        for (int k = 0; k < 4; ++k) {
            // initialize the DG_time_derivative_U[i][j] values to zero
            this->DG_time_derivative_U[i][k] = 0.0; 
            // loop over all the interior nodes of this element
            for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {
                //         dU_i/d_t              = sum_j                M^{-1}_{ij}                             *            R_{j}  
                this->DG_time_derivative_U[i][k] += this->this_element->inverse_mass_matrix_physical_space[i][j] * this->DG_residual_vector[j][k]; 
            }
        }       
    }

    int elem_number = 0;
    if ( this->this_element->number == elem_number ){
        for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {
            for (int k = 0; k < 4; ++k) {
                std::cout << " DG_time_derivative_U["<< i <<"]["<< k <<"] = " << this->DG_time_derivative_U[i][k] << std::endl;
            }    
        }
    }
}

// compute new vectors U and F
void Evolve_element::compute_new_U_and_F(double& time_step){

    double rho, u, v, E, H, p;
    double gamma = 1.4;

    // Euler time stepping
    // loop over all the interior nodes of this element
    for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {
        // loop over hidrodynamics indices
        for (int k = 0; k < 4; ++k) {
            this->this_element->hidrodynamics_vector_U[i][k] += time_step * this->DG_time_derivative_U[i][k]; 
        }

        rho = this->this_element->hidrodynamics_vector_U[i][0]; // density
        u   = this->this_element->hidrodynamics_vector_U[i][1] / rho; // horizontal velocity
        v   = this->this_element->hidrodynamics_vector_U[i][2] / rho; // vertical velocity
        E   = this->this_element->hidrodynamics_vector_U[i][3] / rho; // energy
        p   = rho * ( gamma - 1 ) * ( E - ( pow( u , 2) + pow( v , 2) ) / 2 ); // pressure
        H   = E + p / rho; // Entalpy

        // initialize the hidrodinamic vector f, x component  
        this->this_element->hidrodynamics_vector_F[i][0][0] = rho * u;
        this->this_element->hidrodynamics_vector_F[i][0][1] = rho * pow( u , 2 ) + p;
        this->this_element->hidrodynamics_vector_F[i][0][2] = rho * u * v;
        this->this_element->hidrodynamics_vector_F[i][0][3] = rho * u * H;

        // initialize the hidrodinamic vector f, y component  
        this->this_element->hidrodynamics_vector_F[i][1][0] = rho * v;
        this->this_element->hidrodynamics_vector_F[i][1][1] = rho * u * v;
        this->this_element->hidrodynamics_vector_F[i][1][2] = rho * pow( v , 2 ) + p;
        this->this_element->hidrodynamics_vector_F[i][1][3] = rho * v * H;
    }

    // advance in time 
    this->this_element->time += time_step;
}