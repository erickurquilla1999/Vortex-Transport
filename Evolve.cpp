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

    // this is the number of quadrature point for line integration
    int num_gauss_quad_pnts = this->gau_integ_line.size();

    // evaluate the lagrange polinomial in the quadrature points for line integral
    for (int i = 0; i < num_gauss_quad_pnts; ++i) {

        // Side 1. This is the first side going counterclockwise from the square angle vertex of the element
        std::vector<double> plus_xi_eta_gauss_side_1(2);  // Cordinates in the 2d reference space ( xi(xi') , eta(xi') ) of the gauss quadrature xi' that is inside [0, 1] for this element
        std::vector<double> minus_xi_eta_gauss_side_1(2); // Cordinates in the 2d reference space ( xi(xi') , eta(xi') ) of the gauss quadrature xi' that is inside [0, 1] for the in the boundary of side one
        plus_xi_eta_gauss_side_1[0]  = gau_integ_line[i][0];     // xi = xi'
        plus_xi_eta_gauss_side_1[1]  = 0.0;                      // eta = 0
        minus_xi_eta_gauss_side_1[0] = 1 - gau_integ_line[i][0]; // xi = 1 - xi'
        minus_xi_eta_gauss_side_1[1] = 0.0;                      // eta = 0

        // Side 2. This is the second side going counterclockwise from the square angle vertex of the element
        std::vector<double> plus_xi_eta_gauss_side_2(2);  // Cordinates in the 2d reference space ( xi(xi') , eta(xi') ) of the gauss quadrature xi' that is inside [0, 1] for this element
        std::vector<double> minus_xi_eta_gauss_side_2(2); // Cordinates in the 2d reference space ( xi(xi') , eta(xi') ) of the gauss quadrature xi' that is inside [0, 1] for the in the boundary of side one
        plus_xi_eta_gauss_side_2[0]  = 1 - gau_integ_line[i][0]; // xi = 1 - xi'
        plus_xi_eta_gauss_side_2[1]  = gau_integ_line[i][0];     // eta = xi'
        minus_xi_eta_gauss_side_2[0] = gau_integ_line[i][0];     // xi = xi'
        minus_xi_eta_gauss_side_2[1] = 1 - gau_integ_line[i][0]; // eta = 1 - xi'

        // Side 3. This is the third side going counterclockwise from the square angle vertex of the element
        std::vector<double> plus_xi_eta_gauss_side_3(2);  // Cordinates in the 2d reference space ( xi(xi') , eta(xi') ) of the gauss quadrature xi' that is inside [0, 1] for this element
        std::vector<double> minus_xi_eta_gauss_side_3(2); // Cordinates in the 2d reference space ( xi(xi') , eta(xi') ) of the gauss quadrature xi' that is inside [0, 1] for the in the boundary of side one
        plus_xi_eta_gauss_side_3[0]  = 0.0;                    // xi = 0
        plus_xi_eta_gauss_side_3[1]  = 1-gau_integ_line[i][0]; // eta = 1 - xi'
        minus_xi_eta_gauss_side_3[0] = 0.0;                    // xi = 0
        minus_xi_eta_gauss_side_3[1] = gau_integ_line[i][0];   // eta = xi'
    
        // evaluate the lagrange polinomial in the quadrature points for line integral for this element
        std::vector<double> plus_phi_in_xi_eta_gauss_side_1 = lagrange_basis_reference_space( this->p , plus_xi_eta_gauss_side_1 ); // side 1
        std::vector<double> plus_phi_in_xi_eta_gauss_side_2 = lagrange_basis_reference_space( this->p , plus_xi_eta_gauss_side_2 ); // side 2 
        std::vector<double> plus_phi_in_xi_eta_gauss_side_3 = lagrange_basis_reference_space( this->p , plus_xi_eta_gauss_side_3 ); // side 3 

        // evaluate the lagrange polinomial in the quadrature points for line integral for the bournady element
        std::vector<double> minus_phi_in_xi_eta_gauss_side_1 = lagrange_basis_reference_space( this->p , minus_xi_eta_gauss_side_1 ); // side 1
        std::vector<double> minus_phi_in_xi_eta_gauss_side_2 = lagrange_basis_reference_space( this->p , minus_xi_eta_gauss_side_2 ); // side 2
        std::vector<double> minus_phi_in_xi_eta_gauss_side_3 = lagrange_basis_reference_space( this->p , minus_xi_eta_gauss_side_3 ); // side 3 

        // save the phi(xi,eta) in the class variables
        for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {
            // this vector store the values of the lagrange polinomial in evaluated in the quadrature points
            // first index runs over interior nodes number, that is, the lagrange polinimial that is one on this node
            // second item runs over the evaluation of the lagrange poliniam in the quadrature points
            this->plus_phi_in_quadrature_points_side_1[j][i] = plus_phi_in_xi_eta_gauss_side_1[j];   // this element side 1
            this->plus_phi_in_quadrature_points_side_2[j][i] = plus_phi_in_xi_eta_gauss_side_2[j];   // this element side 2
            this->plus_phi_in_quadrature_points_side_3[j][i] = plus_phi_in_xi_eta_gauss_side_3[j];   // this element side 3
            this->minus_phi_in_quadrature_points_side_1[j][i] = minus_phi_in_xi_eta_gauss_side_1[j]; // boundary element side 1
            this->minus_phi_in_quadrature_points_side_2[j][i] = minus_phi_in_xi_eta_gauss_side_2[j]; // boundary element side 2
            this->minus_phi_in_quadrature_points_side_3[j][i] = minus_phi_in_xi_eta_gauss_side_3[j]; // boundary element side 3
        }
    }
}

// this function compute the numerical flux on the element boundaries, side 1, 2 and 3. 
void Evolve_element::compute_numerical_flux(){

    // get size, that is the number of quadrature point for line integration
    int num_gauss_quad_pnts = this->gau_integ_line.size();

    // interpolate the u from the nodes to the gauss quadrature for side 1, 2 and 3

    // this loop runs over gauss quadrature points
    for (int i = 0; i < num_gauss_quad_pnts; ++i) {

        // initialize all the new state vectors at the quadrature point in zero to start interpolation
        // this loops runs over hidrodynamic index
        for (int k = 0; k < 4; ++k) {
            // initialize all the new state vector at the quadrature point in zero to start interpolation
            this->u_plus_side_1[i][k] = 0.0; // this element side 1
            this->u_plus_side_2[i][k] = 0.0; // this element side 2
            this->u_plus_side_3[i][k] = 0.0; // this element side 3
            // initialize all the new state vector at the quadrature point in zero to start interpolation
            this->u_minus_side_1[i][k] = 0.0; // boundary element side 1
            this->u_minus_side_2[i][k] = 0.0; // boundary element side 2
            this->u_minus_side_3[i][k] = 0.0; // boundary element side 3
        }

        // this loops runs over hidrodynamic index
        for (int m = 0; m < 4; ++m) {
            // this loop runs over interior nodes
            for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {    
                // do interpolation u = sum phi(x,y) * U(t) 
                // this element
                this->u_plus_side_1[i][m] += this->plus_phi_in_quadrature_points_side_1[j][i] * this->this_element->hidrodynamics_vector_u[j][m]; // this element side 1
                this->u_plus_side_2[i][m] += this->plus_phi_in_quadrature_points_side_2[j][i] * this->this_element->hidrodynamics_vector_u[j][m]; // this element side 2
                this->u_plus_side_3[i][m] += this->plus_phi_in_quadrature_points_side_3[j][i] * this->this_element->hidrodynamics_vector_u[j][m]; // this element side 3
                // boundary elements
                // this->this_element->type contains information if the element has the square angle up or down.
                // 0: square angle is down
                // 1: square angle is up
                // This is thinking in a not perturbed mesh, but still valid for perturbed mesh.
                if (this->this_element->type == 0){
                    // if this->this_element->type == 0
                    // the element on the boundary of side 1 is the vertical element
                    // the element on the boundary of side 2 is the rigth element
                    // the element on the boundary of side 3 is the left element
                    this->u_minus_side_1[i][m] += this->minus_phi_in_quadrature_points_side_1[j][i] * this->vertival_element->hidrodynamics_vector_u[j][m];
                    this->u_minus_side_2[i][m] += this->minus_phi_in_quadrature_points_side_2[j][i] * this->right_element->hidrodynamics_vector_u[j][m];
                    this->u_minus_side_3[i][m] += this->minus_phi_in_quadrature_points_side_3[j][i] * this->left_element->hidrodynamics_vector_u[j][m];
                }
                else if (this->this_element->type == 1){
                    // if this->this_element->type == 1
                    // the element on the boundary of side 1 is the vertical element
                    // the element on the boundary of side 2 is the left element
                    // the element on the boundary of side 3 is the right element
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

    // compute the numerical flux: roe flux
    // this loop runs over gauss quadrature points
    for (int i = 0; i < num_gauss_quad_pnts; ++i) {
        // initialize the numerical flux at quadrature point in zero
        // this loops runs over hidrodynamic index
        for (int k = 0; k < 4; ++k) {
            this->numerical_flux_side_1[i][k] = 0.0; // side 1
            this->numerical_flux_side_2[i][k] = 0.0; // side 2
            this->numerical_flux_side_3[i][k] = 0.0; // side 3
        }
        
        // call the function numerical_flux in the Numericalflux.cpp file
        // this is a 4 element array just for hidrodinamic index
        std::vector<double> num_flux_side_1 = numerical_flux(this->u_plus_side_1[i], this->u_minus_side_1[i], this->this_element->units_vectors_perpendicular_to_element_boundary[0]); // side 1
        std::vector<double> num_flux_side_2 = numerical_flux(this->u_plus_side_2[i], this->u_minus_side_2[i], this->this_element->units_vectors_perpendicular_to_element_boundary[1]); // side 2
        std::vector<double> num_flux_side_3 = numerical_flux(this->u_plus_side_3[i], this->u_minus_side_3[i], this->this_element->units_vectors_perpendicular_to_element_boundary[2]); // side 3
        
        // save numerical flux in evolve_element object
        // this loops runs over hidrodynamic index
        for (int k = 0; k < 4; ++k) {
            // contains the numerical flux (roe flux) at the quadrature points for sides 1, 2 and 3
            // first index runs over the gauss quadrature points number, second index runs over 0 and 3 for the hidrodynamics quantities
            this->numerical_flux_side_1[i][k] = num_flux_side_1[k]; // side 1
            this->numerical_flux_side_2[i][k] = num_flux_side_2[k]; // side 2
            this->numerical_flux_side_3[i][k] = num_flux_side_3[k]; // side 3
        }
    }    
}

// this function create the DG vector that results from the integration of the numerical flux ( integral phi_i hat_{F} dl )
void Evolve_element::integrate_numerical_flux(){

    // get size, that is the number of quadrature point for line integration
    int num_gauss_quad_pnts = this->gau_integ_line.size();

    // loop over all the interior nodes of this element
    for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {    
        
        // contains integral over sides 1, 2 and 3
        std::vector<double> integral_side_1(4, 0.0);
        std::vector<double> integral_side_2(4, 0.0);
        std::vector<double> integral_side_3(4, 0.0);

        // loop over hidrodynamics index
        for (int k = 0; k < 4; ++k) {
            // loop over the quadrature point 
            for (int j = 0; j < num_gauss_quad_pnts; ++j) {
                //                 +=           w(xi')           * (                    phi(xi')                      *             \hat_{f}              )
                integral_side_1[k] += this->gau_integ_line[j][1] * ( this->plus_phi_in_quadrature_points_side_1[i][j] * this->numerical_flux_side_1[j][k] );
                integral_side_2[k] += this->gau_integ_line[j][1] * ( this->plus_phi_in_quadrature_points_side_2[i][j] * this->numerical_flux_side_2[j][k] );
                integral_side_3[k] += this->gau_integ_line[j][1] * ( this->plus_phi_in_quadrature_points_side_3[i][j] * this->numerical_flux_side_3[j][k] );
            }
        }

        // loop over hidrodynamics index
        for (int k = 0; k < 4; ++k) {
            // sum the integrals over the three side of the element
            this->DG_numerical_flux_integration[i][k] = this->this_element->sides_lenght[0] * integral_side_1[k] + this->this_element->sides_lenght[1] * integral_side_2[k] + this->this_element->sides_lenght[2] * integral_side_3[k];
            // DG_numerical_flux_integration is the DG vector that results from the integration of the numerical flux ( integral phi_i hat_{F} dl ). 
            // First index runs over interior nodes. 
            // Second index runs between 0 and 3 and represend hidrodynamics variables.
        }
    }
}

// compute stiffness vector (area integral over element of nabla phi_i dot F dOmega)
void Evolve_element::compute_stiffness_vector(){

    // initialize the DG_stiffness_vector[i][j] values to zero
    // loop over all the interior nodes of this element
    for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {
        // loop over hidrodynamics indices
        for (int j = 0; j < 4; ++j) {
            this->DG_stiffness_vector[i][j] = 0; 
        }       
    }

    // loop over all the interior nodes of this element
    for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {
        // loop over hidrodynamics indices
        for (int k = 0; k < 4; ++k) {
            // loop over all the interior nodes of this element    
            for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {
                //             S_i              +=                             S^{x}_{ij}                       *               f^{x}_{j}                             +                             S^{y}_{ij}                       *               f^{y}_{j}
                this->DG_stiffness_vector[i][k] += this->this_element->stiffness_matrix_physical_space[0][i][j] * this->this_element->hidrodynamics_vector_f[j][0][k] + this->this_element->stiffness_matrix_physical_space[1][i][j] * this->this_element->hidrodynamics_vector_f[j][1][k];
                // Stiffness vector: DG vector that results from area integral over element of nabla phi_i dot F dOmega.
                // First index runs over interior nodes. 
                // Second index runs between 0 and 3 and represend hidrodynamics variables.
            }       
        }
    }
}

// compute  residual vector: DG vector that results from stiffness vector minus vector result of the numerical flux integration
void Evolve_element::compute_residual_vector(){

    // initialize the DG_residual_vector[i][j] values to zero
    // loop over all the interior nodes of this element
    for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {
        // loop over hidrodynamics indices
        for (int k = 0; k < 4; ++k) {
            this->DG_residual_vector[i][k] = 0; 
        }       
    }

    // compute DG_residual_vector[i][j]
    // loop over all the interior nodes of this element
    for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {
        // loop over hidrodynamics indices
        for (int k = 0; k < 4; ++k) {
            //           R_i               =               S_i               -            int num F_i
            this->DG_residual_vector[i][k] = this->DG_stiffness_vector[i][k] - this->DG_numerical_flux_integration[i][k]; 
        }       
    }
}

// compute residial vector stiffness vector minus vector result of the numerical flux integration
void Evolve_element::compute_time_derivative_U(){

    // initialize the DG_time_derivative_U[i][j] values to zero
    // loop over all the interior nodes of this element
    for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {
        // loop over hidrodynamics indices
        for (int k = 0; k < 4; ++k) {
            this->DG_time_derivative_U[i][k] = 0; 
        }       
    }

    // compute DG_time_derivative_U[i][j]
    // loop over all the interior nodes of this element
    for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {
        // loop over hidrodynamics indices
        for (int k = 0; k < 4; ++k) {
            // loop over all the interior nodes of this element
            for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {
                //         dU_i/d_t              = sum_j                M^{-1}_{ij}                             *            R_{j}  
                this->DG_time_derivative_U[i][k] += this->this_element->inverse_mass_matrix_physical_space[i][j] * this->DG_residual_vector[j][k]; 
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
            this->this_element->hidrodynamics_vector_u[i][k] += time_step * this->DG_time_derivative_U[i][k]; 
        }

        rho = this->this_element->hidrodynamics_vector_u[i][0]; // density
        u   = this->this_element->hidrodynamics_vector_u[i][1] / rho; // horizontal velocity
        v   = this->this_element->hidrodynamics_vector_u[i][2] / rho; // vertical velocity
        E   = this->this_element->hidrodynamics_vector_u[i][3] / rho; // energy
        p   = rho * ( gamma - 1 ) * ( E - ( pow( u , 2) + pow( v , 2) ) / 2 ); // pressure
        H   = E + p / rho; // Entalpy

        // initialize the hidrodinamic vector f, x component  
        this->this_element->hidrodynamics_vector_f[i][0][0] = rho * u;
        this->this_element->hidrodynamics_vector_f[i][0][1] = rho * pow( u , 2 ) + p;
        this->this_element->hidrodynamics_vector_f[i][0][2] = rho * u * v;
        this->this_element->hidrodynamics_vector_f[i][0][3] = rho * u * H;

        // initialize the hidrodinamic vector f, y component  
        this->this_element->hidrodynamics_vector_f[i][1][0] = rho * v;
        this->this_element->hidrodynamics_vector_f[i][1][1] = rho * u * v;
        this->this_element->hidrodynamics_vector_f[i][1][2] = rho * pow( v , 2 ) + p;
        this->this_element->hidrodynamics_vector_f[i][1][3] = rho * v * H;
    }

    // advance in time 
    this->this_element->time += time_step;
}