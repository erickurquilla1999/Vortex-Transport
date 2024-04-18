#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "Element.H"
#include "Lagrangebasis.H"
#include "Utilities.H"

Element::Element() {}

Element::Element(const int& ele_num, const mesh& mesh_info, const std::vector<std::vector<double>>& nods_ref_spa, const int& p_lagrange):
    
    // Initialize element properties    
    time(0.0), // Simulation initial time
    p(p_lagrange), // lagrange polinomial order
    number(ele_num), // Element number
    type(mesh_info.element_type[ele_num]), // Element type 0 of squere angle is down and 1 if up
    right_element(mesh_info.elements_at_boundary[ele_num][0]), // Elemnt to the right
    left_element(mesh_info.elements_at_boundary[ele_num][1]), // Elemnt to the left
    vertical_element(mesh_info.elements_at_boundary[ele_num][2]), // Elemnt in the vertical direction
    vertices_coords_phys_space(mesh_info.element_coordinates[ele_num]), // Coordinates of the verctices of the element in physical space. First index represent vertices of the element (run between 0 and 2): 0 is the vertice at the square angle, 1 and 2 are the other vertices going counter clock wise. Second index represent x:0 and y:1 position
    nods_coords_refe_space(nods_ref_spa), // Coordinates of the interior nodes of the element in reference space. The first index in the node number, the second inxed runs between 0 and 1. 0: xi and 1: eta.
    nods_coords_phys_space((this->p + 1) * (this->p + 2) / 2, std::vector<double>(2)), // Coordinates of the interior nodes of the element in physical space. The first index is the node number, the second inxed runs between 0 and 1. 0: x and 1: y.
    hidrodynamics_vector_u((this->p + 1) * (this->p + 2) / 2, std::vector<double>(4)), // Hidrodynamics: vector u for each interior node. The first index in the node number, the second inxed runs between 0 and 4 for hidrodynamics quantities. 0: rho, 1: rho u, 2: rho v, 3: rho E.
    hidrodynamics_vector_f((this->p + 1) * (this->p + 2) / 2, std::vector<std::vector<double>>(2, std::vector<double>(4))), // Hidrodynamics: vector f for each interior node. The first index in the node number. The second index runs between 0 and 1 and represente physical space 0:x and 1: y. The third inxed runs between 0 and 4 for hidrogynamics quantities.
    jacobian(2,std::vector<double>(2)), // jacobian between transformation from reference space to physical space d vec{x} / d vec{xi} = [ [ x2 - x1 , x3 - x1 ] , [ y2 - y1 , y3 - y1 ] ]
    inverse_jacobian(2,std::vector<double>(2)), // jacobian between transformation from physical space to reference space d vec{xi} / d vec{x} = ( 1 / det( J ) ) * [ [ y3 - y1 , x1 - x3 ] , [ y1 - y2 , x2 - x1 ] ]
    inverse_mass_matrix_physical_space((this->p + 1) * (this->p + 2) / 2 , std::vector<double>((this->p + 1) * (this->p + 2) / 2)), // mass_ij = int in Omega phi_i phi_j dOmega . Size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2
    stiffness_matrix_physical_space(2, std::vector<std::vector<double>>( (this->p + 1) * (this->p + 2) / 2 , std::vector<double>( (this->p + 1) * (this->p + 2) / 2 ))), // S_ij = integral in Omega of ( Nabla phi_i ) phi_j dOmega . form :  hat{e}_x * matrix[ ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ] + hat{e}_y * matrix[ ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ]. first index run between spacial components in physics space. 0: x and 1 y. second and third index run over matrix inidices of size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ] 
    units_vectors_perpendicular_to_element_boundary(3,std::vector<double>(2)), // contains the units vectors perperdicular to the elements boundary. the first index runs from 0 two 2 and represent the side number. side 1 is the one found going counterclockwise from the initial vextex (square angle vertex), in continuation side 2 and 3 going conter clockwise. the second index runs from 0 to 1, 0 is the x component of the unit vector 1 is the y components.
    sides_lenght(3) // contains the element side lenghts, this array contains just three values. 0: side 1, 1: side 2, 2: side 3

    {
    
    // compute interior nodes coordinate in physical space
    for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
        this->nods_coords_phys_space[i] = reference_to_physical_space(nods_coords_refe_space[i], vertices_coords_phys_space);
    }
}

// compute jacobians to connetc referece space to physical space and viceversa for each element
void Element::build_jacobians(){

    // compute jacobian between transformation from reference space to physical space d vec{x} / d vec{xi} = [ [ x2 - x1 , x3 - x1 ] , [ y2 - y1 , y3 - y1 ] ] 
    this->jacobian[0][0] = this->vertices_coords_phys_space[1][0] - this->vertices_coords_phys_space[0][0];
    this->jacobian[0][1] = this->vertices_coords_phys_space[2][0] - this->vertices_coords_phys_space[0][0];
    this->jacobian[1][0] = this->vertices_coords_phys_space[1][1] - this->vertices_coords_phys_space[0][1];
    this->jacobian[1][1] = this->vertices_coords_phys_space[2][1] - this->vertices_coords_phys_space[0][1];

    // compute determinant of jacobian
    this->determinant_jacobian = this->jacobian[0][0] * this->jacobian[1][1] - this->jacobian[1][0] * this->jacobian[0][1];

    // compute jacobian between transformation from physical space to reference space d vec{xi} / d vec{x} = ( 1 / det( J ) ) * [ [ y3 - y1 , x1 - x3 ] , [ y1 - y2 , x2 - x1 ] ]
    this->inverse_jacobian[0][0] = ( 1 / this->determinant_jacobian ) * ( this->vertices_coords_phys_space[2][1] - this->vertices_coords_phys_space[0][1] );
    this->inverse_jacobian[0][1] = ( 1 / this->determinant_jacobian ) * ( this->vertices_coords_phys_space[0][0] - this->vertices_coords_phys_space[2][0] );
    this->inverse_jacobian[1][0] = ( 1 / this->determinant_jacobian ) * ( this->vertices_coords_phys_space[0][1] - this->vertices_coords_phys_space[1][1] );
    this->inverse_jacobian[1][1] = ( 1 / this->determinant_jacobian ) * ( this->vertices_coords_phys_space[1][0] - this->vertices_coords_phys_space[0][0] );

    // compute determinant of inverse jacobian
    this->determinant_inverse_jacobian = this->inverse_jacobian[0][0] * this->inverse_jacobian[1][1] - this->inverse_jacobian[1][0] * this->inverse_jacobian[0][1];

    // compute unit vectors perpendicular to the element edges
    //                                                      (            y2                     -              y1                 ) / (       (             x2                   -           x1                     )^2   +    (             y2                   -           y1                     )^2   )^1/2             
    this->units_vectors_perpendicular_to_element_boundary[0][0] = ( this->vertices_coords_phys_space[1][1] - this->vertices_coords_phys_space[0][1] ) / pow( pow( this->vertices_coords_phys_space[1][0] - this->vertices_coords_phys_space[0][0] , 2 ) + pow( this->vertices_coords_phys_space[1][1] - this->vertices_coords_phys_space[0][1] , 2 ) , 0.5 );
    this->units_vectors_perpendicular_to_element_boundary[1][0] = ( this->vertices_coords_phys_space[2][1] - this->vertices_coords_phys_space[1][1] ) / pow( pow( this->vertices_coords_phys_space[2][0] - this->vertices_coords_phys_space[1][0] , 2 ) + pow( this->vertices_coords_phys_space[2][1] - this->vertices_coords_phys_space[1][1] , 2 ) , 0.5 );
    this->units_vectors_perpendicular_to_element_boundary[2][0] = ( this->vertices_coords_phys_space[0][1] - this->vertices_coords_phys_space[2][1] ) / pow( pow( this->vertices_coords_phys_space[0][0] - this->vertices_coords_phys_space[2][0] , 2 ) + pow( this->vertices_coords_phys_space[0][1] - this->vertices_coords_phys_space[2][1] , 2 ) , 0.5 );

    //                                                      -1 * (            x2                     -              x1                 ) / (       (             x2                   -           x1                     )^2   +    (             y2                   -           y1                     )^2   )^1/2             
    this->units_vectors_perpendicular_to_element_boundary[0][1] = -1 * ( this->vertices_coords_phys_space[1][0] - this->vertices_coords_phys_space[0][0] ) / pow( pow( this->vertices_coords_phys_space[1][0] - this->vertices_coords_phys_space[0][0] , 2 ) + pow( this->vertices_coords_phys_space[1][1] - this->vertices_coords_phys_space[0][1] , 2 ) , 0.5 );
    this->units_vectors_perpendicular_to_element_boundary[1][1] = -1 * ( this->vertices_coords_phys_space[2][0] - this->vertices_coords_phys_space[1][0] ) / pow( pow( this->vertices_coords_phys_space[2][0] - this->vertices_coords_phys_space[1][0] , 2 ) + pow( this->vertices_coords_phys_space[2][1] - this->vertices_coords_phys_space[1][1] , 2 ) , 0.5 );
    this->units_vectors_perpendicular_to_element_boundary[2][1] = -1 * ( this->vertices_coords_phys_space[0][0] - this->vertices_coords_phys_space[2][0] ) / pow( pow( this->vertices_coords_phys_space[0][0] - this->vertices_coords_phys_space[2][0] , 2 ) + pow( this->vertices_coords_phys_space[0][1] - this->vertices_coords_phys_space[2][1] , 2 ) , 0.5 );

    // compute the element side lenghts
    //                   (     (           x2                    -             x1                   )^2   +    (             y2                   -              y1                  )^2   )^0.5      
    this->sides_lenght[0] = pow( pow( this->vertices_coords_phys_space[1][0] - this->vertices_coords_phys_space[0][0] , 2 ) + pow( this->vertices_coords_phys_space[1][1] - this->vertices_coords_phys_space[0][1] , 2 ) , 0.5 ); // side 1
    this->sides_lenght[1] = pow( pow( this->vertices_coords_phys_space[2][0] - this->vertices_coords_phys_space[1][0] , 2 ) + pow( this->vertices_coords_phys_space[2][1] - this->vertices_coords_phys_space[1][1] , 2 ) , 0.5 ); // side 2
    this->sides_lenght[2] = pow( pow( this->vertices_coords_phys_space[0][0] - this->vertices_coords_phys_space[2][0] , 2 ) + pow( this->vertices_coords_phys_space[0][1] - this->vertices_coords_phys_space[2][1] , 2 ) , 0.5 ); // side 3

}

// builds mass matrix inverse from referece space to physical space for each element
void Element::build_mass_matrix_inverse(const std::vector<std::vector<double>>& inv_mass_matrix){
    // compute inverse mass matrix in physical space
    for (int i = 0; i < (this->p + 1) * (this->p + 2) / 2; ++i) {
        for (int j = 0; j < (this->p + 1) * (this->p + 2) / 2; ++j) {
            this->inverse_mass_matrix_physical_space[i][j] = this->determinant_jacobian * inv_mass_matrix[i][j];     
        }
    }
}

// builds stiffness matrix from referece space to physical space for each element
void Element::build_stiffness_matrix(const std::vector<std::vector<std::vector<double>>>& stiff_matrix){
    // compute stiffness matrix in physical space 
    for (int i = 0; i < (this->p + 1) * (this->p + 2) / 2; ++i) {
        for (int j = 0; j < (this->p + 1) * (this->p + 2) / 2; ++j) {
            this->stiffness_matrix_physical_space[0][i][j] = this->determinant_jacobian * ( stiff_matrix[0][i][j] * this->inverse_jacobian[0][0] + stiff_matrix[1][i][j] * this->inverse_jacobian[1][0] );     
            this->stiffness_matrix_physical_space[1][i][j] = this->determinant_jacobian * ( stiff_matrix[0][i][j] * this->inverse_jacobian[0][1] + stiff_matrix[1][i][j] * this->inverse_jacobian[1][1] );
        }
    }
}

// initialize the hydronimics quantities u and f
void Element::initialize_hydrodinamics(){

    for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {

        double x, y, t;
        x = this->nods_coords_phys_space[i][0]; // x position of node i
        y = this->nods_coords_phys_space[i][1]; // y position of node i
        t = this->time; // initial time

        double rho_infty, rc, epsilon, gamma, M_infty, p_infty, U_infty, V_infty, x0, y0;

        rho_infty = 1.0;
        rc = 1.0;
        epsilon = 0.3;
        gamma = 1.4;
        M_infty = 0.5;
        p_infty = 20.0 / 7.0;
        U_infty = 1.0 / pow( 2.0 , 0.5);
        V_infty = 1.0 / pow( 2.0 , 0.5);
        x0 = 0.0;
        y0 = 0.0;

        double f0, f1, f2;

        f0 = 1 - ( pow ( x - x0 - U_infty * t , 2 ) + pow ( y - y0 - V_infty * t , 2 ) ) / pow ( rc , 2 );
        f1 = 1 - pow ( epsilon , 2 ) * ( gamma -1 ) * pow ( M_infty , 2) * exp( f0 ) / ( 8 * pow ( M_PI , 2 ) );
        f2 = epsilon * ( pow( U_infty , 2) + pow( V_infty , 2) ) * exp( f0 / 2 ) / ( 2 * M_PI * rc );
 
        // hidrodynamic quantities 
        double rho, u, v, p, E, H;

        rho = rho_infty * pow( f1 , 1 / ( gamma - 1 ) ); // density
        u   = U_infty - f2 * ( y - y0 - V_infty * t ); // horizontal velocity
        v   = V_infty - f2 * ( x - x0 - U_infty * t ); // vertical velocity
        p   = p_infty * pow( f1 , gamma / ( gamma - 1 ) ); // pressure
        E   = p / ( rho * ( gamma -1 ) ) + ( pow( u , 2) + pow( v , 2) ) / 2; // Energy
        H   = E + p / rho; // Entalpy

        // initialize the hidrodinamic vector u  
        this->hidrodynamics_vector_u[i][0] = rho;
        this->hidrodynamics_vector_u[i][1] = rho * u;
        this->hidrodynamics_vector_u[i][2] = rho * v;
        this->hidrodynamics_vector_u[i][3] = rho * E;

        // initialize the hidrodinamic vector f, x component  
        this->hidrodynamics_vector_f[i][0][0] = rho * u;
        this->hidrodynamics_vector_f[i][0][1] = rho * pow( u , 2 ) + p;
        this->hidrodynamics_vector_f[i][0][2] = rho * u * v;
        this->hidrodynamics_vector_f[i][0][3] = rho * u * H;

        // initialize the hidrodinamic vector f, y component  
        this->hidrodynamics_vector_f[i][1][0] = rho * v;
        this->hidrodynamics_vector_f[i][1][1] = rho * u * v;
        this->hidrodynamics_vector_f[i][1][2] = rho * pow( v , 2 ) + p;
        this->hidrodynamics_vector_f[i][1][3] = rho * v * H;

    }

}

// write element data in output directory
void Element::write_data(const int& step_num){

    if ( this->number == 0 ){
        std::cout << "\nnumber : " << this->number << std::endl;
        std::cout << "time : " << this->time << std::endl;
        std::cout << "type : " << this->type << std::endl;
        std::cout << "right_element : " << this->right_element << std::endl;
        std::cout << "left_element : " << this->left_element << std::endl;
        std::cout << "vertical_element : " << this->vertical_element << std::endl;
        std::cout << "vertices_coords_phys_space[0][0] : " << this->vertices_coords_phys_space[0][0] << std::endl;
        std::cout << "vertices_coords_phys_space[0][1] : " << this->vertices_coords_phys_space[0][1] << std::endl;
        std::cout << "vertices_coords_phys_space[1][0] : " << this->vertices_coords_phys_space[1][0] << std::endl;
        std::cout << "vertices_coords_phys_space[1][1] : " << this->vertices_coords_phys_space[1][1] << std::endl;
        std::cout << "vertices_coords_phys_space[2][0] : " << this->vertices_coords_phys_space[2][0] << std::endl;
        std::cout << "vertices_coords_phys_space[2][1] : " << this->vertices_coords_phys_space[2][1] << std::endl;
        for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
        std::cout << "nods_coords_refe_space["<<i<<"][0]: " << this->nods_coords_refe_space[i][0] << std::endl;
        std::cout << "nods_coords_refe_space["<<i<<"][1]: " << this->nods_coords_refe_space[i][1] << std::endl;
        }
        for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
        std::cout << "nods_coords_phys_space["<<i<<"][0]: " << this->nods_coords_phys_space[i][0] << std::endl;
        std::cout << "nods_coords_phys_space["<<i<<"][1]: " << this->nods_coords_phys_space[i][1] << std::endl;
        }
        std::cout << "units_vectors_perpendicular_to_element_boundary[0][0] : " << this->units_vectors_perpendicular_to_element_boundary[0][0] << std::endl;
        std::cout << "units_vectors_perpendicular_to_element_boundary[0][1] : " << this->units_vectors_perpendicular_to_element_boundary[0][1] << std::endl;
        std::cout << "units_vectors_perpendicular_to_element_boundary[1][0] : " << this->units_vectors_perpendicular_to_element_boundary[1][0] << std::endl;
        std::cout << "units_vectors_perpendicular_to_element_boundary[1][1] : " << this->units_vectors_perpendicular_to_element_boundary[1][1] << std::endl;
        std::cout << "units_vectors_perpendicular_to_element_boundary[2][0] : " << this->units_vectors_perpendicular_to_element_boundary[2][0] << std::endl;
        std::cout << "units_vectors_perpendicular_to_element_boundary[2][1] : " << this->units_vectors_perpendicular_to_element_boundary[2][1] << std::endl;
        std::cout << "sides_lenght[0] : " << this->sides_lenght[0] << std::endl;
        std::cout << "sides_lenght[1] : " << this->sides_lenght[1] << std::endl;
        std::cout << "sides_lenght[2] : " << this->sides_lenght[2] << std::endl;
        std::cout << "jacobian[0][0] : " << this->jacobian[0][0] << std::endl;
        std::cout << "jacobian[0][1] : " << this->jacobian[0][1] << std::endl;
        std::cout << "jacobian[1][0] : " << this->jacobian[1][0] << std::endl;
        std::cout << "jacobian[1][1] : " << this->jacobian[1][1] << std::endl;
        std::cout << "determinant_jacobian : " << this->determinant_jacobian << std::endl;
        std::cout << "inverse_jacobian[0][0] : " << this->inverse_jacobian[0][0] << std::endl;
        std::cout << "inverse_jacobian[0][1] : " << this->inverse_jacobian[0][1] << std::endl;
        std::cout << "inverse_jacobian[1][0] : " << this->inverse_jacobian[1][0] << std::endl;
        std::cout << "inverse_jacobian[1][1] : " << this->inverse_jacobian[1][1] << std::endl;
        std::cout << "determinant_inverse_jacobian : " << this->determinant_inverse_jacobian << std::endl;
        for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
            for (int j = 0; j < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++j) {
                std::cout << "inverse_mass_matrix_physical_space["<<i<<"]["<<j<<"] : " << this->inverse_mass_matrix_physical_space[i][j] << std::endl;
            }
        }
        for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
            for (int j = 0; j < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++j) {
                std::cout << "stiffness_matrix_physical_space[0]["<<i<<"]["<<j<<"] : " << this->stiffness_matrix_physical_space[0][i][j] << std::endl;
            }
        }
        for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
            for (int j = 0; j < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++j) {
                std::cout << "stiffness_matrix_physical_space[1]["<<i<<"]["<<j<<"] : " << this->stiffness_matrix_physical_space[1][i][j] << std::endl;
            }
        }
        for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
        std::cout << "hidrodynamics_vector_u["<<i<<"][0]: " << this->hidrodynamics_vector_u[i][0] << std::endl;
        std::cout << "hidrodynamics_vector_u["<<i<<"][1]: " << this->hidrodynamics_vector_u[i][1] << std::endl;
        std::cout << "hidrodynamics_vector_u["<<i<<"][2]: " << this->hidrodynamics_vector_u[i][2] << std::endl;
        std::cout << "hidrodynamics_vector_u["<<i<<"][3]: " << this->hidrodynamics_vector_u[i][3] << std::endl;
        }
        for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
        std::cout << "hidrodynamics_vector_f["<<i<<"][0][0]: " << this->hidrodynamics_vector_f[i][0][0] << std::endl;
        std::cout << "hidrodynamics_vector_f["<<i<<"][0][1]: " << this->hidrodynamics_vector_f[i][0][1] << std::endl;
        std::cout << "hidrodynamics_vector_f["<<i<<"][0][2]: " << this->hidrodynamics_vector_f[i][0][2] << std::endl;
        std::cout << "hidrodynamics_vector_f["<<i<<"][0][3]: " << this->hidrodynamics_vector_f[i][0][3] << std::endl;
        std::cout << "hidrodynamics_vector_f["<<i<<"][1][0]: " << this->hidrodynamics_vector_f[i][1][0] << std::endl;
        std::cout << "hidrodynamics_vector_f["<<i<<"][1][1]: " << this->hidrodynamics_vector_f[i][1][1] << std::endl;
        std::cout << "hidrodynamics_vector_f["<<i<<"][1][2]: " << this->hidrodynamics_vector_f[i][1][2] << std::endl;
        std::cout << "hidrodynamics_vector_f["<<i<<"][1][3]: " << this->hidrodynamics_vector_f[i][1][3] << std::endl;
        }
    }
    
    // prepare data to be saved
    std::vector<std::string> lines( 1 + ( this->p + 1 ) *( this->p + 2 ) / 2 );
    lines[0]="node_number time x y u0 u1 u2 u3 fx0 fx1 fx2 fx3 fy0 fy1 fy2 fy3";

    for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {

        lines[i+1]=std::to_string(i)+" "; // node number
        lines[i+1]+=std::to_string(this->time)+" "; // time
        lines[i+1]+=std::to_string(this->nods_coords_phys_space[i][0])+" "; // x
        lines[i+1]+=std::to_string(this->nods_coords_phys_space[i][1])+" "; // y 
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_u[i][0])+" "; // u0
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_u[i][1])+" "; // u1
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_u[i][2])+" "; // u2
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_u[i][3])+" "; // u3
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_f[i][0][0])+" "; // fx0
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_f[i][0][1])+" "; // fx1
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_f[i][0][2])+" "; // fx2
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_f[i][0][3])+" "; // fx3
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_f[i][1][0])+" "; // fy0
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_f[i][1][1])+" "; // fy1
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_f[i][1][2])+" "; // fy2
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_f[i][1][3])+" "; // fy3

    }

    writeToFile("output/step_" + std::to_string( step_num ) + "/element_" + std::to_string( this->number) + ".txt", lines);

    // if step is zero it writes the jacobians and their determinants, inverse mass and stiffness matrix
    if ( step_num == 0 ){

        // prepare jacobians and their determinants, inverse mass and stiffness matrix to be saved

        std::vector<std::string> lines_( 13 + 3 * ( this->p + 1 ) *( this->p + 2 ) / 2 + 3 );
        
        lines_[0]="jacobian";
        lines_[1]="[ [ " + std::to_string( this->jacobian[0][0] ) + " , " + std::to_string( this->jacobian[0][1] ) + " ] , ";
        lines_[2]="  [ " + std::to_string( this->jacobian[1][0] ) + " , " + std::to_string( this->jacobian[1][1] ) + " ] ]";
        lines_[3]="determinant of jacobian";
        lines_[4]=std::to_string( this->determinant_jacobian );
        lines_[5]="inverse jacobian";
        lines_[6]="[ [ " + std::to_string( this->inverse_jacobian[0][0] ) + " , " + std::to_string( this->inverse_jacobian[0][1] ) + " ] , ";
        lines_[7]="  [ " + std::to_string( this->inverse_jacobian[1][0] ) + " , " + std::to_string( this->inverse_jacobian[1][1] ) + " ] ]";
        lines_[8]="determinant of inverse jacobian";
        lines_[9]=std::to_string( this->determinant_inverse_jacobian );
        lines_[10]="inverse mass matrix";
        lines_[11]="[ ";
        for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
            lines_[11+i]+="[ " + std::to_string(this->inverse_mass_matrix_physical_space[i][0]);
            for (int j = 1; j < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++j) {
                lines_[11+i]+=", " + std::to_string(this->inverse_mass_matrix_physical_space[i][j]);
            }
            lines_[11+i]+=" ] ,";
        }
        lines_[11 + ( this->p + 1 ) *( this->p + 2 ) / 2 -1 ].erase(lines_[11 + ( this->p + 1 ) *( this->p + 2 ) / 2 -1 ].size() - 1);
        lines_[11 + ( this->p + 1 ) *( this->p + 2 ) / 2 -1 ] += "]";
        lines_[11 + ( this->p + 1 ) *( this->p + 2 ) / 2]="stiffness matrix : x component";
        lines_[12 + ( this->p + 1 ) *( this->p + 2 ) / 2]="[ ";
        for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
            lines_[12 + ( this->p + 1 ) *( this->p + 2 ) / 2 + i]+="[ " + std::to_string(this->stiffness_matrix_physical_space[0][i][0]);
            for (int j = 1; j < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++j) {
                lines_[12 + ( this->p + 1 ) *( this->p + 2 ) / 2 + i]+=", " + std::to_string(this->stiffness_matrix_physical_space[0][i][j]);
            }
            lines_[12 + ( this->p + 1 ) *( this->p + 2 ) / 2 + i]+=" ] ,";
        }
        lines_[12 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ].erase(lines_[ 12 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ].size() - 1);
        lines_[12 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ] += "]";
        lines_[12 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2]="stiffness matrix : y component";
        lines_[13 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2]="[ ";
        for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
            lines_[13 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2 + i]+="[ " + std::to_string(this->stiffness_matrix_physical_space[1][i][0]);
            for (int j = 1; j < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++j) {
                lines_[13 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2 + i]+=", " + std::to_string(this->stiffness_matrix_physical_space[1][i][j]);
            }
            lines_[13 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2 + i]+=" ] ,";
        }
        lines_[13 + 3 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ].erase(lines_[ 13 + 3 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ].size() - 1);
        lines_[13 + 3 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ] += "]";

        lines_[1 + 13 + 3 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ] = "Unit vector perpendicular to edge 1 : [ " + std::to_string( this->units_vectors_perpendicular_to_element_boundary[0][0] ) + " , " + std::to_string( this->units_vectors_perpendicular_to_element_boundary[0][1] ) + " ]";
        lines_[2 + 13 + 3 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ] = "Unit vector perpendicular to edge 2 : [ " + std::to_string( this->units_vectors_perpendicular_to_element_boundary[1][0] ) + " , " + std::to_string( this->units_vectors_perpendicular_to_element_boundary[1][1] ) + " ]";
        lines_[3 + 13 + 3 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ] = "Unit vector perpendicular to edge 3 : [ " + std::to_string( this->units_vectors_perpendicular_to_element_boundary[2][0] ) + " , " + std::to_string( this->units_vectors_perpendicular_to_element_boundary[2][1] ) + " ]";

        writeToFile("output/step_" + std::to_string( step_num ) + "/JMS_element_" + std::to_string( this->number) + ".txt", lines_);

    }

}