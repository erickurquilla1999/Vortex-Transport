#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "Element.H"
#include "Lagrangebasis.H"
#include "Utilities.H"

Element::Element() {}

Element::Element(const int& ele_num, const mesh& mesh_info, const std::vector<std::vector<double>>& nods_ref_spa, const parameters& parms):
    
    // Initialize element properties    
    time(0.0), // Simulation initial time
    number(ele_num), // Element number
    type(mesh_info.element_type[ele_num]), // Element type 0 of squere angle is down and 1 if up
    right_element(mesh_info.elements_at_boundary[ele_num][0]), // Elemnt to the right
    left_element(mesh_info.elements_at_boundary[ele_num][1]), // Elemnt to the left
    vertical_element(mesh_info.elements_at_boundary[ele_num][2]), // Elemnt in the vertical direction
    vertices_coords_phys_space(mesh_info.element_coordinates[ele_num]), // Coordinates of the verctices of the element in physical space. First index represent vertices of the element (run between 0 and 2): 0 is the vertice at the square angle, 1 and 2 are the other vertices going counter clock wise. Second index represent x:0 and y:1 position
    nods_coords_refe_space(nods_ref_spa), // Coordinates of the interior nodes of the element in reference space. The first index in the node number, the second inxed runs between 0 and 1. 0: xi and 1: eta.
    nods_coords_phys_space((parms.p + 1) * (parms.p + 2) / 2, std::vector<double>(2)), // Coordinates of the interior nodes of the element in physical space. The first index is the node number, the second inxed runs between 0 and 1. 0: x and 1: y.
    hidrodynamics_vector_u((parms.p + 1) * (parms.p + 2) / 2, std::vector<double>(4)), // Hidrodynamics: vector u for each interior node. The first index in the node number, the second inxed runs between 0 and 4 for hidrodynamics quantities. 0: rho, 1: rho u, 2: rho v, 3: rho E.
    hidrodynamics_vector_f((parms.p + 1) * (parms.p + 2) / 2, std::vector<std::vector<double>>(2, std::vector<double>(4))), // Hidrodynamics: vector f for each interior node. The first index in the node number. The second index runs between 0 and 1 and represente physical space 0:x and 1: y. The third inxed runs between 0 and 4 for hidrogynamics quantities.
    jacobian(2,std::vector<double>(2)), // jacobian between transformation from reference space to physical space d vec{x} / d vec{xi} = [ [ x2 - x1 , x3 - x1 ] , [ y2 - y1 , y3 - y1 ] ]    
    inverse_jacobian(2,std::vector<double>(2)), // jacobian between transformation from physical space to reference space d vec{xi} / d vec{x} = ( 1 / det( J ) ) * [ [ y3 - y1 , x1 - x3 ] , [ y1 - y2 , x2 - x1 ] ]    
    inverse_mass_matrix_physical_space((parms.p + 1) * (parms.p + 2) / 2 , std::vector<double>((parms.p + 1) * (parms.p + 2) / 2)), // mass_ij = int in T phi_i phi_j dT . size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2
    stiffness_matrix_physical_space(2, std::vector<std::vector<double>>( (parms.p + 1) * (parms.p + 2) / 2 , std::vector<double>( (parms.p + 1) * (parms.p + 2) / 2 ))) // S_ij = integral in T of ( Nabla phi_i ) phi_j dT . form :  hat{e}_x * matrix[ ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ] + hat{e}_y * matrix[ ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ] 

    {
    
    // compute interior nodes coordinate in physical space
    for (int i = 0; i < ( parms.p + 1 ) *( parms.p + 2 ) / 2 ; ++i) {
        this->nods_coords_phys_space[i] = reference_to_physical_space(nods_coords_refe_space[i], vertices_coords_phys_space);
    }

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

}

// initialize the hydronimics quantities u and f
void Element::initialize_hydrodinamics(const parameters& parms){

    for (int i = 0; i < ( parms.p + 1 ) *( parms.p + 2 ) / 2 ; ++i) {

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
        p_infty = 20 / 7;
        U_infty = 1 / pow( 2 , 0.5);
        V_infty = 1 / pow( 2 , 0.5);
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
void Element::write_data(const parameters& parms, const int& step_num){

    // prepare data to be saved
    std::vector<std::string> lines( 1 + ( parms.p + 1 ) *( parms.p + 2 ) / 2 );
    lines[0]="node_number time x y u0 u1 u2 u3 fx0 fx1 fx2 fx3 fy0 fy1 fy2 fy3";

    for (int i = 0; i < ( parms.p + 1 ) *( parms.p + 2 ) / 2 ; ++i) {

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

    writeToFile("output/step_" + std::to_string( step_num ) + "/element_" + std::to_string( this->number) + ".txt", lines);

    }

}