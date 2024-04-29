#include <vector>
#include <cmath>

#include "Element.H"
#include "Evolve.H"

// time stepping with forward euler
void forward_euler(Element* elemts, Evolve_element* evo_elemts, double& t_step, const int& n_elements){

    // Compute require quantities for time evolution of Evolve_element objects in the array evolve_elements
    for (int i = 0; i < n_elements ; ++i) {
        evo_elemts[i].compute_U_plus_minus();      // compute U on the element boundaries
        evo_elemts[i].compute_numerical_flux();    // compute numerical flux
        evo_elemts[i].integrate_numerical_flux();  // compute the vector result of integrating the numerical flux
        evo_elemts[i].compute_stiffness_vector();  // compute sitffness vector
        evo_elemts[i].compute_residual_vector();   // compute residual vector
        evo_elemts[i].compute_time_derivative_U(); // compute time derivative of U
    }

    // define hidrodynamics quantities
    double rho, u, v, E, H, p;
    double gamma = 1.4;

    // loop over elements
    for (int n = 0; n < n_elements ; ++n) {
        // loop over all the interior nodes of this element
        for (int i = 0; i < ( elemts[n].p + 1 ) * ( elemts[n].p + 2 ) / 2; ++i) {
            // loop over hidrodynamics indices
            for (int k = 0; k < 4; ++k) {
                elemts[n].hidrodynamics_vector_U[i][k] += t_step * evo_elemts[n].DG_time_derivative_U[i][k]; 
            }

            rho = elemts[n].hidrodynamics_vector_U[i][0]; // density
            u   = elemts[n].hidrodynamics_vector_U[i][1] / rho; // horizontal velocity
            v   = elemts[n].hidrodynamics_vector_U[i][2] / rho; // vertical velocity
            E   = elemts[n].hidrodynamics_vector_U[i][3] / rho; // energy
            p   = rho * ( gamma - 1 ) * ( E - ( pow( u , 2) + pow( v , 2) ) / 2 ); // pressure
            H   = E + p / rho; // Entalpy

            // initialize the hidrodinamic vector f, x component  
            elemts[n].hidrodynamics_vector_F[i][0][0] = rho * u;
            elemts[n].hidrodynamics_vector_F[i][0][1] = rho * pow( u , 2 ) + p;
            elemts[n].hidrodynamics_vector_F[i][0][2] = rho * u * v;
            elemts[n].hidrodynamics_vector_F[i][0][3] = rho * u * H;

            // initialize the hidrodinamic vector f, y component  
            elemts[n].hidrodynamics_vector_F[i][1][0] = rho * v;
            elemts[n].hidrodynamics_vector_F[i][1][1] = rho * u * v;
            elemts[n].hidrodynamics_vector_F[i][1][2] = rho * pow( v , 2 ) + p;
            elemts[n].hidrodynamics_vector_F[i][1][3] = rho * v * H;
        }
        
        // advance in time
        elemts[n].time += t_step;
    }
}

// time stepping with rk4
void rk4(Element* elemts, Evolve_element* evo_elemts, double& t_step, const int& n_elements){



}