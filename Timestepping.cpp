#include <vector>

#include "Element.H"
#include "Evolve.H"

// time stepping with forward euler
void forward_euler(Evolve_element* evo_elemts, double& t_step, const int& n_elements){

    // Compute require quantities for time evolution of Evolve_element objects in the array evolve_elements
    for (int i = 0; i < n_elements ; ++i) {
        evo_elemts[i].compute_U_plus_minus();      // compute U on the element boundaries
        evo_elemts[i].compute_numerical_flux();    // compute numerical flux
        evo_elemts[i].integrate_numerical_flux();  // compute the vector result of integrating the numerical flux
        evo_elemts[i].compute_stiffness_vector();  // compute sitffness vector
        evo_elemts[i].compute_residual_vector();   // compute residual vector
        evo_elemts[i].compute_time_derivative_U(); // compute time derivative of U
    }

    // loop over elements
    for (int i = 0; i < n_elements ; ++i) {
        // Compute new state vectors U and F
        evo_elemts[i].compute_new_U_and_F(t_step);
    }

}

// time stepping with rk4
void rk4(Element* elemts, Evolve_element* evo_elemts, double& t_step, const int& n_elements){



}