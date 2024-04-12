#ifndef EVOLVE_H
#define EVOLVE_H

#include <vector>

#include "Element.H"

class Evolve_element {
public:

    // Default constructor
    Evolve_element();
 
    // Element default constructor with arguments
    Evolve_element(Element* this_elem, Element* right_elem, Element* left_elem, Element* vertical_elem, const std::vector<std::vector<double>>& gau_int_l,  const int& p_parms);

    // Evolve_element properties
    int p;
    Element* this_element;     // this is the element to be evolved in time
    Element* right_element;    // this the element to the right to the element to be evolved in time
    Element* left_element;     // this the element to the left to the element to be evolved in time
    Element* vertival_element; // this the element in the vertical direction to the element to be evolved in time
    std::vector<std::vector<double>> gau_integ_line; // contains gauss quadrature information for line integral, first index in the gauss quadrature coordinate number, second index runs from 0 to 1. 0 is xi coordinate. 1 is weight.

    // save the lagrange polinomial function evaluated in the quadrature point for the line integral for side 1, 2 and 3
    // this vector store the values of the lagrange polinomial in evaluated in the quadrature points
    // first index runs over interior nodes number, that is, the lagrange polinimial that is one on this node
    // second item runs over the evaluation of the lagrange poliniam in the quadrature points
    std::vector<std::vector<double>> phi_in_quadrature_points_side_1; 
    std::vector<std::vector<double>> phi_in_quadrature_points_side_2;
    std::vector<std::vector<double>> phi_in_quadrature_points_side_3;

    // to following vectors store the interpolation of the hidrodynamic vector u from the interior nodes of the elements to the quadrature points for the line integrals
    // first index runs over the quadrature points number, second index runs over 0 and 3 for the hidrodynamics quantities
    // plus for current element
    std::vector<std::vector<double>> u_plus_side_1; 
    std::vector<std::vector<double>> u_plus_side_2;
    std::vector<std::vector<double>> u_plus_side_3;
    // minus for elements in the boundaries
    std::vector<std::vector<double>> u_minus_side_1;
    std::vector<std::vector<double>> u_minus_side_2;
    std::vector<std::vector<double>> u_minus_side_3;

    // contains the numerical flux (roe flux) at the quadrature points for sides 1, 2 and 3
    // first index runs over the gauss quadrature points number, second index runs over 0 and 3 for the hidrodynamics quantities
    std::vector<std::vector<double>> numerical_flux_side_1;
    std::vector<std::vector<double>> numerical_flux_side_2;
    std::vector<std::vector<double>> numerical_flux_side_3;

    // Evolve_element functions
    void compute_numerical_flux(); // this function compute the numerical flux on the element boundaries, side 1, 2 and 3. 
    // void integrate_numerical_flux(); // initialize the hydronimics quantities u and f
    // void compute_residual_vector(); // initialize the hydronimics quantities u and f
    // void evolve_time_step(); // initialize the hydronimics quantities u and f

};

#endif // EVOLVE_H