#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "Element.H"
#include "Evolve.H"

Evolve_element::Evolve_element() {}

Evolve_element::Evolve_element(Element* this_elem, Element* right_elem, Element* left_elem, Element* vertical_elem, const std::vector<std::vector<double>>& gau_int_l):
    
    // Initialize Evolve_element properties    
    this_element(this_elem),            // this is the element to be evolved in time
    right_element(right_elem),          // this the element to the right to the element to be evolved in time
    left_element(left_elem),            // this the element to the left to the element to be evolved in time
    vertival_element(vertical_elem),    // this the element in the vertical direction to the element to be evolved in time
    gau_integ_line(gau_int_l)           // 

    {

    std::cout << "Element : " << this->this_element->hidrodynamics_vector_u[0][0] << std::endl;
    this->this_element->hidrodynamics_vector_u[0][0] = 9999999.9;

}