#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "Element.H"
#include "Evolve.H"

Evolve_element::Evolve_element() {}

Evolve_element::Evolve_element(Element& this_elem, Element& right_elem, Element& left_elem, Element& vertical_elem):
    
    // Initialize Evolve_element properties    
    this_element(&this_elem),
    right_element(&right_elem),
    left_element(&left_elem),
    vertival_element(&vertical_elem)
 
    {

    std::cout << "Element : " << this->this_element->hidrodynamics_vector_u[0][0] << std::endl;
    this->this_element->hidrodynamics_vector_u[0][0] = 9999999.9;

}