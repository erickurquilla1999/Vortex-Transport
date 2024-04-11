#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "Element.H"
#include "Evolve.H"

Evolve_element::Evolve_element() {}

Evolve_element::Evolve_element(const Element& this_elem, const Element& right_elem, const Element& left_elem, const Element& vertical_elem):
    
    // Initialize Evolve_element properties    
    this_element(this_elem),
    right_element(right_elem),
    left_element(left_elem),
    vertival_element(vertical_elem)
 
    {
}