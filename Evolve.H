#ifndef EVOLVE_H
#define EVOLVE_H

#include <vector>

#include "Element.H"

class Evolve_element {
public:

    // Default constructor
    Evolve_element();
 
    // Element default constructor with arguments
    Evolve_element(Element* this_elem, Element* right_elem, Element* left_elem, Element* vertical_elem, const std::vector<std::vector<double>>& gau_int_l);

    // Evolve_element properties
    Element* this_element;     // this is the element to be evolved in time
    Element* right_element;    // this the element to the right to the element to be evolved in time
    Element* left_element;     // this the element to the left to the element to be evolved in time
    Element* vertival_element; // this the element in the vertical direction to the element to be evolved in time
    std::vector<std::vector<double>> gau_integ_line; // contains gauss quadrature information for line integral, first index in the gauss quadrature coordinate number, second index runs from 0 to 1. 0 is xi coordinate. 1 is weight.

};

#endif // EVOLVE_H