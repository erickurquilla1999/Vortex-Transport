#ifndef EVOLVE_H
#define EVOLVE_H

#include <vector>

#include "Element.H"

class Evolve_element {
public:

    // Default constructor
    Evolve_element();
 
    // Element default constructor with arguments
    Evolve_element(Element& this_elem, Element& right_elem, Element& left_elem, Element& vertical_elem);

    // Evolve_element properties
    Element* this_element;
    Element* right_element;
    Element* left_element;
    Element* vertival_element;

};

#endif // EVOLVE_H