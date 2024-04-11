#ifndef EVOLVE_H
#define EVOLVE_H

#include <vector>

#include "Element.H"


class Evolve_element {
public:

    // Default constructor
    Evolve_element();
 
    // Element default constructor with arguments
    Evolve_element(const Element& this_elem, const Element& right_elem, const Element& left_elem, const Element& vertical_elem);

    // Evolve_element properties
    Element this_element;
    Element right_element;
    Element left_element;
    Element vertival_element;

};

#endif // EVOLVE_H