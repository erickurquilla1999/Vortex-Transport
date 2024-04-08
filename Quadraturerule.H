#ifndef QUADRATURERULE_H
#define QUADRATURERULE_H

#include <vector>

// return gauss quadrature information for line integral, first index in the gauss quadrature coordinate number, second index runs from 0 to 1. 0 is xi coordinate. 1 is weight.
std::vector<std::vector<double>> gauss_line_integral(const int& IntOrder);

// save gauss quadrature information for area integral, first index in the gauss quadrature coordinate number, second index runs from 0 to 2. 0 is xi coordinate. 1 is eta coordinate. 2 is weight.
std::vector<std::vector<double>> gauss_area_integral(const int& IntOrder);

#endif // QUADRATURERULE_H
