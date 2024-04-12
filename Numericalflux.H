#ifndef NUMERICALFLUX_H
#define NUMERICALFLUX_H

// compute numerical flux given the state vector u in the left and right and the normal vector
std::vector<double> numerical_flux(const std::vector<double>& u_left, const std::vector<double>& u_right, const std::vector<double>& normal_vector);

#endif // NUMERICALFLUX_H
