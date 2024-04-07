#ifndef LAGRANGEBASIS_H
#define LAGRANGEBASIS_H

// compute the lagrange basis function in reference space ( xi , eta )
std::vector<double> lagrange_basis_reference_space(int& p, const std::vector<double>& coords_ref_spa);

// compute the gradiente of the lagrange basis function in reference space ( xi , eta )
std::vector<std::vector<double>> lagrange_basis_gradient_reference_space(int& p, const std::vector<double>& coords_ref_spa);

#endif // LAGRANGEBASIS_H
