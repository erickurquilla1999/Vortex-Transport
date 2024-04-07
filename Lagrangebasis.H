#ifndef LAGRANGEBASIS_H
#define LAGRANGEBASIS_H

// compute the lagrange basis function in reference space ( xi , eta )
std::vector<double> lagrange_basis_reference_space(const int& p, const std::vector<double>& coords_ref_spa);

// compute the gradiente of the lagrange basis function in reference space ( xi , eta )
std::vector<std::vector<double>> lagrange_basis_gradient_reference_space(const int& p, const std::vector<double>& coords_ref_spa);

// compute coordinates in physical space given ( xi , eta ) in reference space
std::vector<double> reference_to_physical_space(const std::vector<double>& coord_ref_spa,  const std::vector<std::vector<double>>& vertex_phys_spa);

#endif // LAGRANGEBASIS_H
