#ifndef PREEVOLVE_H
#define PREEVOLVE_H

#include <vector>

// compute the inverse of the mass matrix in reference space 
// mass_ij = int in T phi_i phi_j dT
// T is an triangle in reference space with vertex (0,0), (1,0) and (0,1) in reference space.
// return a matrix of size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2
std::vector<std::vector<double>> inverse_mass_matrix_reference_space(const int& p, const std::vector<std::vector<double>>& gauss_area_int);

// compute stiffness matrix in reference space 
// S_ij = integral in T of ( Nabla phi_i ) phi_j dT
// T is an triangle in reference space with vertex (0,0), (1,0) and (0,1) in reference space.
// return a two dimenional vector in reference space xi and eta, the components of xi and eta are matrices of size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2
std::vector<std::vector<std::vector<double>>> sitffness_matrix_reference_space(const int& p, const std::vector<std::vector<double>>& gauss_area_int);

#endif // PREEVOLVE_H