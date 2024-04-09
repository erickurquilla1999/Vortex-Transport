#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "Lagrangebasis.H"

// compute the inverse of the mass matrix in reference space 
// mass_ij = int in T phi_i phi_j dT
// T is an triangle in reference space with vertex (0,0), (1,0) and (0,1) in reference space.
// return a matrix of size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2
std::vector<std::vector<double>> inverse_mass_matrix_reference_space(const int& p, const std::vector<std::vector<double>>& gauss_area_int){

    // get size, that is the number of quadrature point for integration
    int size = gauss_area_int.size();

    // this vector store the values of the lagrange polinomial in evaluated in the quadrature points
    // first index runs over interior nodes number, that is, the lagrange polinimial that is one on this node
    // second item runs over the evaluation of the lagrange poliniam in the quadrature points
    std::vector<std::vector<double>> phi_in_quadrature_points( ( p + 1 ) * ( p + 2 ) / 2 , std::vector<double>( size ) );

    // initialize counter
    int counter = 0;

    // evaluate the lagrange polinomial in the quadrature points
    for (int i = 0; i < size; ++i) {

        std::vector<double> xi_eta_gauss(2);
        xi_eta_gauss[0] = gauss_area_int[i][0];
        xi_eta_gauss[1] = gauss_area_int[i][1];

        std::vector<double> phi_in_xi_eta_gauss = lagrange_basis_reference_space( p , xi_eta_gauss ); 

        for (int j = 0; j < ( p + 1 ) * ( p + 2 ) / 2; ++j) {
            phi_in_quadrature_points[j][counter] = phi_in_xi_eta_gauss[j];
        }
        counter++;
    }

    // define mass matrix
    std::vector<std::vector<double>> mass_matrix( ( p + 1 ) * ( p + 2 ) / 2 , std::vector<double>( ( p + 1 ) * ( p + 2 ) / 2 ) );

    // integrate mass matrix using gauss methods
    // mass_ij = int in T phi_i phi_j dT , T is an triangle with vertex (0,0), (1,0) and (0,1) in reference space.
    for (int i = 0; i < ( p + 1 ) * ( p + 2 ) / 2; ++i) {
        for (int j = 0; j < ( p + 1 ) * ( p + 2 ) / 2; ++j) {
            mass_matrix[i][j] = 0;
            for (int n = 0; n < size; ++n) {
                mass_matrix[i][j] += phi_in_quadrature_points[i][n] * phi_in_quadrature_points[j][n] * gauss_area_int[n][2];
            }
        }
    } 

    // Convert the mass matrix to an Eigen matrix just to compute the inverse matrix
    Eigen::MatrixXd eigen_mass_matrix( ( p + 1 ) * ( p + 2 ) / 2 ,  ( p + 1 ) * ( p + 2 ) / 2 );
    for (int i = 0; i <  ( p + 1 ) * ( p + 2 ) / 2 ; ++i) {
        for (int j = 0; j <  ( p + 1 ) * ( p + 2 ) / 2 ; ++j) {
            eigen_mass_matrix(i, j) = mass_matrix[i][j];
        }
    }

    // Compute the inverse of the mass matrix using Eigen
    Eigen::MatrixXd inv_matrix = eigen_mass_matrix.inverse();

    // Output mass matrix and its inverse
    std::cout << "Mass matrix:" << std::endl << eigen_mass_matrix << std::endl;
    std::cout << "Inverse mass matrix:" << std::endl << inv_matrix << std::endl;

    // Convert the inv_matrix to an vector matrix to return the right function type
    std::vector<std::vector<double>> mass_matrix_inverse( ( p + 1 ) * ( p + 2 ) / 2 , std::vector<double>( ( p + 1 ) * ( p + 2 ) / 2 ) );
    for (int i = 0; i <  ( p + 1 ) * ( p + 2 ) / 2 ; ++i) {
        for (int j = 0; j <  ( p + 1 ) * ( p + 2 ) / 2 ; ++j) {
            mass_matrix_inverse[i][j] = inv_matrix(i, j);
        }
    }

    // return inverse mass matrix
    return mass_matrix_inverse;

}