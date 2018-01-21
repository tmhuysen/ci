#ifndef DOCI_EXTRAS_HPP
#define DOCI_EXTRAS_HPP


#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>


struct State {
    double eigenValue;          // The energy of the solution, a.k.a. the eigenvalue
    Eigen::VectorXd eigenVector; // The coefficients of the solution with respect to the given basis, a.k.a.
    // the eigenvector corresponding to the eigenvalue
};



#endif // DOCI_EXTRAS_HPP
