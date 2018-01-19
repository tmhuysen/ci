//
// Created by Wulfix on 19/01/2018.
//

#ifndef DOCI_HEAD_DOCI_UTILITY_HPP
#define DOCI_HEAD_DOCI_UTILITY_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

struct CI_basis {
    Eigen::MatrixXd one_int; //one electron integral matrix
    Eigen::Tensor<double, 4> two_int; //two electron integral tensor
    double nuc; //nuclear energy
    unsigned long n_bf; //number of basis functions
    unsigned long n_electrons; //number of electrons
};



#endif //DOCI_HEAD_DOCI_UTILITY_HPP
