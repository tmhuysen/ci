//
// Created by Wulfix on 19/01/2018.
//

#ifndef DOCI_HEAD_DOCI_UTILITY_HPP
#define DOCI_HEAD_DOCI_UTILITY_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <hf.hpp>
#include <libwrp.hpp>

struct CI_basis {
    Eigen::MatrixXd one_int; //one electron integral matrix
    Eigen::Tensor<double, 4> two_int; //two electron integral tensor
    double nuc; //nuclear energy
    unsigned long n_bf; //number of basis functions
    unsigned long n_electrons; //number of electrons
};

CI_basis rhf_to_CI_basis(hf::rhf::RHF& rhf_basis);

#endif //DOCI_HEAD_DOCI_UTILITY_HPP
