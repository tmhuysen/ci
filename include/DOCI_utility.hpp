#ifndef DOCI_DOCI_UTILITY_HPP
#define DOCI_DOCI_UTILITY_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <hf.hpp>
#include <libwint.hpp>

struct CI_basis {
    Eigen::MatrixXd one_int; //one electron integral matrix
    Eigen::Tensor<double, 4> two_int; //two electron integral tensor
    double nuc; //nuclear energy
    size_t n_bf; //number of basis functions
    size_t n_electrons; //number of electrons
};

CI_basis rhf_to_CI_basis(hf::rhf::RHF &rhf_basis);
CI_basis file_to_CI_basis(const std::string& filename);


#endif // DOCI_DOCI_UTILITY_HPP
