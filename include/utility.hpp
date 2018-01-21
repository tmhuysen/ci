#ifndef DOCI_DOCI_UTILITY_HPP
#define DOCI_DOCI_UTILITY_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <hf.hpp>
#include <libwint.hpp>



void symmatu(Eigen::MatrixXd &mat);
void symmatu_reverse(Eigen::MatrixXd &mat);
void symmatu_tensor_reverse(Eigen::Tensor<double, 4> tensor);

bool compareMat(Eigen::MatrixXd &mat1,Eigen::MatrixXd &mat2);
bool compareState(const State &o1, const State &o2);
bool areSame(const State &o1, const State &o2);


#endif // DOCI_DOCI_UTILITY_HPP
