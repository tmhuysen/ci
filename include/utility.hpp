#ifndef DOCI_DOCI_UTILITY_HPP
#define DOCI_DOCI_UTILITY_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>


namespace doci {

void symmatu(Eigen::MatrixXd& mat);  // FIXME: this function is never used
void symmatu_reverse(Eigen::MatrixXd& mat);

void symmatu_tensor_reverse(Eigen::Tensor<double, 4> tensor);  // FIXME: this function is never used and not even implemented

} // namespace doci

#endif // DOCI_DOCI_UTILITY_HPP
