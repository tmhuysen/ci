#ifndef CI_UTILITY_HPP
#define CI_UTILITY_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>


#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <Eigen/Sparse>

namespace doci {

void symmatu(Eigen::MatrixXd& mat);  //fills lower triagonal of upper triagonal matrix
void symmatu_reverse(Eigen::MatrixXd& mat); //fills upper triagonal of lower triagonal matrix
void symmatu_tensor_reverse(Eigen::Tensor<double, 4> &tensor); //fills upper part of tensor based on some kind of symmetry.






/*
 * PARSING UTILITY
 */
// trim from start (in place)
void ltrim(std::string &s);
// trim from end (in place)
void rtrim(std::string &s);

} // namespace doci

#endif // CI_UTILITY_HPP
