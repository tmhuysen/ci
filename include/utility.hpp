#ifndef DOCI_DOCI_UTILITY_HPP
#define DOCI_DOCI_UTILITY_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>


#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>

namespace doci {

void symmatu(Eigen::MatrixXd& mat);  // FIXME: this function is never used
void symmatu_reverse(Eigen::MatrixXd& mat);

void symmatu_tensor_reverse(Eigen::Tensor<double, 4> tensor);  // FIXME: this function is never used and not even implemented








/*
 * PARSING UTILITY
 */
// trim from start (in place)
void ltrim(std::string &s);
// trim from end (in place)
void rtrim(std::string &s);
// trim from both ends (in place)
//void trim(std::string &s);
//copies trimmed string and returns it.
//static inline std::string ltrim_copy(std::string s);
//static inline std::string rtrim_copy(std::string s);
//static inline std::string trim_copy(std::string s);
} // namespace doci

#endif // DOCI_DOCI_UTILITY_HPP
