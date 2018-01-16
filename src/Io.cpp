//
// Created by Wulfix on 21/12/2017.
//

#include <DataWrapper.h>
#include <Io.h>

DataBasis read_hopping_matrix_from_file(const std::string &filename, char delimiter_char) {
    Eigen::MatrixXd one_int(18, 18);
    Eigen::Tensor<double, 4> two_int(18,18,18,18);
    std::ifstream is (filename);

}
