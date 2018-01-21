#include "utility.hpp"

CI_basis file_to_CI_basis(const std::string &filename) {
    Eigen::Tensor<double, 4> tei;
    tei.setZero();
    double nuc;
    size_t nbf;
    size_t nec;
    Eigen::MatrixXd one_int;
    one_int = Eigen::MatrixXd::Zero(nbf,nbf);

    std::ifstream is (filename);

    if (is.is_open()) {
        // The dimension of the resulting matrix is found as the first entry in the file name
        std::string startline;
        std::getline(is, startline);
        std::stringstream linestream(startline);
    }


    return CI_basis {one_int,tei,nuc,nbf,nec};
}
