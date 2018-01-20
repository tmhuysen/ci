#include "DOCI_utility.hpp"

CI_basis rhf_to_CI_basis(hf::rhf::RHF &rhf_basis) {
    Eigen::MatrixXd SO_V = libwint::transform_AO_to_SO(rhf_basis.basis.V,rhf_basis.C_canonical);
    Eigen::MatrixXd SO_T = libwint::transform_AO_to_SO(rhf_basis.basis.T,rhf_basis.C_canonical);
    Eigen::Tensor<double, 4> SO_tei = libwint::transform_AO_to_SO(rhf_basis.basis.tei,rhf_basis.C_canonical);

    Eigen::MatrixXd one_int = SO_V + SO_T;
    double nuc = rhf_basis.basis.molecule.internuclear_repulsion();
    auto nbf = rhf_basis.basis.nbf();
    auto nec = rhf_basis.basis.molecule.nelec;

    return CI_basis {one_int,SO_tei,nuc,nbf,nec};
}

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
