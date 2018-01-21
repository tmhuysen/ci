#include "CI_basis.hpp"

/** Constructor based on a given RHF instance
 */
doci::CI_basis::CI_basis(hf::rhf::RHF& rhf) {
    // Transform the one- and two-electron integrals from the AO basis to the SO basis
    Eigen::MatrixXd SO_V = libwint::transform_AO_to_SO(rhf.basis.V,rhf.C_canonical);
    Eigen::MatrixXd SO_T = libwint::transform_AO_to_SO(rhf.basis.T,rhf.C_canonical);
    this->one_ints = SO_V + SO_T;
    this->two_ints = libwint::transform_AO_to_SO(rhf.basis.tei,rhf.C_canonical);

    // Set some more parameters
    this->internuclear_repulsion = rhf.basis.molecule.internuclear_repulsion();
    this->nbf = rhf.basis.nbf();
    this->nelec = rhf.basis.molecule.nelec;
}

/** Constructor based on a given filename
 */
doci::CI_basis::CI_basis(const std::string &filename) {
    Eigen::Tensor<double, 4> tei;
    tei.setZero();
    double nuc = 0;
    size_t nbf = 0;
    size_t nelec = 0;
    Eigen::MatrixXd one_int = Eigen::MatrixXd::Zero (nbf,nbf);

    std::ifstream is (filename);

    if (is.is_open()) {
        // The dimension of the resulting matrix is found as the first entry in the file name
        std::string startline;
        std::getline(is, startline);
        std::stringstream linestream(startline);
    }

    // FIXME: this function doesn't seem to be finished yet
}