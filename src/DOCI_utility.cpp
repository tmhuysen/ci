//
// Created by Wulfix on 19/01/2018.
//

#include "DOCI_utility.hpp"

CI_basis rhf_to_CI_basis(hf::rhf::RHF &rhf_basis) {
    Eigen::MatrixXd SO_V = libwrp::transform_AO_to_SO(rhf_basis.basis.V,rhf_basis.C_canonical);
    Eigen::MatrixXd SO_T = libwrp::transform_AO_to_SO(rhf_basis.basis.T,rhf_basis.C_canonical);
    Eigen::Tensor<double, 4> SO_tei = libwrp::transform_AO_to_SO(rhf_basis.basis.tei,rhf_basis.C_canonical);

    Eigen::MatrixXd one_int = SO_V + SO_T;
    double nuc = rhf_basis.basis.molecule.internuclear_repulsion();
    auto nbf = rhf_basis.basis.nbf();
    auto nec = rhf_basis.basis.molecule.nelec;

    return CI_basis {one_int,SO_tei,nuc,nbf,nec};
}
