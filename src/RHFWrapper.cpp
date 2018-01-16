//
// Created by Wulfix on 10/12/2017.
//

#include "RHFWrapper.h"

RHFWrapper::RHFWrapper(hf::rhf::RHF &rhf_basis):rhf_basis(rhf_basis) {
    SO_V = libwrp::transform_AO_to_SO(rhf_basis.basis.V,rhf_basis.C_canonical);
    SO_T = libwrp::transform_AO_to_SO(rhf_basis.basis.T,rhf_basis.C_canonical);
    SO_tei = libwrp::transform_AO_to_SO(rhf_basis.basis.tei,rhf_basis.C_canonical);




}

double RHFWrapper::calculateOverlap(unsigned long site1, unsigned long site2) {
    double kin = SO_T(site1,site2);
    double nuc = SO_V(site1,site2);
    return kin+nuc;
}

double RHFWrapper::calculateOverlap(unsigned long site1, unsigned long site2, unsigned long site3, unsigned long site4) {
    return SO_tei(site1,site2,site3,site4);
}

unsigned long RHFWrapper::getN_bf() {
    return rhf_basis.basis.nbf();
};

unsigned long RHFWrapper::getN_electrons() {
    return rhf_basis.basis.molecule.nelec;
}