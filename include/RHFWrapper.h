//
// Created by Wulfix on 10/12/2017.
//

#ifndef DOCIPROJECT_RHFWRAPPER_H
#define DOCIPROJECT_RHFWRAPPER_H

#include <StaticWrapper.h>
#include <hf.hpp>
#include <libwrp.hpp>
class RHFWrapper : public StaticWrapper {
private:
    hf::rhf::RHF& rhf_basis;
    Eigen::MatrixXd SO_V;
    Eigen::MatrixXd SO_T;
    Eigen::Tensor<double, 4> SO_tei;
public:
    explicit RHFWrapper(hf::rhf::RHF& rhf_basis);
    double calculateOverlap(unsigned long site1, unsigned long site2) override;
    double calculateOverlap(unsigned long site1, unsigned long site2, unsigned long site3, unsigned long site4) override;
    unsigned long getN_bf() override;
    unsigned long getN_electrons() override;
};



#endif //DOCIPROJECT_TESTWRAPPER_H
