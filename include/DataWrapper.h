//
// Created by Wulfix on 21/12/2017.
//

#ifndef DOCI_HEAD_DATAWRAPPER_H
#define DOCI_HEAD_DATAWRAPPER_H

#include <StaticWrapper.h>
#include <hf.hpp>
#include <libwrp.hpp>

struct DataBasis {
    Eigen::MatrixXd one_int;
    Eigen::Tensor<double, 4> two_int;
    double nuc;
    unsigned long n_bf;
    unsigned long n_electrons;
};

class DataWrapper : public StaticWrapper {
private:
DataBasis energies;
public:
explicit DataWrapper(DataBasis energies1);
double calculateOverlap(unsigned long site1, unsigned long site2) override;
double calculateOverlap(unsigned long site1, unsigned long site2, unsigned long site3, unsigned long site4) override;
unsigned long getN_bf() override;
unsigned long getN_electrons() override;
};

#endif //DOCI_HEAD_DATAWRAPPER_H
