//
// Created by Wulfix on 14/12/2017.
//

#ifndef DOCIPROJECT_DOCI_DENSE_H
#define DOCIPROJECT_DOCI_DENSE_H


#include <DOCI.hpp>

class DOCI_DENSE : public DOCI {
public:
    DOCI_DENSE(CI_basis basis);
private:
    Eigen::MatrixXd hamiltonian;
protected:
    void addToHamiltonian(double value, size_t index1, size_t index2) override;

public:
    Eigen::MatrixXd getHam();
    void print() override;
};


#endif //DOCIPROJECT_DOCI_DENSE_H
