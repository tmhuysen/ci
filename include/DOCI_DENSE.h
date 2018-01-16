//
// Created by Wulfix on 14/12/2017.
//

#ifndef DOCIPROJECT_DOCI_DENSE_H
#define DOCIPROJECT_DOCI_DENSE_H


#include <DOCI.h>

class DOCI_DENSE : public DOCI {
public:
    DOCI_DENSE(StaticWrapper &calculator);
private:
    Eigen::MatrixXd hamiltonian;
protected:
    void addToHamiltonian(double value, unsigned long index1, unsigned long index2) override;

public:
    Eigen::MatrixXd getHam();
    void print() override;
};


#endif //DOCIPROJECT_DOCI_DENSE_H
