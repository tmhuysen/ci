//
// Created by Wulfix on 04/12/2017.
//

#ifndef DOCIPROJECT_DOCI_H
#define DOCIPROJECT_DOCI_H



#include <iostream>
#include "AddressingMatrix.h"
#include <Eigen/Dense>
#include "Extras.h"
#include "DOCI_utility.hpp"


class DOCI {
public:
    explicit DOCI(CI_basis ciBasis);
    const std::vector<State> &getGroundstates() const;
    //Virtuals
    virtual void print()=0;
protected:
    unsigned long sites;
    unsigned long electrons;
    unsigned long nbf;
    CI_basis basis;
    AddressingMatrix ad_mat;
    std::vector<State> groundstates;

protected:
    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;
    void calculateDoci(double start, double end);
    void groundStates(State state);
    //Virtuals
    virtual void addToHamiltonian(double value, size_t index1, size_t index2)=0;
};


#endif //DOCIPROJECT_DOCI_H
