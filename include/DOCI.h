//
// Created by Wulfix on 04/12/2017.
//

#ifndef DOCIPROJECT_DOCI_H
#define DOCIPROJECT_DOCI_H


#include <StaticWrapper.h>
#include <iostream>
#include "AddressingMatrix.h"
#include <Eigen/Dense>
#include "Extras.h"
struct State {
    double eigenValue;          // The energy of the solution, a.k.a. the eigenvalue
    Eigen::VectorXd eigenVector; // The coefficients of the solution with respect to the given basis, a.k.a. the eigenvector corresponding to the eigenvalue
};

class DOCI {
public:
    explicit DOCI(StaticWrapper& calculator);
    const std::vector<State> &getGroundstates() const;
    //Virtuals
    virtual void print()=0;
protected:
    unsigned long sites;
    unsigned long electrons;
    unsigned long nbf;
    AddressingMatrix ad_mat;
    StaticWrapper* integralCalculator;
    std::vector<State> groundstates;

protected:
    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;
    void calculateDoci(double start, double end);
    void groundStates(State state);
    //Virtuals
    virtual void addToHamiltonian(double value, unsigned long index1, unsigned long index2)=0;
};

bool compareState(const State &o1, const State &o2);
bool areSame(const State &o1, const State &o2);

#endif //DOCIPROJECT_DOCI_H
