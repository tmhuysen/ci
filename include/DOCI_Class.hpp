#ifndef DOCI_DOCI_HPP
#define DOCI_DOCI_HPP


#include <iostream>
#include "AddressingMatrix.hpp"
#include <Eigen/Dense>
#include "Extras.hpp"
#include "DOCI_utility.hpp"


class DOCI_Class {
public:
    explicit DOCI_Class(CI_basis ciBasis);
    const std::vector<State> &getGroundstates() const;
    //Virtuals
    virtual void print()=0;
protected: //variables
    unsigned long sites;
    unsigned long electrons;
    unsigned long nbf;
    CI_basis basis;
    AddressingMatrix ad_mat;
    std::vector<State> groundstates;
    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;

protected://methods
    void calculateDoci(double start, double end);
    void groundStates(State state);
    //Virtuals
    virtual void addToHamiltonian(double value, size_t index1, size_t index2)=0;
};


#endif // DOCI_DOCI_HPP
