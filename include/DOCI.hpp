#ifndef DOCIPROJECT_DOCI_H
#define DOCIPROJECT_DOCI_H



#include <iostream>
#include "AddressingMatrix.hpp"
#include <Eigen/Dense>
#include "Extras.hpp"
#include "DOCI_utility.hpp"


class DOCI {
public:
    explicit DOCI(CI_basis ciBasis);
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


#endif //DOCIPROJECT_DOCI_H
