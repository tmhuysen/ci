#ifndef DOCI_DOCI_HPP
#define DOCI_DOCI_HPP

#include "AddressingMatrix.hpp"
#include "CI_basis.hpp"
#include "State.hpp"

#include <iostream>
#include <Eigen/Dense>


class DOCI {

// Protected variables
protected:
    size_t K;  // number of spatial orbitals
    size_t npairs;
    size_t nbf;

    AddressingMatrix ad_mat;

    std::vector<State> groundstates;

    doci::CI_basis basis;

    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;


// Protected methods
protected:
    void calculateDoci(double start, double end);  // FIXME: add comments/documentation
    void groundStates(State state);  // FIXME:: add comments/documentation

    // Virtuals
    virtual void addToHamiltonian(double value, size_t index1, size_t index2)=0;  // FIXME: add comments/documentation


public:
    /** Constructor based on a given CI_basis
     */
    explicit DOCI(doci::CI_basis ciBasis);  // FIXME: add comments/documentation
    const std::vector<State> &getGroundstates() const;  // FIXME: add comments/documentation

    // Virtuals
    virtual void print()=0;  // FIXME: add comments/documentation
};


#endif // DOCI_DOCI_HPP
