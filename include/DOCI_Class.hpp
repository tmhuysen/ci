#ifndef DOCI_DOCI_HPP
#define DOCI_DOCI_HPP

#include "AddressingMatrix.hpp"
#include "CI_basis.hpp"
#include "State.hpp"

#include <iostream>
#include <Eigen/Dense>


namespace doci {

class DOCI {

// Protected variables
protected:
    size_t K;  // number of spatial orbitals
    size_t npairs; // number of electron pairs
    size_t nbf; // number of basis functions

    AddressingMatrix ad_mat; // FIXME: scheme

    std::vector<doci::State> groundstates; //vector of the all states with the ground energy.

    doci::CI_basis basis; //contains all information required to do a DOCI calculation

    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;


// Protected methods
protected:
    /**
     * calculate hamiltonian elements (the lower triagonal).
     * @param start,end : (for parallellization?) calculates only the fraction between start and end
     * example: start:O.5 to end:0.75. currently excludes based on nbf (iterates over fraction of bf)
     */
    void calculateDoci(double start, double end);

    /**
     * Adds a @param state to the groundstates vector of our DOCI,
     * if this state's eigenvalue is equal to the eigenvalue of current groundstates, then it is added.
     * if this state's eigenvalue is lower than it replaced the current groundstates.
     */
    void groundStates(doci::State state);

    // Virtuals
    /**
     * Adds a value as an element of the hamiltonian matrix.
     * Virtual: when implementing this class one can opt for many ways to represent the hamiltonian.
     */
    virtual void addToHamiltonian(double value, size_t index1, size_t index2)=0;


public:

    /** Constructor based on a given CI_basis
     * Initializes the base variables of the class, and calculates the number of basis functions (nbf)
     */
    explicit DOCI(doci::CI_basis ciBasis);


    /**
     * Getters
     */
    const std::vector<doci::State>& getGroundstates() const;

    // Virtuals
    /**
     * Virtual print function, ment for printing the representation of the hamiltonian.
     */
    virtual void print()=0;
};

} // namespace doci

#endif // DOCI_DOCI_HPP
