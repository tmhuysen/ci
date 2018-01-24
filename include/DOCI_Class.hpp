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
     * calculate all hamiltonian elements.
     * @param start : for parallellization
     * @param end
     */
    void calculateDoci(double start, double end);  // FIXME: add comments/documentation
    void groundStates(doci::State state);  // FIXME:: add comments/documentation

    // Virtuals
    virtual void addToHamiltonian(double value, size_t index1, size_t index2)=0;  // FIXME: add comments/documentation


public:
    /** Constructor based on a given CI_basis
     */
    explicit DOCI(doci::CI_basis ciBasis);  // FIXME: add comments/documentation
    const std::vector<doci::State>& getGroundstates() const;  // FIXME: add comments/documentation

    // Virtuals
    virtual void print()=0;  // FIXME: add comments/documentation
};

} // namespace doci

#endif // DOCI_DOCI_HPP
