#ifndef CI_DENSEDOCI_HPP
#define CI_DENSEDOCI_HPP

#include "DOCI_Class.hpp"


namespace doci {
/**
 * Dense DOCI for calculations where the hamiltonian is stored in a dense matrix from the eigen lib
 * High memory requirements but fast diagonalization.
 */
class DenseDOCI : public doci::DOCI {


private:
    Eigen::MatrixXd hamiltonian;


// Protected methods
protected:
    /**
     * function that stores a calculated value in the Hamiltonian
     */
    void addToHamiltonian(double value, size_t index1, size_t index2) override;

// Public methods
public:
    /** Constructor based on a given CI_basis
     * Applies the base DOCI_class constructor and calls the DOCI calculation and
     * solves the eigenvalues of the hamiltonian with the EigenSolver.
     */
    DenseDOCI(doci::CI_basis ciBasis);


    /**
     * Helper function for printing the hamiltonian to the console
     */
    void print() override;

    /**
     * Getters
     */
    Eigen::MatrixXd getHamiltonian();
};

} // namespace doci

#endif // CI_DENSEDOCI_HPP
