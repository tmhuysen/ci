#ifndef DOCI_SPARSEDOCI_HPP
#define DOCI_SPARSEDOCI_HPP




#include "DOCI_Class.hpp"
#include <Eigen/Sparse>
#include <SymEigsSolver.h>
#include <MatOp/SparseSymMatProd.h>


namespace doci {
/**
 * Sparse DOCI for calculations where the hamiltonian is stored in a sparse matrix from the eigen lib
 * Lower requirements but slower diagonalization.
 */
class SparseDOCI : public doci::DOCI {


private:
    Eigen::SparseMatrix<double> hamiltonian;


// Protected methods
protected:
    /**
     * function that stores a calculated value in the Hamiltonian
     */
    void addToHamiltonian(double value, size_t index1, size_t index2) override;


public:
    /** Constructor based on a given CI_basis
     * Applies the base DOCI_class constructor and calls the DOCI calculation and
         * solves the eigenvalues of the hamiltonian with the EigenSolver.
     */
    SparseDOCI(doci::CI_basis ciBasis);


    /**
     * Getters
     */
    Eigen::SparseMatrix<double> getHamiltonian();
    /**
     * Helper function for printing the hamiltonian to the console
     */
    void print() override;
};

} // namespace doci



#endif //DOCI_SPARSEDOCI_HPP
