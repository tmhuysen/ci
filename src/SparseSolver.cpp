#include "SparseSolver.hpp"

#include <SymEigsSolver.h>
#include <MatOp/SparseSymMatProd.h>



namespace ci {
namespace solver {


/*
 *  CONSTRUCTORS
 */

/**
 *   Constructor based on the dimension @param dim of the Hamiltonian matrix representation.
 */
SparseSolver::SparseSolver(size_t dim) :
    BaseMatrixSolver(dim),
    hamiltonian (Eigen::SparseMatrix<double> (this->dim, this->dim))  // Eigen::Sparse is always initiated to zeros
{}



/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  Solve the sparse eigenvalue problem of the Hamiltonian matrix.
 */
void SparseSolver::solve() {

    // Solve the Sparse eigenvalue problem of the Hamiltonian matrix.
    Spectra::SparseSymMatProd<double> matrixVectorProduct (this->hamiltonian);
    Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, Spectra::SparseSymMatProd<double>> sparse_eigensolver (&matrixVectorProduct, 1, 128);
    sparse_eigensolver.init();
    sparse_eigensolver.compute();

    // Set the eigenvalue and eigenvector as the lowest-energy eigenpair. We can use index 0 because
    // we have specified Spectra::SMALLEST_ALGE, which selects eigenvalues with smallest algebraic value. Furthermore,
    // we have only requested 1 eigenpair.
    if (sparse_eigensolver.info() == Spectra::SUCCESSFUL) {
        this->is_solved = true;
        this->eigenvalue = sparse_eigensolver.eigenvalues()(0);
        this->eigenvector = sparse_eigensolver.eigenvectors(0);
    }
}


/**
 *  Add @param value to the matrix representation of the Hamiltonian at (@param index1, @param index2).
 */
void SparseSolver::addToMatrix(double value, size_t index1, size_t index2) {
    this->hamiltonian.coeffRef(index1,index2) += value;
}


}  // namespace solver
}  // namespace ci
