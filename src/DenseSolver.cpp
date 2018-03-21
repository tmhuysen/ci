#include "DenseSolver.hpp"



namespace ci {
namespace solver {


/*
 *  CONSTRUCTORS
 */

/**
 *   Constructor based on the dimension @param dim of the Hamiltonian matrix representation.
 */
DenseSolver::DenseSolver(size_t dim) :
    BaseMatrixSolver(dim),
    hamiltonian (Eigen::MatrixXd::Zero(this->dim, this->dim))
{}



/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  Solve the full dense eigenvalue problem of the Hamiltonian matrix.
 */
void DenseSolver::solve() {

    // Solve the dense eigenvalue problem of the Hamiltonian matrix.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> self_adjoint_eigensolver (this->hamiltonian);


    // Set the eigenvalue and eigenvector as the lowest-energy eigenpair. We can use index 0 because
    // SelfAdjointEigenSolver gives the eigenvalues (and corresponding eigenvalues) in ascending order.
    this->is_solved = true;
    this->eigenvalue = self_adjoint_eigensolver.eigenvalues()(0);
    this->eigenvector = self_adjoint_eigensolver.eigenvectors().col(0);
}


/**
 *  Add @param value to the matrix representation of the Hamiltonian at (@param index1, @param index2).
 */
void DenseSolver::addToMatrix(double value, size_t index1, size_t index2) {

    this->hamiltonian(index1, index2) += value;
}


}  // namespace solver
}  // namespace ci
