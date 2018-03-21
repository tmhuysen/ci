#include "DenseHamiltonianSolver.hpp"



namespace ci {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param dimension.
 */
DenseHamiltonianSolver::DenseHamiltonianSolver(size_t dimension):
    dimension (dimension),
    hamiltonian (Eigen::MatrixXd::Zero(this->dimension, this->dimension))
{}



/*
 *  PUBLIC OVERRIDDEN VIRTUAL METHODS
 */

/**
 *  Solves the eigenvalue problem of the Hamiltonian.
 *
 *  For this derived class, this means diagonalizing the full dense Hamiltonian.
 */
void DenseHamiltonianSolver::solve() override {

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> dense_eigensolver (this->hamiltonian);
    this->eigenvalues = dense_eigensolver.eigenvalues();  // eigenvalues and their associated eigenvectors are stored in ascending order
    this->eigenvectors = dense_eigensolver.eigenvectors();

    this->is_solved = true;


    // Extract the ground state (and its degeneracies) into this->ground_states
    this->addToGroundStates(ci::State(eigenvalues(0), eigenvectors.col(0)));  // add the first eigenvector anyways
    for (size_t i = 1; i < this->dimension; i++) {
        ci::State state (eigenvalues(i), eigenvectors.col(i));

        if (!this->addToGroundStates(state)) {  // if the state isn't added to this->ground_states, we can quite since eigenvalues are stored in ascending order
            break;
        }
    }
}


/**
 *  Adds @param value to the current value of the Hamiltonian matrix at position (@param index1, @param index2).
 */
void DenseHamiltonianSolver::addElement(double value, size_t index1, size_t index2) override {
    this->hamiltonian(index1,index2) += value;
}


}  // namespace ci
