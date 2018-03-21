#include "BaseCI.hpp"

#include <numopt.hpp>

#include "DenseSolver.hpp"
#include "SparseSolver.hpp"



namespace ci {


/*
 *  PRIVATE METHODS
 */

/**
 *  Initialize and subsequently solve the eigenvalue problem associated to the derived CI class, using a pointer
 *  to a @param matrix_solver and to speed up the searches for the addresses of spin strings using a @param
 *  addressing_scheme
 */
void BaseCI::solveMatrixEigenvalueProblem(ci::solver::BaseMatrixSolver* matrix_solver, const bmqc::AddressingScheme& addressing_scheme) {

    // Initialize the Hamiltonian matrix and solve the eigenvalue problem associated to it.
    this->constructHamiltonian(matrix_solver, addressing_scheme);
    matrix_solver->solve();


    // Set the eigenvalues and eigenvectors
    this->is_solved = true;
    this->eigenvalue = matrix_solver->get_eigenvalue();
    this->eigenvector = matrix_solver->get_eigenvector();
}



/*
 *  PUBLIC METHODS
 */

/**
 *  Find the lowest energy eigenpair of the Hamiltonian.
 */
void BaseCI::solve(ci::solver::SolverType solver_type, const bmqc::AddressingScheme& addressing_scheme) {

    switch (solver_type) {

        case ci::solver::SolverType::DENSE: {
            auto *dense_solver = new ci::solver::DenseSolver(this->dim);
            this->solveMatrixEigenvalueProblem(dense_solver, addressing_scheme);
        }

        case ci::solver::SolverType::SPARSE: {
            auto *sparse_solver = new ci::solver::SparseSolver(this->dim);
            this->solveMatrixEigenvalueProblem(sparse_solver, addressing_scheme);
        }

        case ci::solver::SolverType::DAVIDSON: {
            numopt::VectorFunction matrixVectorProduct = [this, addressing_scheme] (const Eigen::VectorXd& x) { return this->matrixVectorProduct(addressing_scheme, x); };
            Eigen::VectorXd diagonal = this->calculateDiagonal();

            Eigen::VectorXd t_0 = Eigen::VectorXd::Zero(this->dim);
            t_0(this->dim) = 1;  // in reverse lexical notation, the Hartree-Fock determinant has the highest address

            numopt::DavidsonSolver davidson_solver (matrixVectorProduct, t_0, diagonal);

            // Set the eigenvalues and eigenvectors
            this->is_solved = true;
            this->eigenvalue = davidson_solver.get_eigenvalue();
            this->eigenvector = davidson_solver.get_eigenvector();
        }

    }

}



/*
 *  GETTERS
 */

double BaseCI::get_eigenvalue () const {

    if (this->is_solved) {
        return this->eigenvalue;
    } else {
        throw std::runtime_error("The eigenvalue problem hasn't been solved yet and you are trying to get the eigenvalue.");
    }
}

Eigen::VectorXd BaseCI::get_eigenvector () const {

    if (this->is_solved) {
        return this->eigenvector;
    } else {
        throw std::runtime_error("The eigenvalue problem hasn't been solved yet and you are trying to get the eigenvector.");
    }
}


}  // namespace ci
