#include "BaseCI.hpp"

#include <numopt.hpp>

#include "DenseSolver.hpp"
#include "SparseSolver.hpp"



namespace ci {



/*
 *  PROTECTED CONSTRUCTORS
 */

/**
 *  Protected constructor to initialize the reference @member so_basis by @param so_basis
 */
BaseCI::BaseCI(libwint::SOBasis& so_basis) :
    so_basis (so_basis)
{}



/*
 *  PROTECTED METHODS
 */

/**
 *  Initialize and subsequently solve the eigenvalue problem associated to the derived CI class, using a pointer
 *  to a @param matrix_solver.
 */
void BaseCI::solveMatrixEigenvalueProblem(ci::solver::BaseMatrixSolver* matrix_solver) {

    // Initialize the Hamiltonian matrix and solve the eigenvalue problem associated to it.
    this->constructHamiltonian(matrix_solver);
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
 *  Find the lowest energy eigenpair of the Hamiltonian, using a @param solver_type.
 */
void BaseCI::solve(ci::solver::SolverType solver_type) {

    switch (solver_type) {

        case ci::solver::SolverType::DENSE: {
            auto *dense_solver = new ci::solver::DenseSolver(this->dim);
            this->solveMatrixEigenvalueProblem(dense_solver);
        }

        case ci::solver::SolverType::SPARSE: {
            auto *sparse_solver = new ci::solver::SparseSolver(this->dim);
            this->solveMatrixEigenvalueProblem(sparse_solver);
        }

        case ci::solver::SolverType::DAVIDSON: {
            numopt::VectorFunction matrixVectorProduct = [this] (const Eigen::VectorXd& x) { return this->matrixVectorProduct(x); };
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
