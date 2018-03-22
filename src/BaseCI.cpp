#include "BaseCI.hpp"

#include <numopt.hpp>

#include "DenseSolver.hpp"
#include "SparseSolver.hpp"



namespace ci {



/*
 *  PROTECTED CONSTRUCTORS
 */

/**
 *  Protected constructor to initialize the reference @member so_basis by @param so_basis.
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
void BaseCI::solveMatrixEigenvalueProblem() {

    // Check if we gave the right pointer: it should be derived from BaseMatrixSolver
    assert(dynamic_cast<numopt::eigenproblem::BaseMatrixSolver*>(this->eigensolver_ptr));

    // Initialize the Hamiltonian matrix and solve the eigenvalue problem associated to it.
    this->constructHamiltonian();
    this->eigensolver_ptr->solve();
}



/*
 *  DESTRUCTOR
 */

BaseCI::~BaseCI() {
    delete[] this->eigensolver_ptr;
}



/*
 *  PUBLIC METHODS
 */

/**
 *  Find the lowest energy eigenpair of the Hamiltonian, using a @param solver_type.
 */
void BaseCI::solve(numopt::eigenproblem::SolverType solver_type) {

    // Depending on how the user wants to solve the eigenvalue problem, construct the appropriate solver
    switch (solver_type) {

        case numopt::eigenproblem::SolverType::DENSE: {
            this->eigensolver_ptr = new numopt::eigenproblem::DenseSolver(this->dim);
            this->solveMatrixEigenvalueProblem();
        }

        case numopt::eigenproblem::SolverType::SPARSE: {
            this->eigensolver_ptr= new numopt::eigenproblem::SparseSolver(this->dim);
            this->solveMatrixEigenvalueProblem();
        }

        case numopt::eigenproblem::SolverType::DAVIDSON: {
            numopt::VectorFunction matrixVectorProduct = [this] (const Eigen::VectorXd& x) { return this->matrixVectorProduct(x); };
            Eigen::VectorXd diagonal = this->calculateDiagonal();

            Eigen::VectorXd t_0 = Eigen::VectorXd::Zero(this->dim);
            t_0(this->dim) = 1;  // in reverse lexical notation, the Hartree-Fock determinant has the highest address

            this->eigensolver_ptr = new numopt::eigenproblem::DavidsonSolver(matrixVectorProduct, t_0, diagonal);
            this->eigensolver_ptr->solve();
        }

    }

}


}  // namespace ci
