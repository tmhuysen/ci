#include "BaseCI.hpp"

#include <numopt.hpp>

#include "DenseSolver.hpp"
#include "SparseSolver.hpp"



namespace ci {



/*
 *  PROTECTED CONSTRUCTORS
 */

/**
 *  Protected constructor given a @param so_basis and a dimension @dim.
 */
BaseCI::BaseCI(libwint::SOBasis& so_basis, size_t dim) :
    so_basis (so_basis),
    dim (dim)
{}



/*
 *  PROTECTED METHODS
 */

/**
 *  Initialize and subsequently solve the eigenvalue problem associated to the derived CI class, using a pointer
 *  to a @param matrix_solver_ptr.
 */
void BaseCI::solveMatrixEigenvalueProblem(numopt::eigenproblem::BaseMatrixSolver* matrix_solver_ptr) {

    // Initialize the Hamiltonian matrix and solve the eigenvalue problem associated to it.
    this->constructHamiltonian(matrix_solver_ptr);
    matrix_solver_ptr->solve();
}



/*
 *  DESTRUCTOR
 */

BaseCI::~BaseCI() {
    delete this->eigensolver_ptr;
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
            auto dense_solver = new numopt::eigenproblem::DenseSolver(this->dim);
            this->solveMatrixEigenvalueProblem(dense_solver);
            this->eigensolver_ptr = dense_solver;  // prevent data from going out of scope
                                                   // we are only assigning this->eigensolver_ptr now, because
                                                   // this->solveMatrixEigenvalueProblem only accepts BaseMatrixSolver*
            break;
        }

        case numopt::eigenproblem::SolverType::SPARSE: {
            auto sparse_solver = new numopt::eigenproblem::SparseSolver(this->dim);
            this->solveMatrixEigenvalueProblem(sparse_solver);
            this->eigensolver_ptr = sparse_solver;  // prevent data from going out of scope
                                                    // we are only assigning this->eigensolver_ptr now, because
                                                    // this->solveMatrixEigenvalueProblem only accepts BaseMatrixSolver*
            break;
        }

        case numopt::eigenproblem::SolverType::DAVIDSON: {
            numopt::VectorFunction matrixVectorProduct = [this] (const Eigen::VectorXd& x) { return this->matrixVectorProduct(x); };
            Eigen::VectorXd diagonal = this->calculateDiagonal();

            Eigen::VectorXd t_0 = Eigen::VectorXd::Zero(this->dim);
            t_0(this->dim) = 1;  // in reverse lexical notation, the Hartree-Fock determinant has the highest address

            this->eigensolver_ptr = new numopt::eigenproblem::DavidsonSolver(matrixVectorProduct, t_0, diagonal);
            this->eigensolver_ptr->solve();
            break;
        }

    }

}



/*
 * GETTERS
 */

Eigen::MatrixXd BaseCI::get_one_rdm_aa() const {
    if(!this->are_computed_one_rdm){
        throw std::logic_error("The requested reduced density matrix is not computed yet");
    }
    return this->one_rdm_aa;
}


Eigen::MatrixXd BaseCI::get_one_rdm_bb() const {
    if(!this->are_computed_one_rdm){
        throw std::logic_error("The requested reduced density matrix is not computed yet");
    }
    return this-> one_rdm_bb;
}


Eigen::Tensor<double, 4> BaseCI::get_two_rdm_aaaa() const {
    if(!this->are_computed_two_rdm){
        throw std::logic_error("The requested reduced density matrix is not computed yet");
    }
    return this->two_rdm_aaaa;
}


Eigen::Tensor<double, 4> BaseCI::get_two_rdm_abba() const {
    if(!this->are_computed_two_rdm){
        throw std::logic_error("The requested reduced density matrix is not computed yet");
    }
    return this->two_rdm_abba;
}


Eigen::Tensor<double, 4> BaseCI::get_two_rdm_baab() const {
    if(!this->are_computed_two_rdm){
        throw std::logic_error("The requested reduced density matrix is not computed yet");
    }
    return this->two_rdm_baab;
}


Eigen::Tensor<double, 4> BaseCI::get_two_rdm_bbbb() const {
    if(!this->are_computed_two_rdm){
        throw std::logic_error("The requested reduced density matrix is not computed yet");
    }
    return this->two_rdm_bbbb;
}


}  // namespace ci
