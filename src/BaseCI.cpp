#include "BaseCI.hpp"

#include <numopt.hpp>
#include <iomanip>
#include <limits>
#include <DavidsonSolverLemmens.hpp>

#include "DenseSolver.hpp"
#include "SparseSolver.hpp"
#include "ArmaDenseSolver.hpp"



namespace ci {



/*
 *  PROTECTED CONSTRUCTORS
 */

/**
 *  Protected constructor given a @param so_basis and a dimension @dim.
 */
BaseCI::BaseCI(libwint::SOMullikenBasis &so_basis, size_t dim) :
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
            auto dense_solver = new numopt::eigenproblem::DenseSolver(this->dim);
            this->constructHamiltonian(dense_solver);

            Eigen::VectorXd t_0 = Eigen::VectorXd::Zero(this->dim);
            t_0(0) = 1; //  lexical notation, the Hartree-Fock determinant has the highest address
            this->eigensolver_ptr = new numopt::eigenproblem::DavidsonSolver(dense_solver->get_matrix(), t_0);


            this->eigensolver_ptr->solve();
            break;
        }

        case numopt::eigenproblem::SolverType::DAVIDSONLEMMENS: {
            auto dense_solver = new numopt::eigenproblem::DenseSolver(this->dim);
            this->constructHamiltonian(dense_solver);

            Eigen::VectorXd t_0 = Eigen::VectorXd::Zero(this->dim);
            t_0(0) = 1; //  lexical notation, the Hartree-Fock determinant has the highest address
            this->eigensolver_ptr = new numopt::eigenproblem::DavidsonSolverLemmens(dense_solver->get_matrix(), t_0);


            this->eigensolver_ptr->solve();
            break;
        }

        case numopt::eigenproblem::SolverType::ARMADENSE: {
            auto dense_solver = new numopt::eigenproblem::ArmaDenseSolver(this->dim);
            this->solveMatrixEigenvalueProblem(dense_solver);
            this->eigensolver_ptr = dense_solver;  // prevent data from going out of scope
            // we are only assigning this->eigensolver_ptr now, because
            // this->solveMatrixEigenvalueProblem only accepts BaseMatrixSolver*
            break;
        }



    }

}

/**
 *  Solves the eigenvalue problem with a mulliken constraint and returns the energy.
 */
double BaseCI::findConstrained(numopt::eigenproblem::SolverType solver_type, std::vector<size_t> AO_set,
                               double constraint){

    //init the procedure
    this->so_basis.calculateMullikenMatrix(AO_set);
    this->solve(solver_type);
    this->compute1RDM();
    double threshold = 1e-6;
    double population = so_basis.mullikenPopulationCI(this->one_rdm_aa,this->one_rdm_bb);
    double error = std::abs(population - constraint);
    // Iteration variables
    size_t iterations = 0;
    size_t max_iterations = 4;
    double lagrange_multiplier = 0;
    double range = 3;
    double max_bounds = lagrange_multiplier + range;
    double min_bounds = lagrange_multiplier - range;
    double interval = range/10;

    while(error > threshold && iterations<max_iterations){
        iterations++;
        for(double i = max_bounds; i>min_bounds;i -= interval) {
            so_basis.set_lagrange_multiplier(i);
            this->solve(solver_type);
            this->compute1RDM();
            double new_population = so_basis.mullikenPopulationCI(this->one_rdm_aa, this->one_rdm_bb);
            double new_error = std::abs(new_population - constraint);
            if (new_error < error) {
                lagrange_multiplier = i;
                error = new_error;
                population = new_population;
            }

        }
        range = range/10;
        max_bounds = lagrange_multiplier + range;
        min_bounds = lagrange_multiplier - range;
        interval = range/10;

    }
    so_basis.set_lagrange_multiplier(lagrange_multiplier);
    this->solve(solver_type);
    this->compute1RDM();
    double energy = this->get_eigenvalue()+lagrange_multiplier*population;
    if(error > threshold){
        std::cout<<std::endl<<" WARNING THE ERROR = "<<error<<std::endl;
        std::cout<<std::endl<<" with multiplier = "<<lagrange_multiplier<<std::endl;
    }
    return energy;
}

/**
 *  Solves the eigenvalue problem with a for a langrange multiplier mulliken constraint and returns the energy.
 */
double BaseCI::solveConstrained(numopt::eigenproblem::SolverType solver_type, std::vector<size_t> AO_set, double multiplier){
    this->so_basis.calculateMullikenMatrix(AO_set);
    this->so_basis.set_lagrange_multiplier(multiplier);
    this->solve(solver_type);
    this->compute1RDM();
    this->population_set = this->so_basis.mullikenPopulationCI(this->one_rdm_aa, this->one_rdm_bb);
    double energy = this->get_eigenvalue()+multiplier*population_set;
    return energy;

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
