#include "BaseCI.hpp"

#include <numopt.hpp>

#include "SolverType.hpp"
#include "DenseSolver.hpp"
#include "SparseSolver.hpp"



namespace ci {


/*
 *  PUBLIC METHODS
 */

/**
 *  Find the lowest energy eigenpair of the Hamiltonian.
 */
void BaseCI::solve(ci::solver::SolverType solver_type, const bmqc::AddressingScheme& addressing_scheme) {

    switch (solver_type) {

        case ci::solver::SolverType::DENSE:
            auto* dense_solver = new ci::solver::DenseSolver(this->dim);

        case ci::solver::SolverType::SPARSE:
            auto* sparse_solver = new ci::solver::SparseSolver(this->dim);

        case ci::solver::SolverType::DAVIDSON:
            numopt::DavidsonSolver
    }

}


}  // namespace ci
