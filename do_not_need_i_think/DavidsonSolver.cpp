#include "DavidsonSolver.hpp"



namespace ci {
namespace solver {


/*
 *  CONSTRUCTORS
 */
DavidsonSolver::DavidsonSolver(const numopt::VectorFunction& matrixVectorProduct, const Eigen::VectorXd& t_0, const Eigen::VectorXd& diagonal, double residue_tolerance, double correction_threshold, size_t maximum_subspace_dimension) :
    numopt_davidson_solver (numopt::DavidsonSolver(matrixVectorProduct, t_0, diagonal, residue_tolerance, correction_threshold, maximum_subspace_dimension))
{}




/*
 *  OVERRIDDEN PUBLIC FUNCTIONS
 */
/**
 *  Find the lowest energy eigenpair using Davidson's algorithm.
 */
DavidsonSolver::solve() {

    // Since this class is
}





}  // namespace solver
}  // namespace ci
