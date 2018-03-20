#ifndef CI_DAVIDSONSOLVER_HPP
#define CI_DAVIDSONSOLVER_HPP



#include <numopt.hpp>

#include "BaseSolver.hpp"



namespace ci {
namespace solver {

/**
 *  A wrapper class around numopt::DavidsonSolver.
 *
 *  The reason a wrapper class is needed is because this class is derived from BaseSolver and therefore leads to cleaner implementation.
 *  See for example BaseCI, FCI and DOCI: their functions constructHamiltonian accept pointers to
 */
class DavidsonSolver : ci::solver::BaseSolver {
private:
    numopt::DavidsonSolver numopt_davidson_solver;

public:
    /*
     *  CONSTRUCTORS
     */
    DavidsonSolver(const numopt::VectorFunction& matrixVectorProduct, const Eigen::VectorXd& t_0, const Eigen::VectorXd& diagonal, double residue_tolerance = 1.0e-08, double correction_threshold = 1.0e-03, size_t maximum_subspace_dimension = 15);

    /*
     *  DESTRUCTOR
     */
    ~DavidsonSolver() override = default;


    /*
     *  OVERRIDDEN PUBLIC FUNCTIONS
     */
    /**
     *  Find the lowest energy eigenpair using Davidson's algorithm.
     */
    void solve() override;


    /*
     *  GETTERS
     */
    double get_eigenvalue() const { return this->numopt_davidson_solver.get_eigenvalue(); }
    Eigen::VectorXd get_eigenvector() const { return this->numopt_davidson_solver.get_eigenvector(); }


};





}  // namespace solver
}  // namespace ci





#endif  // CI_DAVIDSONSOLVER_HPP
