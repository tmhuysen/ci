#include "BaseMatrixSolver.hpp"



namespace ci {
namespace solver {


/*
 *  GETTERS
 */
double BaseMatrixSolver::get_eigenvalue() const {

    if (this->is_solved) {
        return this->eigenvalue;
    } else {
        throw std::runtime_error("The eigenvalue problem hasn't been solved yet and you are trying to get the eigenvalue.");
    }
}

Eigen::VectorXd BaseMatrixSolver::get_eigenvector() const {

    if (this->is_solved) {
        return this->eigenvector;
    } else {
        throw std::runtime_error("The eigenvalue problem hasn't been solved yet and you are trying to get the eigenvector.");
    }
}


}  // namespace solver
}  // namespace ci
