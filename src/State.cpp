#include "State.hpp"



namespace ci {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given eigenvalue and corresponding eigenvector
 */
State::State(double eigenvalue, const Eigen::VectorXd& eigenvector) :
    eigenvalue (eigenvalue),
    eigenvector (eigenvector)
{}



/*
 *  PUBLIC METHODS
 */

/**
  *  @return if this is degenerate to @param other, i.e. the eigenvalues are the same within the given @param tolerance.
  *
  *  @param tolerance: defaults to 1.0e-06
  */
bool State::isDegenerate(const ci::State& other, double tolerance) const {

    return std::abs(this->eigenvalue - other.eigenvalue) < tolerance;
}


/**
 *  @return if this has an eigenvalue lower than @param other, within a given @param tolerance.
 *
 *  @param tolerance: defaults to 1.0e-06
 */
bool State::hasLowerEnergy(const ci::State& other, double tolerance) const {

    return (other.eigenvalue - this->eigenvalue) > tolerance;
}


}  // namespace ci



// doci::State::State() : doci::State(std::numeric_limits<double>::max(), Eigen::VectorXd()){};
