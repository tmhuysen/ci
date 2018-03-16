#include "State.hpp"

<<<<<<< HEAD


namespace ci {

=======
/** Constructor based on a given eigenvalue and corresponding eigenvector
 */
doci::State::State(double eval, Eigen::VectorXd evec) : eigenvalue(eval), eigenvector(evec) {}


bool doci::State::operator<(const doci::State& rhs) {
    return (this->eigenvalue < rhs.eigenvalue);
}

bool doci::State::operator>(const doci::State& rhs) {
    return (this->eigenvalue > rhs.eigenvalue);
}
>>>>>>> a1be18ed473b2634f7d4cca3483eb610e27257b0

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


<<<<<<< HEAD

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
=======
    double ELIPSON = (this->eigenvalue > rhs.eigenvalue) ?  rhs.eigenvalue/precision : this->eigenvalue/precision;

    return fabs(this->eigenvalue - rhs.eigenvalue) < fabs(ELIPSON);
}

double doci::State::get_eigenvalue() const {
    return eigenvalue;
}

const Eigen::VectorXd &doci::State::get_eigenvector() const {
    return eigenvector;
>>>>>>> a1be18ed473b2634f7d4cca3483eb610e27257b0
}


}  // namespace ci



// doci::State::State() : doci::State(std::numeric_limits<double>::max(), Eigen::VectorXd()){};
