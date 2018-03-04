#include "State.hpp"

/** Constructor based on a given eigenvalue and corresponding eigenvector
 */
doci::State::State(double eval, Eigen::VectorXd evec) : eigenvalue(eval), eigenvector(evec) {}


bool doci::State::operator<(const doci::State& rhs) {
    return (this->eigenvalue < rhs.eigenvalue);
}

bool doci::State::operator>(const doci::State& rhs) {
    return (this->eigenvalue > rhs.eigenvalue);
}

bool doci::State::operator==(const doci::State& rhs) {
    // Does not check for eigenvectors, operator is meant to check for degeneracy.

    double precision = 1000000;

    double ELIPSON = (this->eigenvalue > rhs.eigenvalue) ?  rhs.eigenvalue/precision : this->eigenvalue/precision;

    return fabs(this->eigenvalue - rhs.eigenvalue) < fabs(ELIPSON);
}

double doci::State::get_eigenvalue() const {
    return eigenvalue;
}

const Eigen::VectorXd &doci::State::get_eigenvector() const {
    return eigenvector;
}

doci::State::State() : doci::State(std::numeric_limits<double>::max(), Eigen::VectorXd()){};
