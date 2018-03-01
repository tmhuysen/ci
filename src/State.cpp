#include "State.hpp"

/** Constructor based on a given eigenvalue and corresponding eigenvector
 */
doci::State::State(double eval, Eigen::VectorXd evec) : eval(eval), evec(evec) {}


bool doci::State::operator<(const doci::State& rhs) {
    return (this->eval < rhs.eval);
}

bool doci::State::operator>(const doci::State& rhs) {
    return (this->eval > rhs.eval);
}

bool doci::State::operator==(const doci::State& rhs) {
    // Does not check for eigenvectors, operator is meant to check for degeneracy.

    double precision = 1000000;

    double ELIPSON = (this->eval > rhs.eval) ?  rhs.eval/precision : this->eval/precision;

    return fabs(this->eval - rhs.eval) < fabs(ELIPSON);
}

double doci::State::getEval() const {
    return eval;
}

const Eigen::VectorXd &doci::State::getEvec() const {
    return evec;
}

doci::State::State() : doci::State(std::numeric_limits<double>::max(), Eigen::VectorXd()){};
