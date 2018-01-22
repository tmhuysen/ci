#include "State.hpp"

/** Constructor based on a given eigenvalue and corresponding eigenvector
 */
State::State(double eval, Eigen::VectorXd evec) : eval(eval), evec(evec) {}


bool State::operator<(const State& rhs) {
    return (this->eval < rhs.eval);
}

bool State::operator==(const State& rhs) {
    // FIXME: this implementation only checks eigenvalues, should it also check eigenvectors?

    double precision = 10000000; //

    double ELIPSON = (this->eval > rhs.eval) ?  rhs.eval/precision : this->eval/precision;

    return fabs(this->eval - rhs.eval) < fabs(ELIPSON);
}
