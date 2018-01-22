#include "State.hpp"

/** Constructor based on a given eigenvalue and corresponding eigenvector
 */
State::State(double eval, Eigen::VectorXd evec) : eval(eval), evec(evec) {}



bool State::operator<(const State& l, const State& r) {
    return (l.eval < r.eval);
}

bool State::operator==(const State& l, const State& r) {
    // FIXME: this implementation only checks eigenvalues, it should also check eigenvectors

    double precision = 10000000; //

    double ELIPSON = (l.eval > r.eval) ?  r.eval/precision : l.eval/precision;

    return fabs(l.eval - r.eval) < fabs(ELIPSON);
}
