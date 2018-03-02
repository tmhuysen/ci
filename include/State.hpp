#ifndef CI_STATE_HPP
#define CI_STATE_HPP


#include <Eigen/Dense>

namespace doci {

class State {
private:
    double eval;  // The energy of the solution, i.e. the eigenvalue
    Eigen::VectorXd evec;  // The coefficients of the solution with respect to the given basis, i.e. the eigenvector corresponding to the eigenvalue

public:
    // Constructors
    /**
     * Default constructor
     */
    State();
    /** Constructor based on a given eigenvalue and corresponding eigenvector
     */
    State(double eval, Eigen::VectorXd evec);

    // Getters
    double getEval() const;
    const Eigen::VectorXd &getEvec() const;

    // Operator overloading
    bool operator<(const doci::State& rhs); //only compares eval.
    bool operator>(const doci::State& rhs); //only compares eval.
    bool operator==(const doci::State& rhs); //only compares eval.

};

} // namespace doci


#endif // CI_STATE_HPP
