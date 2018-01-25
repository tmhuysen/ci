#ifndef DOCI_STATE_HPP
#define DOCI_STATE_HPP


#include <Eigen/Dense>

namespace doci {

class State {
private:
public:
    double getEval() const;
    const Eigen::VectorXd &getEvec() const;

private:
    double eval;  // The energy of the solution, i.e. the eigenvalue
    Eigen::VectorXd evec;  // The coefficients of the solution with respect to the given basis, i.e. the eigenvector corresponding to the eigenvalue

public:


    // Constructors
    /** Constructor based on a given eigenvalue and corresponding eigenvector
     */
    State(double eval, Eigen::VectorXd evec);


    // Operator overloading
    bool operator<(const doci::State& rhs); //only compares eval.
    bool operator>(const doci::State& rhs); //only compares eval.

    bool operator==(const doci::State& rhs); //only compares eval.

};

} // namespace doci

#endif // DOCI_STATE_HPP
