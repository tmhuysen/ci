#ifndef DOCI_STATE_HPP
#define DOCI_STATE_HPP


#include <Eigen/Dense>


class State {  // FIXME: this class should be tested
private:
    double eval;  // The energy of the solution, i.e. the eigenvalue
    Eigen::VectorXd evec;  // The coefficients of the solution with respect to the given basis
                           // The eigenvector corresponding to the eigenvalue

public:

    // Constructors
    /** Constructor based on a given eigenvalue and corresponding eigenvector
     */
    State(double eval, Eigen::VectorXd evec);


    // Operator overloading
    bool operator<(const State& rhs);
    bool operator==(const State& rhs);

};

#endif // DOCI_STATE_HPP
