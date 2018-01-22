#ifndef DOCI_STATE_HPP
#define DOCI_STATE_HPP


#include <Eigen/Dense>


class State {  // FIXME: this class should be tested
public:
    // FIXME: use private members and getters
    double eval;  // The energy of the solution, i.e. the eigenvalue
    Eigen::VectorXd evec;  // The coefficients of the solution with respect to the given basis, i.e. the eigenvector corresponding to the eigenvalue

    // Constructors
    /** Constructor based on a given eigenvalue and corresponding eigenvector
     */
    State(double eval, Eigen::VectorXd evec);


    // Operator overloading
    bool operator<(const State& rhs);
    bool operator==(const State& rhs);

};

#endif // DOCI_STATE_HPP
