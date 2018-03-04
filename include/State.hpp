#ifndef CI_STATE_HPP
#define CI_STATE_HPP


#include <Eigen/Dense>

namespace doci {

class State {
private:
    double eigenvalue;  // The energy of the solution, i.e. the eigenvalue
    Eigen::VectorXd eigenvector;  // The coefficients of the solution with respect to the given basis, i.e. the eigenvector corresponding to the eigenvalue

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
    double get_eigenvalue() const;
    const Eigen::VectorXd &get_eigenvector() const;

    // Operator overloading
    bool operator<(const doci::State& rhs); //only compares eigenvalue.
    bool operator>(const doci::State& rhs); //only compares eigenvalue.
    bool operator==(const doci::State& rhs); //only compares eigenvalue.

};

} // namespace doci


#endif // CI_STATE_HPP
