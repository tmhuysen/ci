#ifndef CI_STATE_HPP
#define CI_STATE_HPP



#include <Eigen/Dense>



namespace ci {


class State {
private:
    const double eigenvalue;  // the energy of the solution, i.e. the eigenvalue
    const Eigen::VectorXd eigenvector;  // the coefficients of the solution with respect to the given basis, i.e. the eigenvector corresponding to the eigenvalue

public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given eigenvalue and corresponding eigenvector
     */
    State(double eigenvalue, const Eigen::VectorXd& eigenvector);


    // GETTERS
    double get_eigenvalue() const { return this->eigenvalue; }
    Eigen::VectorXd get_eigenvector() const { return this->eigenvector; }


    // PUBLIC METHODS
    /**
     *  @return if this is degenerate to @param other, i.e. the eigenvalues are the same within the given @param tolerance.
     *
     *  @param tolerance: defaults to 1.0e-06
     */
    bool isDegenerate(const ci::State& other, double tolerance = 1.0e-06) const;

    /**
     *  @return if this has an eigenvalue lower than @param other, within a given @param tolerance.
     *
     *  @param tolerance: defaults to 1.0e-06
     */
    bool hasLowerEnergy(const ci::State& other, double tolerance = 1.0e-06) const;
};


} // namespace ci


#endif  // CI_STATE_HPP
