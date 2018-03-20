#ifndef CI_DENSESOLVER_HPP
#define CI_DENSESOLVER_HPP



#include <Eigen/Dense>

#include "BaseMatrixSolver.hpp"



namespace ci {
namespace solver {


class DenseSolver : public ci::solver::BaseMatrixSolver {
private:
    Eigen::MatrixXd hamiltonian;


public:
    /*
     *  CONSTRUCTORS
     */
    /**
     *   Constructor based on the dimension @param dim of the Hamiltonian matrix representation.
     */
    explicit DenseSolver(size_t dim);


    /*
     *  DESTRUCTOR
     */
    ~DenseSolver() override = default;


    /*
     *  PUBLIC OVERRIDDEN FUNCTIONS
     */
    /**
     *  Solve the full dense eigenvalue problem of the Hamiltonian matrix.
     */
    void solve() override;

    /**
     *  Add @param value to the matrix representation of the Hamiltonian at (@param index1, @param index2).
     */
    void addToMatrix(double value, size_t index1, size_t index2) override;
};



}  // namespace solver
}  // namespace ci



#endif  // CI_DENSESOLVER_HPP
