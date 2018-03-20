#ifndef CI_SPARSESOLVER_HPP
#define CI_SPARSESOLVER_HPP



#include <Eigen/Sparse>

#include "BaseMatrixSolver.hpp"



namespace ci {
namespace solver {


class SparseSolver : public ci::solver::BaseMatrixSolver {
private:
    Eigen::SparseMatrix<double> hamiltonian;

public:
    /*
     *  DESTRUCTOR
     */
    ~SparseSolver() override = default;


    /*
     *  PUBLIC OVERRIDDEN FUNCTIONS
     */
    /**
     *  Solve the sparse eigenvalue problem of the Hamiltonian matrix.
     */
    void solve() override;

    /**
     *  Add @param value to the matrix representation of the Hamiltonian at (@param index1, @param index2).
     */
    void addToMatrix(double value, size_t index1, size_t index2) override;
};



}  // namespace solver
}  // namespace ci



#endif  // CI_SPARSESOLVER_HPP
