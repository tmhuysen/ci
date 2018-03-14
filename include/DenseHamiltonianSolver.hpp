#ifndef CI_DENSEHAMILTONIANSOLVER_HPP
#define CI_DENSEHAMILTONIANSOLVER_HPP



#include "BaseHamiltonianSolver.hpp"



namespace ci {


/**
 *
 */
class DenseHamiltonianSolver: public ci::BaseHamiltonianSolver {
private:
    Eigen::MatrixXd hamiltonian;  // the dense Hamiltonian matrix in ONV basis (i.e. in Fock space)

    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param dimension.
     */
    DenseHamiltonianSolver(size_t dimension);


    // GETTERS
    Eigen::MatrixXd get_hamiltonian() const { return this->hamiltonian; }


    // PUBLIC OVERRIDDEN VIRTUAL METHODS
    /**
     *  Solves the eigenvalue problem of the Hamiltonian.
     *
     *  For this derived class, this means diagonalizing the full dense Hamiltonian.
     */
    void solve() override;

    /**
     *  Adds @param value to the current value of the Hamiltonian matrix at position (@param index1, @param index2).
     */
    void addElement(double value, size_t index1, size_t index2) override;
};


}  // namespace ci



#endif  // CI_DENSEHAMILTONIANSOLVER_HPP
