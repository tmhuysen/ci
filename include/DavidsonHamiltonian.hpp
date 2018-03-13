#ifndef CI_DAVIDSONHAMILTONIAN_HPP
#define CI_DAVIDSONHAMILTONIAN_HPP


#include "Hamiltonian.hpp"
#include "DavidsonSolver.hpp"

namespace doci {


/**
* Davidson Hamiltonian for calculations
* where the hamiltonian cannot be stored in memory
* recommended to only retrieve a limited amount of eigenvectors.
* will always retrieve lowest.
*/
class DavidsonHamiltonian: public doci::Hamiltonian {
private:
    StorageType type = StorageType::DAVIDSON;
    numopt::DavidsonSolver* davidson_solver;
    size_t nbf;
    Eigen::VectorXd diagonal_matrix;


public:
    /**
     * Constructor
     * @param nbf the dimensions of hamiltonian(number of basis functions)
     */
    DavidsonHamiltonian(size_t nbf);
    /**
     *  Solves the eigenvalue problem of the hamiltonian with the EigenSolver.
     */
    bool solve() override;
    /**
     * Adds a value as an element of the hamiltonian matrix.
     */
    void add(double value, size_t index1, size_t index2) override;

};


} // namespace doci


#endif // CI_DENSEHAMILTONIAN_HPP
