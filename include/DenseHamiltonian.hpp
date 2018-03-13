#ifndef CI_DENSEHAMILTONIAN_HPP
#define CI_DENSEHAMILTONIAN_HPP


#include "Hamiltonian.hpp"

namespace doci {


/**
 * Dense Hamiltonian for calculations
 * where the hamiltonian is stored in a dense matrix from the eigen lib
 * High memory requirements but fast diagonalization.
 */
class DenseHamiltonian: public doci::Hamiltonian {
private:
	Eigen::MatrixXd hamiltonian;

public:
	/**
	 * Constructor
	 * @param nbf the dimensions of hamiltonian(number of basis functions)
	 */
	DenseHamiltonian(size_t nbf);
	/**
	 *  Solves the eigenvalue problem of the hamiltonian with the EigenSolver.
	 */
    bool solve() override;
	/**
	 * Adds a value as an element of the hamiltonian matrix.
	 */
	void add(double value, size_t index1, size_t index2) override;

	/**
	 * Getters
	 */
	const Eigen::MatrixXd& getHamiltonian() const;
};

} // namespace doci


#endif // CI_DENSEHAMILTONIAN_HPP
