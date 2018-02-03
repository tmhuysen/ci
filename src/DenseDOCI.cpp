#include "DenseDOCI.hpp"

#include "utility.hpp"

/** Constructor based on a given CI_basis
 * Applies the base DOCI_class constructor and calls the DOCI calculation and
 * solves the eigenvalues of the hamiltonian with the EigenSolver.
 */
doci::DenseDOCI::DenseDOCI(doci::CI_basis ciBasis) : doci::DOCI(ciBasis) {

    // Construct the Hamiltonian matrix
    this->hamiltonian = Eigen::MatrixXd::Zero(this->nbf, this->nbf);
    calculateDoci(0,1);

    // Diagonalize the Hamiltonian matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver (this->hamiltonian);
    this->eigenvalues = solver.eigenvalues().real().cast<double>();
    this->eigenvectors = solver.eigenvectors().real().cast<double>();

    // Extract only the ground state
    for (size_t i = 0; i < this->eigenvalues.size(); i++) {
        groundStates(doci::State(eigenvalues[i], eigenvectors.col(i)));
    }
}

//Protected overridden function.
void doci::DenseDOCI::addToHamiltonian(double value, unsigned long index1, unsigned long index2) {
    this->hamiltonian(index1,index2) += value;

}
//Public overridden function.
void doci::DenseDOCI::print() {
    std::cout << hamiltonian;
}

/**
 * Getters
 */

Eigen::MatrixXd doci::DenseDOCI::getHamiltonian() {
    return this->hamiltonian;
}