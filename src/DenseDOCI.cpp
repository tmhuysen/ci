#include "DenseDOCI.hpp"

#include "utility.hpp"


DenseDOCI::DenseDOCI(doci::CI_basis calculator) : DOCI(calculator) {

    // Construct the Hamiltonian matrix
    this->hamiltonian = Eigen::MatrixXd::Zero(this->nbf, this->nbf);
    calculateDoci(0,1);
    symmatu_reverse(this->hamiltonian);

    // Diagonalize the Hamiltonian matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver (this->hamiltonian);
    this->eigenvalues = solver.eigenvalues().real().cast<double>();
    this->eigenvectors = solver.eigenvectors().real().cast<double>();

    // Extract only the ground state
    for (size_t i = 0; i < this->eigenvalues.size(); i++) {
        groundStates(State (eigenvalues[i], eigenvectors.col(i)));
    }
}

void DenseDOCI::addToHamiltonian(double value, unsigned long index1, unsigned long index2) {
    this->hamiltonian(index1,index2) += value;

}

Eigen::MatrixXd DenseDOCI::getHam() {
    return this->hamiltonian;
}

void DenseDOCI::print() {
    std::cout << hamiltonian;
}
