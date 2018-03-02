#include "DenseHamiltonian.hpp"

//Constructor
doci::DenseHamiltonian::DenseHamiltonian(size_t nbf) {
	this->nbf = nbf;
	this->hamiltonian = Eigen::MatrixXd::Zero(this->nbf, this->nbf);
}
//Public overridden functions.
void doci::DenseHamiltonian::solve() {
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver (this->hamiltonian);
	this->eigenvalues = solver.eigenvalues().real().cast<double>();
	this->eigenvectors = solver.eigenvectors().real().cast<double>();

	// Extract only the ground states
	for (size_t i = 0; i < this->eigenvalues.size(); i++) {
		groundStates(doci::State(eigenvalues[i], eigenvectors.col(i)));
	}

}

void doci::DenseHamiltonian::add(double value, size_t index1, size_t index2) {
	this->hamiltonian(index1,index2) += value;

}

const Eigen::MatrixXd &doci::DenseHamiltonian::getHamiltonian() const {
	return hamiltonian;
}
