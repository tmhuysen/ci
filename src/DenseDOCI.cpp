#include "DenseDOCI.hpp"

DenseDOCI::DenseDOCI(CI_basis calculator) : DOCI_Class(calculator) {
    hamiltonian = Eigen::MatrixXd::Zero(nbf,nbf);
    calculateDoci(0,1);
    symmatu_reverse(hamiltonian);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonian);
    eigenvalues = solver.eigenvalues().real().cast<double>();
    eigenvectors = solver.eigenvectors().real().cast<double>();
    //We extract only the groundstate
    for (int i = 0; i<eigenvalues.size(); i++) {
        groundStates(State {eigenvalues[i], eigenvectors.col(i)});
    }
}

void DenseDOCI::addToHamiltonian(double value, unsigned long index1, unsigned long index2) {
    hamiltonian(index1,index2) += value;

}

Eigen::MatrixXd DenseDOCI::getHam() {
    return hamiltonian;
}

void DenseDOCI::print(){
    std::cout<<hamiltonian;
}
