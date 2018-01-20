#include "DOCI_dense.hpp"

DOCI_dense::DOCI_dense(CI_basis calculator) : DOCI_Class(calculator) {
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

void DOCI_dense::addToHamiltonian(double value, unsigned long index1, unsigned long index2) {
    hamiltonian(index1,index2) += value;

}

Eigen::MatrixXd DOCI_dense::getHam() {
    return hamiltonian;
}

void DOCI_dense::print(){
    std::cout<<hamiltonian;
}
