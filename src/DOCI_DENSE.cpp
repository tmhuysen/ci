//
// Created by Wulfix on 14/12/2017.
//

#include "DOCI_DENSE.h"

DOCI_DENSE::DOCI_DENSE(CI_basis calculator) : DOCI(calculator) {
    hamiltonian = Eigen::MatrixXd::Zero(nbf,nbf);
    calculateDoci(0,1);
    symmatu(hamiltonian);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonian);
    eigenvalues = solver.eigenvalues().real().cast<double>();
    eigenvectors = solver.eigenvectors().real().cast<double>();
    //We extract only the groundstate
    for (int i = 0; i<eigenvalues.size(); i++) {
        groundStates(State {eigenvalues[i], eigenvectors.col(i)});
    }
}

void DOCI_DENSE::addToHamiltonian(double value, unsigned long index1, unsigned long index2) {
    hamiltonian(index1,index2) += value;

}

Eigen::MatrixXd DOCI_DENSE::getHam() {
    return hamiltonian;
}

void DOCI_DENSE::print(){
    std::cout<<hamiltonian;
}
