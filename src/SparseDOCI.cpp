#include <MatOp/SparseSymMatProd.h>
#include "SparseDOCI.hpp"
#include "utility.hpp"


doci::SparseDOCI::SparseDOCI(doci::CI_basis ciBasis) : doci::DOCI(ciBasis) {

    // Construct the Hamiltonian matrix
    this->hamiltonian = Eigen::SparseMatrix<double>(this->nbf,this->nbf);
    this->hamiltonian.setZero();
    calculateDoci(0,1);

    // Diagonalize the Hamiltonian matrix

    // Construct matrix operation object using the wrapper class SparseGenMatProd
    Spectra::SparseSymMatProd<double> op(hamiltonian);

    // Construct eigen solver object, requesting the largest x eigenvalues (in magnitude)
    Spectra::SymEigsSolver< double, Spectra::LARGEST_MAGN, Spectra::SparseSymMatProd<double> > eigs(&op,1,6);

    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute();
    // Retrieve results
    if(eigs.info() == Spectra::SUCCESSFUL){
        this->eigenvalues = eigs.eigenvalues();
        this->eigenvectors = eigs.eigenvectors();
    }


    // Extract only the ground state
    for (size_t i = 0; i < this->eigenvalues.size(); i++) {
        groundStates(doci::State(eigenvalues[i], eigenvectors.col(i)));
    }
}


void doci::SparseDOCI::addToHamiltonian(double value, unsigned long index1, unsigned long index2) {
    this->hamiltonian.coeffRef(index1,index2) += value;

}


Eigen::SparseMatrix<double> doci::SparseDOCI::getHamiltonian() {
    return this->hamiltonian;
}


void doci::SparseDOCI::print() {
    std::cout << hamiltonian;
}
