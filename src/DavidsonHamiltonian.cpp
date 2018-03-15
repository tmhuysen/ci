#include "DavidsonHamiltonian.hpp"

doci::DavidsonHamiltonian::DavidsonHamiltonian(size_t nbf) {
    this->nbf = nbf;
    diagonal_matrix = Eigen::VectorXd::Zero(nbf);

}

void doci::DavidsonHamiltonian::add(double value, size_t index1, size_t index2) {
    if(index1==index2){
        diagonal_matrix(index1) += value;
        if(index1 == nbf-1){
            davidson_solver = new numopt::DavidsonSolver(diagonal_matrix,1,1.0e-06);
        }
    }else{
        davidson_solver->enterMatrixElement(value,index1,index2);
    }

}

bool doci::DavidsonHamiltonian::solve() {
    ++iterations;
    if(davidson_solver->solve()){

        this->eigenvalues = davidson_solver->get_eigenvalues();
        this->eigenvectors = davidson_solver->get_eigenvalues();
        for (size_t i = 0; i < this->eigenvalues.size(); i++) {
            groundStates(doci::State(eigenvalues[i], eigenvectors.col(i)));
        }
        std::cout<<" we have "<<iterations<<" itterations!";
        std::cout<<davidson_solver->get_num_guess_vectors();

        return true;

    }else{

        return false;
    }
}
