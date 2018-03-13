#include <DenseHamiltonian.hpp>
#include <DavidsonHamiltonian.hpp>



/**
 * Adds a @param state to the groundstates vector of our Hamiltonian,
 * if this state's eigenvalue is equal to the eigenvalue of current groundstates, it is added.
 * if this state's eigenvalue is lower then it replaces the current groundstates.
 */

void doci::Hamiltonian::groundStates(doci::State state) {
    if (state == this->groundstates.at(0)) {
        this->groundstates.push_back(state);
    } else {
        if (state < this->groundstates.at(0)) {
            this->groundstates = std::vector<State> {state};
        }

    }
}


/*
 * Getters
 */
const std::vector<doci::State>& doci::Hamiltonian::getGroundstates() const {
    return this->groundstates;
}

/*
 * STATIC METHODS
 */

doci::Hamiltonian *doci::Hamiltonian::make_hamiltonian(size_t nbf) {
    return new doci::DenseHamiltonian(nbf);
}


doci::Hamiltonian *doci::Hamiltonian::make_dense(size_t nbf) {
    return new doci::DenseHamiltonian(nbf);
}

doci::Hamiltonian *doci::Hamiltonian::make_davidson(size_t nbf) {
    return new doci::DavidsonHamiltonian(nbf);
}

doci::Hamiltonian *doci::Hamiltonian::make_hamiltonian(size_t nbf, StorageType storageType) {
    return pFactoryMethods[(int)storageType](nbf);
}

doci::Hamiltonian* (*doci::Hamiltonian::pFactoryMethods[])(size_t nbf) = {doci::Hamiltonian::make_dense,doci::Hamiltonian::make_davidson};
