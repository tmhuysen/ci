#include "CI_Abstract_Class.hpp"




/** Default Constructor
 */
doci::CI::CI() {}

doci::CI::CI(doci::CI_basis *ciBasis){
    basis = ciBasis;
}

/** Destructor
 * clears the hamiltonian pointer.
 */
doci::CI::~CI() {
    delete hamiltonian;
}

/** Getters
 */
const doci::State &doci::CI::getLowestEigenState() const {
    return lowestEigenState;
}


doci::Hamiltonian *doci::CI::getHamiltonian() const {
    return hamiltonian;
}

