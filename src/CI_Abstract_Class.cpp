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
const doci::State &doci::CI::get_lowest_eigenstate() const {
    return lowest_eigenstate;
}


doci::Hamiltonian *doci::CI::get_hamiltonian() const {
    return hamiltonian;
}

