#include "CI_Abstract_Class.hpp"




/** Constructors
 */
doci::CI::CI() {}

doci::CI::CI(doci::CI_basis *ciBasis){
	basis = ciBasis;
}


const doci::State &doci::CI::getLowestEigenState() const {
	return lowestEigenState;
}


doci::Hamiltonian *doci::CI::getHamiltonian() const {
	return hamiltonian;
}

doci::CI::CI(doci::CI_basis_mulliken *ciBasis, double constraint, std::vector<size_t> set_of_AO) {
    basis = ciBasis;


}

