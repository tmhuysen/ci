#include "CI_Abstract_Class.hpp"

const doci::State &doci::CI::getLowestEigenState() const {
	return lowestEigenState;
}


/** Constructors
 */
doci::CI::CI(doci::CI_basis *ciBasis){
	basis = ciBasis;
	//construct();
	//this->hamiltonian = Hamiltonian::make_hamiltonian(nbf);
	//calculateCI(0,this->nbf);


}

doci::CI::CI() {}

doci::Hamiltonian *doci::CI::getHamiltonian() const {
	return hamiltonian;
}
