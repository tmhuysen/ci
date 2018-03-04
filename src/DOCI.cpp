#include "DOCI.hpp"

/**
 * Helper function for the constructors
 */
void doci::DOCI::construct() {
	// Set the number of spatial orbitals and electrons.
	size_t K_ = this->basis->get_K(); //spatial
	size_t nelec_ = this->basis->get_N();
	if (K_ < nelec_/2) {
		throw std::overflow_error("Invalid argument: too many electrons to place into the given number of spatial orbitals");
	}
	if(nelec_%2){
		throw std::invalid_argument("Your basis contains an odd amount of electrons and is not suitable for DOCI");
	}
	this->npairs = nelec_/2;
	this->K = K_;


	// Set the number of basis functions
	auto nbf_ = boost::math::binomial_coefficient<double>(this->K, this->npairs);
	if (nbf_ > 4294967295.0) {
		// before casting into unsigned long, we have to make sure that it can fit
		throw std::overflow_error("The number of basis functions for the sector is too high to be cast into unsigned long.");
	}
	this->nbf = static_cast<unsigned long>(nbf_);

	this->ad_mat = bmqc::AddressingScheme(this->K, this->npairs); //constructing Addressing Scheme


}

/**
 * calculate hamiltonian elements.
 * @param start,end : indicates the bf you want to start with and where the iteration ends(excluded).
 */
void doci::DOCI::calculateCI(size_t start, size_t end) {
	boost::dynamic_bitset<> bf_base = this->ad_mat.generateBitVector_bitset(start); //first basis function
	for (size_t i = 0; i < end; i++) {
		if(i>0){
			bmqc::next_bitset_permutation(bf_base);
		}

		for (size_t j = 0; j < this->K; j++) { //First iteration over SO's.
			if (bmqc::annihilation(bf_base,j)){ //single excitation
				//A single excitation in doci can only be done in place.
				//Exciting only one electron to a vacant SO, will break the double occupancy(not part of the basis).
				double one_int = this->basis->get_one_int(j, j);
				this->hamiltonian->add(2 * one_int, i, i); //Twice : alpha and beta.
				//There are also two in-place double excitations of abba en baab combination.
				double two_int = this->basis->get_two_int(j, j, j, j);
				this->hamiltonian->add(two_int, i, i);

				for(size_t l = 0; l < j; l++){//creation l=j is covered in the first loop and since we can't annihilate twice this combination would be redundant.
					if (bmqc::creation(bf_base,l)){ //we can never excite a single electron to a new site we have to do it in pairs.
						size_t address = this->ad_mat.fetchAddress(bf_base);
						//integrals parameters are entered in chemical notation!
						// This means that first 2 parameters are for the first electrons and subsequent ones are for the second
						//Multiply by 2 getting rid of 1/2 two electron term because we have 2 equal combinations:
						//abba and baab. We do not correct for the truncated SO iteration because we only calculate the lower triagonal.
						//and then copy accordingly
						double mix_spin_two_int = this->basis->get_two_int(j, l, j, l);
						this->hamiltonian->add(mix_spin_two_int, i, address);
						this->hamiltonian->add(mix_spin_two_int, address, i);

						bf_base.flip(l);//flip back (so we don't need to copy the set)
					}else{ //if we can't create we can annihilated (but only on the diagonal, no transform required

						// Integral parameters are entered in chemical notation!
						double same_spin_two_int = this->basis->get_two_int(j, j, l, l); //=mixed_spin_two_int (exciting a beta and an alpha in-place)
						double same_spin_two_int_negative = -this->basis->get_two_int(j, l, l, j); //mixed_spin does not have this because it would result in 0 term (integral of alpha-beta)

						//We don't iterate over all the SO's the second time so multiply by 2 getting rid of 1/2 two electron term.
						//multiply by 2 again because alpha,alpha is the same as beta,beta combinations.
						//same_spin (positive) = mixed, so multiply that by 2 again.
						this->hamiltonian->add((4 * same_spin_two_int + 2 * same_spin_two_int_negative), i, i);

					}

				}
				bf_base.flip(j);//flip back (so we don't need to copy the set)
			}


		}
	}
}

doci::DOCI::DOCI( doci::CI_basis *ciBasis) : CI(ciBasis) {
	construct();
	this->hamiltonian = Hamiltonian::make_hamiltonian(nbf);
	calculateCI(0,this->nbf);
	this->hamiltonian->solve();
	this->lowest_eigenstate = this->hamiltonian->get_groundstates().at(0);



}
