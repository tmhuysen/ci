#include "FCI.hpp"

/**
 * Helper function for the constructors
 */
void doci::FCI::construct() {
	// Set the number of spatial orbitals and electrons.
	size_t K_ = this->basis->getK();
	size_t nelec_ = this->basis->getNelec();
	if (K_ < nelec_/2) {
		throw std::overflow_error("Invalid argument: too many electrons to place into the given number of spatial orbitals");
	}
	this->nelec_b = nelec_ / 2;
	this->nelec_a = nelec_ - nelec_b;
	this->K = K_;


	// Set the number of spatial orbitals
	auto nbf_a_ = boost::math::binomial_coefficient<double>(this->K, this->nelec_a);
	if (nbf_a_ > 4294967295.0) {
		// before casting into unsigned long, we have to make sure that it can fit
		throw std::overflow_error("The number of basis functions for the separated sector is too high to be cast into size_t. - ALPHA");
	}
	this->nbf_a = nbf_a_;

	auto nbf_b_ = boost::math::binomial_coefficient<double>(this->K, this->nelec_b);
	if (nbf_b_ > 4294967295.0) {
		// before casting into unsigned long, we have to make sure that it can fit
		throw std::overflow_error("The number of basis functions for the separated sector is too high to be cast into size_t. - BETA");
	}
	this->nbf_b = nbf_b_;

	auto nbf_ = nbf_b_ * nbf_a_;
	if (nbf_ > 4294967295.0) {
		// before casting into unsigned long, we have to make sure that it can fit
		throw std::overflow_error("The number of basis functions for the sector is too high to be cast into size_t.");
	}
	this->nbf = static_cast<unsigned long>(nbf_);


	this->ad_mat_a = bmqc::AddressingScheme(this->K, this->nelec_a);
	this->ad_mat_b = bmqc::AddressingScheme(this->K, this->nelec_b);

}

/**
 * calculate hamiltonian elements.
 * @param start,end : indicates the bf you want to start with and where the iteration ends(excluded).
 */
void doci::FCI::calculateCI(size_t start, size_t end) {
	boost::dynamic_bitset<> bf_base_b = this->ad_mat_b.generateBitVector_bitset(0); //first basis function alpa
	boost::dynamic_bitset<> bf_base_a;
	for (size_t j = 0; j < this->nbf_b; j++) { //Iterate alpha
		if (j > 0) {
			bmqc::next_bitset_permutation(bf_base_b);
		}
		bf_base_a = this->ad_mat_a.generateBitVector_bitset(0); //first basis function beta
		for (size_t i = 0; i < this->nbf_a; i++) { //Iterate beta
			if (i > 0) {
				bmqc::next_bitset_permutation(bf_base_a);
			}


			for (size_t a = 0; a < this->K; a++) { //SO ITTERATION 1
				//one int A
				int sign1 = 1;
				if (bmqc::annihilation_s(bf_base_a,a,sign1)){
					for (size_t b = 0; b < this->K; b++) { //SO ITTERATION 2
						int sign2 = sign1;
						if(bmqc::creation_s(bf_base_a,b,sign2)){
							size_t address = this->ad_mat_a.fetchAddress(bf_base_a);
							double one_int = sign2 * this->basis->getOne_int(a, b);
							this->hamiltonian->add(one_int, i+j*nbf_a, address+j*nbf_a);



							//undo
							bf_base_a.flip(b);
						}
						sign2 = sign1;
						//second anni A
						if(bmqc::annihilation_s(bf_base_a,b,sign2)){
							for (size_t c = 0; c < this->K; c++) {

								int sign3 = sign2;
								if(bmqc::creation_s(bf_base_a,c,sign3)){
									for (size_t d = 0; d < this->K; d++) {

										int sign4 = sign3;
										if(bmqc::creation_s(bf_base_a,d,sign4)){

											size_t address = this->ad_mat_a.fetchAddress(bf_base_a);
											double two_int = sign4 *this->basis->getTwo_int(a,d,b,c)/2;
											this->hamiltonian->add(two_int, i+j*nbf_a, address+j*nbf_a);


											bf_base_a.flip(d);
										}
									}
									bf_base_a.flip(c);
								}


							}
							bf_base_a.flip(b);

						}
						sign2 = sign1;
						//second anni B
						if(bmqc::annihilation_s(bf_base_b,b,sign2)){
							for (size_t c = 0; c < this->K; c++) {
								int sign3 = sign2;
								if(bmqc::creation_s(bf_base_b,c,sign3)){
									for (size_t d = 0; d < this->K; d++) {
										int sign4 = sign3;
										if(bmqc::creation_s(bf_base_a,d,sign4)){
											size_t addressa = this->ad_mat_a.fetchAddress(bf_base_a);
											size_t addressb = this->ad_mat_b.fetchAddress(bf_base_b);
											double two_int = sign4*this->basis->getTwo_int(a,d,b,c); //for switcharoo
											this->hamiltonian->add(two_int, i+j*nbf_a, addressa+addressb*nbf_a);


											bf_base_a.flip(d);
										}
									}
									bf_base_b.flip(c);
								}


							}
							bf_base_b.flip(b);

						}

					}



					bf_base_a.flip(a);
				}
				sign1 = 1;
				if (bmqc::annihilation_s(bf_base_b,a,sign1)) {
					for (size_t b = 0; b < this->K; b++) { //SO ITTERATION 2
						int sign2 = sign1;
						if (bmqc::creation_s(bf_base_b, b, sign2)) {
							size_t address = this->ad_mat_b.fetchAddress(bf_base_b);
							double one_int = sign2*this->basis->getOne_int(a, b);
							this->hamiltonian->add(one_int, i + j * nbf_a, i + address * nbf_a);


							bf_base_b.flip(b);
						}
						sign2 = sign1;
						if (bmqc::annihilation_s(bf_base_b, b, sign2)) {
							for (size_t c = 0; c < this->K; c++) {
								int sign3 = sign2;
								if (bmqc::creation_s(bf_base_b, c, sign3)) {
									for (size_t d = 0; d < this->K; d++) {
										int sign4 = sign3;
										if (bmqc::creation_s(bf_base_b, d, sign4)) {
											size_t address = this->ad_mat_b.fetchAddress(bf_base_b);
											double two_int = sign4*this->basis->getTwo_int(a, d, b, c)/2;
											this->hamiltonian->add(two_int, i + j * nbf_a, i + address * nbf_a);

											bf_base_b.flip(d);
										}
									}
									bf_base_b.flip(c);
								}


							}
							bf_base_b.flip(b);

						}
					}
					bf_base_b.flip(a);
				}
			}
		}
	}
}






doci::FCI::FCI(doci::CI_basis *ciBasis) : CI(ciBasis) {
	construct();
	this->hamiltonian = Hamiltonian::make_hamiltonian(nbf);
	calculateCI(0,this->nbf);
	this->hamiltonian->solve();
	this->lowestEigenState = this->hamiltonian->getGroundstates().at(0);

}
