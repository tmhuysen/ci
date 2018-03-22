#include <RDMdoci.hpp>
#include <NormalGenerator.hpp>
#include "DOCI.hpp"

/**
 * Helper function for the constructors
 */
void doci::DOCI::construct() {
	// Set the number of spatial orbitals and electrons.
	size_t K_ = this->basis->getK(); //spatial
	size_t nelec_ = this->basis->getNelec();
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
				double one_int = this->basis->getOne_int(j, j);
				this->hamiltonian->add(2 * one_int, i, i); //Twice : alpha and beta.
				//There are also two in-place double excitations of abba en baab combination.
				double two_int = this->basis->getTwo_int(j, j, j, j);
				this->hamiltonian->add(two_int, i, i);

				for(size_t l = 0; l < j; l++){//creation l=j is covered in the first loop and since we can't annihilate twice this combination would be redundant.
					if (bmqc::creation(bf_base,l)){ //we can never excite a single electron to a new site we have to do it in pairs.
						size_t address = this->ad_mat.fetchAddress(bf_base);
						//integrals parameters are entered in chemical notation!
						// This means that first 2 parameters are for the first electrons and subsequent ones are for the second
						//Multiply by 2 getting rid of 1/2 two electron term because we have 2 equal combinations:
						//abba and baab. We do not correct for the truncated SO iteration because we only calculate the lower triagonal.
						//and then copy accordingly
						double mix_spin_two_int = this->basis->getTwo_int(j, l, j, l);
						this->hamiltonian->add(mix_spin_two_int, i, address);
						this->hamiltonian->add(mix_spin_two_int, address, i);

						bf_base.flip(l);//flip back (so we don't need to copy the set)
					}else{ //if we can't create we can annihilated (but only on the diagonal, no transform required

						// Integral parameters are entered in chemical notation!
						double same_spin_two_int = this->basis->getTwo_int(j, j, l, l); //=mixed_spin_two_int (exciting a beta and an alpha in-place)
						double same_spin_two_int_negative = -this->basis->getTwo_int(j, l, l, j); //mixed_spin does not have this because it would result in 0 term (integral of alpha-beta)

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
	this->lowestEigenState = this->hamiltonian->getGroundstates().at(0);



}

doci::DOCI::DOCI(
        doci::CI_basis_mulliken *ciBasis,
         double constraint,
        std::vector<size_t> set_of_AO)
        : CI(ciBasis, constraint, set_of_AO){

    // evaluate all mulliken operator values for the AO set.
    ciBasis->calculateMullikenMatrix(set_of_AO);
    double threshold = 1e-6; // set a threshold
    construct(); // define the CI dims

    double error = this->npairs*2; //maximum value
    std::cout<<error<<"start";
    double langrange_multiplier = 0; //no constraint
    double max = 1;
    double min = -1;
    double std_dev = 0.3;
    double transform_mag = 0.9999;
    double interval = 0.00002;
    size_t iterations = 0;
    NormalGenerator langrange_generator(langrange_multiplier,std_dev,min,max);

    while(std::abs(error) > threshold){
        std::cout<<std::endl<<" multi : "<<langrange_multiplier;
        ciBasis->set_lagrange_multiplier(langrange_multiplier);
        this->basis=ciBasis;
        this->hamiltonian = Hamiltonian::make_hamiltonian(nbf);
        calculateCI(0,this->nbf);
        this->hamiltonian->solve();
        this->lowestEigenState = this->hamiltonian->getGroundstates().at(0);
        delete(hamiltonian);
        rdm::RDMdoci rdm_doci = rdm::RDMdoci(lowestEigenState.getEvec(),this->K,this->npairs);
        double population = ciBasis->mullikenPopulationCI(&rdm_doci);
        //std::cout<<std::endl<<" pop: "<<population;
        double new_error = constraint - population;
        //std::cout<<std::endl<<"error: "<<new_error<< " v.s " <<error;
        if(std::abs(error)>std::abs(new_error)){
            std::cout<<std::endl<<" winner : "<<"error: "<<new_error<< " v.s " <<error<< " with lambda : "<<langrange_multiplier;
            langrange_generator.setMean(langrange_multiplier);
            langrange_generator.transform(transform_mag);
            error = new_error;

        }
        langrange_multiplier = langrange_generator.generate();
        //std::cout<<std::endl<<" multiplier: "<<langrange_multiplier;
        iterations++;
        if(iterations > 1000000){
            throw std::overflow_error("to many iterations");

        }
    }





}
