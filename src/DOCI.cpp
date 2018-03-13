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
    calculateDiagonal();
    do{
        calculateOffDiagonal();
    }while(!hamiltonian->solve());


}

doci::DOCI::DOCI( doci::CI_basis *ciBasis) : CI(ciBasis) {
    construct();
    this->hamiltonian = Hamiltonian::make_hamiltonian(nbf);
    calculateCI(0,this->nbf);
    this->lowestEigenState = this->hamiltonian->getGroundstates().at(0);



}

void doci::DOCI::calculateOffDiagonal() {
    boost::dynamic_bitset<> bf_base = this->ad_mat.generateBitVector_bitset(0); //first basis function
    for (size_t i = 0; i < this->nbf; i++) {
        if(i>0){
            bmqc::next_bitset_permutation(bf_base);
        }
        for (size_t j = 0; j < this->K; j++) {
            if (bmqc::annihilation(bf_base,j)){
                for(size_t l = 0; l < j; l++){
                    if (bmqc::creation(bf_base,l)) { //we can never excite a single electron to a new site we have to do it in pairs.
                        size_t address = this->ad_mat.fetchAddress(bf_base);
                        //integrals parameters are entered in chemical notation!
                        // This means that first 2 parameters are for the first electrons and subsequent ones are for the second
                        //Multiply by 2 getting rid of 1/2 two electron term because we have 2 equal combinations:
                        //abba and baab.
                        //and then copy accordingly
                        double mix_spin_two_int = this->basis->getTwo_int(j, l, j, l);
                        this->hamiltonian->add(mix_spin_two_int, i, address);
                        this->hamiltonian->add(mix_spin_two_int, address, i);

                        bf_base.flip(l);//flip back (so we don't need to copy the set)
                    }
                }
                bf_base.flip(j);//flip back (so we don't need to copy the set)
            }


        }
    }

}

void doci::DOCI::calculateDiagonal() {
    boost::dynamic_bitset<> bf_base = this->ad_mat.generateBitVector_bitset(0); //first basis function
    for (size_t i = 0; i < this->nbf; i++) {
        if(i>0){
            bmqc::next_bitset_permutation(bf_base);
        }
        double hamiltonian_value = 0;
        for (size_t j = 0; j < this->K; j++) {
            if (bf_base.test(j)){
                hamiltonian_value += 2*this->basis->getOne_int(j, j);
                hamiltonian_value += this->basis->getTwo_int(j, j, j, j);
                for(size_t l = 0; l < j; l++){
                    if(bf_base.test(l)){
                        hamiltonian_value += 4*this->basis->getTwo_int(j, j, l, l); //=mixed_spin_two_int (exciting a beta and an alpha in-place)
                        hamiltonian_value -= 2*this->basis->getTwo_int(j, l, l, j); //mixed_spin does not have this because it would result in 0 term (integral of alpha-beta)

                    }
                }
            }


        }
        this->hamiltonian->add(hamiltonian_value,i,i);
    }

}
