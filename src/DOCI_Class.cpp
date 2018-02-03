#include "DOCI_Class.hpp"


/** Constructor based on a given CI_basis
 */
doci::DOCI::DOCI(CI_basis &ciBasis) : basis(ciBasis) {

    // Set the number of spatial orbitals and electron pairs
    size_t K_ = basis.getK();
    size_t npairs_ = basis.getNelec() / 2;

    if (K_ < npairs_) {
        throw std::overflow_error("Invalid argument: too many electrons to place into the given number of spatial orbitals");
    }

    this->K = K_;
    this->npairs = npairs_;

    // Set the number of spatial orbitals
    auto nbf_ = boost::math::binomial_coefficient<double>(this->K, this->npairs);
    if (nbf_ > 4294967295.0) {
        // before casting into unsigned long, we have to make sure that it can fit
        throw std::overflow_error("The number of basis functions for the sector is too high to be cast into unsigned long.");
    }
    this->nbf = static_cast<unsigned long>(nbf_);


    this->ad_mat = bmqc::AddressingScheme(this->K, this->npairs); //constructing Addressing Scheme
    //initializing max State as starting groundstate.
    this->groundstates = { doci::State (std::numeric_limits<double>::max(), Eigen::VectorXd()) };
}


/**
* calculate hamiltonian elements (the lower triagonal+copy to upper triagonal).
* @param start,end : (for parallellization?) calculates only the fraction between start and end
* example: start:O.5 to end:0.75. currently excludes based on nbf (iterates over fraction of bf)
*/

void doci::DOCI::calculateDoci(double start, double end) {
    boost::dynamic_bitset<> bf_base = this->ad_mat.generateBitVector_bitset(start * this->nbf); //first basis function
    for (size_t i = 0; i < this->nbf * end; i++) {
        if(i>0){
            bmqc::next_bitset_permutation(bf_base);
        }

        for (size_t j = 0; j < this->K; j++) { //First iteration over SO's.
            if (bmqc::annihilation(bf_base,j)){ //single excitation
                //A single excitation in doci can only be done in place.
                //Exciting only one electron to a vacant SO, will break the double occupancy(not part of the basis).
                double one_int = this->basis.getOne_int(j, j);
                addToHamiltonian(2 * one_int, i, i); //Twice : alpha and beta.
                //There are also two in-place double excitations of abba en baab combination.
                double two_int = this->basis.getTwo_int(j, j, j, j);
                addToHamiltonian(two_int, i, i);

                for(size_t l = 0; l < j; l++){//creation l=j is covered in the first loop and since we can't annihilate twice this combination would be redundant.
                    if (bmqc::creation(bf_base,l)){ //we can never excite a single electron to a new site we have to do it in pairs.
                        size_t address = this->ad_mat.fetchAddress(bf_base);
                        //integrals parameters are entered in chemical notation!
                        //Multiply by 2 getting rid of 1/2 two electron term because we have 2 equal combinations:
                        //abba and baab. We do not correct for the truncated SO iteration because we only calculate the lower triagonal.
                        //and then copy accordingly
                        double mix_spin_two_int = this->basis.getTwo_int(j, l, j, l);
                        addToHamiltonian(mix_spin_two_int, i, address);
                        addToHamiltonian(mix_spin_two_int, address, i);

                        bf_base.flip(l);//flip back (so we don't need to copy the set)
                    }else{ //if we can't create we can annihilated (but only on the diagonal, no transform required

                        // Integral parameters are entered in chemical notation!
                        // This means that first 2 parameters are for the first electrons and subsequent ones are for the second
                        double same_spin_two_int = this->basis.getTwo_int(j, j, l, l); //=mixed_spin_two_int (exciting a beta and an alpha in-place)
                        double same_spin_two_int_negative = -this->basis.getTwo_int(j, l, l, j); //mixed_spin does not have this because it would result in 0 term (integral of alpha-beta)

                        //We don't iterate over all the SO's the second time so multiply by 2 getting rid of 1/2 two electron term.
                        //multiply by 2 again because alpha,alpha is the same as beta,beta combinations.
                        //same_spin (positive) = mixed, so multiply that by 2 again.
                        addToHamiltonian((4 * same_spin_two_int + 2 * same_spin_two_int_negative), i, i);

                    }

                }
                bf_base.flip(j);//flip back (so we don't need to copy the set)
            }


        }
    }
}

/**
 * Adds a @param state to the groundstates vector of our DOCI,
 * if this state's eigenvalue is equal to the eigenvalue of current groundstates, then it is added.
 * if this state's eigenvalue is lower than it replaced the current groundstates.
 */

void doci::DOCI::groundStates(doci::State state) {
    if (state == this->groundstates.at(0)) {
        this->groundstates.push_back(state);
    }

    else {
        if (state < this->groundstates.at(0)) {
            this->groundstates = std::vector<State> {state};
        }
    }
}

/**
 * Getters
 */
const std::vector<doci::State>& doci::DOCI::getGroundstates() const {
    return this->groundstates;
};





