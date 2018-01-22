#include "DOCI_Class.hpp"


/** Constructor based on a given CI_basis
 */
doci::DOCI::DOCI(doci::CI_basis ciBasis) {

    // Set the number of spatial orbitals and electron pairs
    size_t K_ = ciBasis.nbf;
    size_t npairs_ = ciBasis.nelec / 2;

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


    this->ad_mat = AddressingMatrix(this->K, this->npairs);
    this->groundstates = { doci::State (std::numeric_limits<double>::max(), Eigen::VectorXd()) };
    this->basis = ciBasis;
}




void doci::DOCI::calculateDoci(double start, double end) {
    boost::dynamic_bitset<> basic_bit = this->ad_mat.generateBinaryVector(start * this->nbf);
    for (size_t i = 0; i < this->nbf * end; i++) {
        for (size_t j = 0; j < this->K; j++) {
            if (basic_bit[j]){
                double one_int = this->basis.one_ints(j,j);
                addToHamiltonian(2 * one_int, i, i);
            }
            for(size_t l = 0; l < j+1; l++){
                if( j!=l) {
                    boost::dynamic_bitset<> two_target_dia = basic_bit;
                    if (annihilation(two_target_dia, j) && annihilation(two_target_dia, l)){
                        // Integral parameters are entered in chemical notation!
                        // This means that first 2 parameters are for the first electrons and subsequent ones are for the second
                        double same_spin_two_int = this->basis.two_ints(j,j,l,l);
                        double mix_spin_two_int = same_spin_two_int; //just illustrative
                        double same_spin_two_int_negative = -this->basis.two_ints(j,l,l,j);
                        addToHamiltonian((4 * same_spin_two_int + 2 * same_spin_two_int_negative), i, i);
                    }
                }
                boost::dynamic_bitset<> two_target = basic_bit;
                if (annihilation(two_target, j) && creation(two_target, l)) {
                    size_t address = this->ad_mat.fetchAddress(two_target);
                    //integrals parameters are entered in chemical notation!
                    double mix_spin_two_int = this->basis.two_ints(j,l,j,l);
                    addToHamiltonian(mix_spin_two_int, i, address);
                }
            }
        }
        basic_bit = next_bitset_permutation(basic_bit);
    }
}


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


const std::vector<doci::State>& doci::DOCI::getGroundstates() const {
    return this->groundstates;
};





