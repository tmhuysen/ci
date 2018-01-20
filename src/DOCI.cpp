#include <DOCI_utility.hpp>
#include "DOCI.hpp"

DOCI::DOCI(CI_basis ciBasis) {
    this->sites = ciBasis.n_bf;
    this->electrons =  ciBasis.n_electrons/2;

    if (sites < electrons) {
        throw std::overflow_error("Invalid argument: to many electrons");
    }
    auto nbf_ = boost::math::binomial_coefficient<double>(sites, electrons);
    if (nbf_ > 4294967295.0) {
        // before casting into unsigned long, we have to make sure that it can fit
        throw std::overflow_error("The number of basis functions for the sector is too high to be cast into unsigned long.");
    }
    this->nbf = static_cast<unsigned long>(nbf_);


    ad_mat = AddressingMatrix(sites,electrons);
    groundstates = { State {std::numeric_limits<double>::max(),Eigen::VectorXd()} };
    basis = ciBasis;


}




void DOCI::calculateDoci(double start, double end) {
    boost::dynamic_bitset<> basic_bit = ad_mat.generateBinaryVector(start * nbf);
    for (size_t i = 0; i < nbf * end; i++) {
        for (size_t j = 0; j < sites; j++) {
            if (basic_bit[j]){
                double one_int = basis.one_int(j,j);
                addToHamiltonian(2*one_int,i,i);
            }
            for(size_t l = 0; l < j+1; l++){
                if(j!=l){
                    boost::dynamic_bitset<> two_target_dia = basic_bit;
                    if (annihilation(two_target_dia, j) && annihilation(two_target_dia, l)){
                        //integral parameters are entered in chemical notation!
                        //This means that first 2 parameters are for the first electrons and subsequent ones are for the second
                        double same_spin_two_int = basis.two_int(j,j,l,l);
                        double mix_spin_two_int = same_spin_two_int; //just illustrative
                        double same_spin_two_int_negative = -1*basis.two_int(j,l,l,j);
                        addToHamiltonian((4*same_spin_two_int+2*same_spin_two_int_negative),i,i);
                    }
                }
                boost::dynamic_bitset<> two_target = basic_bit;
                if (annihilation(two_target, j) && creation(two_target, l)){
                    size_t address = ad_mat.fetchAddress(two_target);
                    //integrals parameters are entered in chemical notation!
                    double mix_spin_two_int = basis.two_int(j,l,j,l);
                    addToHamiltonian(mix_spin_two_int,i,address);
                }
            }
        }
        basic_bit = next_bitset_permutation(basic_bit);
    }
}




void DOCI::groundStates(State state) {
    if(areSame(state,this->groundstates.at(0))){
        groundstates.push_back(state);
    }
    else{
        if(compareState(state,groundstates.at(0))){
            groundstates = std::vector<State> {state};
        }
    }

}

const std::vector<State> &DOCI::getGroundstates() const {
    return groundstates;
};





