//
// Created by Wulfix on 04/12/2017.
//


#include "DOCI.h"

DOCI::DOCI(StaticWrapper& calculator) {
    this->sites = calculator.getN_bf();
    this->electrons =  calculator.getN_electrons()/2;

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
    integralCalculator = &calculator;


}




void DOCI::calculateDoci(double start, double end) {
    boost::dynamic_bitset<> basic_bit = ad_mat.generateBinaryVector(start * nbf);
    for (unsigned long i = 0; i < nbf * end; i++) {
        for (unsigned long j = 0; j < sites; j++) {
            if (basic_bit[j]){
                double overlap = integralCalculator->calculateOverlap(j,j);
                addToHamiltonian(2*overlap,i,i);
            }
            for(unsigned long l = j; l < sites; l++){
                if(j!=l){
                    boost::dynamic_bitset<> two_target_dia = basic_bit;
                    if (annihilation(two_target_dia, j) && annihilation(two_target_dia, l)){
                        //Overlap parameters are entered in chemical notation!
                        //This means that first 2 parameters are for the first electrons and subsequent ones are for the second
                        double overlap1 = integralCalculator->calculateOverlap(j,j,l,l);
                        double overlap4 = integralCalculator->calculateOverlap(j,l,l,j);
                        overlap4 *= -1;
                        addToHamiltonian((4*overlap1+2*overlap4),i,i);
                    }
                }
                boost::dynamic_bitset<> two_target = basic_bit;
                if (annihilation(two_target, j) && creation(two_target, l)){
                    unsigned long address = ad_mat.fetchAddress(two_target);
                    //Overlap parameters are entered in chemical notation!
                    double overlap = integralCalculator->calculateOverlap(j,l,j,l);
                    addToHamiltonian(overlap,i,address);
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





bool compareState(const State &o1, const State &o2) {
    return o1.eigenValue < o2.eigenValue;
}


bool areSame(const State &o1, const State &o2) {
    double precision = 10000000; //

    double ELIPSON = (o1.eigenValue > o2.eigenValue) ?  o2.eigenValue/precision : o1.eigenValue/precision   ;

    return fabs(o1.eigenValue - o2.eigenValue) < fabs(ELIPSON);
}
