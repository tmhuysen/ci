#include "CI_basis_mulliken.hpp"
#include <functional>

doci::CI_basis_mulliken::CI_basis_mulliken(hf::rhf::RHF &rhf) : CI_basis(rhf) {
    this->S = rhf.basis.S;
    this->S_inverse = this->S.inverse();
    this->C = rhf.C_canonical;
    this->mulliken_matrix = Eigen::MatrixXd::Zero(this->K,this->K);



}

double doci::CI_basis_mulliken::evaluateMullikenOperator(size_t molecular_orbital1, size_t molecular_orbital2,
                                                         size_t atomic_orbital) {
    //FIXME vector products should work?
    double mulliken_AO_vector_sum= 0;
    for(int i = 0; i<this->K;i++){
        for(int j = 0;j<this->K;j++){
            mulliken_AO_vector_sum += this->S_inverse(atomic_orbital,i)*this->C(j,molecular_orbital2)*this->S(i,j);
        }
    }
    double mulliken_evaluation = 0;
    for(int i = 0; i<this->K;i++){
        mulliken_evaluation += this->C.transpose()(molecular_orbital1,i)*this->S(i,atomic_orbital)*mulliken_AO_vector_sum;

    }
    return mulliken_evaluation;
}

void doci::CI_basis_mulliken::calculateMullikenMatrix(std::vector<size_t> set_of_AO) {
    
    this->mulliken_matrix = Eigen::MatrixXd::Zero(this->K,this->K);
    for(size_t ao : set_of_AO){
        for(size_t i = 0; i<this->K;i++) {
            double mulliken_evaluation_diagonal = evaluateMullikenOperator(i,i,ao);
            mulliken_matrix(i,i) += mulliken_evaluation_diagonal;
            for (size_t j = 0; j < i; j++) {
                double mulliken_evaluation = evaluateMullikenOperator(i,j,ao)/2 + evaluateMullikenOperator(j,i,ao)/2; // take the hermitian evaluation
                mulliken_matrix(i,j) += mulliken_evaluation;
                mulliken_matrix(j,i) += mulliken_evaluation;

            }
        }

    }




}

double doci::CI_basis_mulliken::getOne_int(size_t index1, size_t index2) const {
    return CI_basis::getOne_int(index1, index2) + lagrange_multiplier*mulliken_matrix(index1,index2);
}

void doci::CI_basis_mulliken::set_lagrange_multiplier(double lagrange_multiplier) {
    CI_basis_mulliken::lagrange_multiplier = lagrange_multiplier;
}

double doci::CI_basis_mulliken::mullikenPopulationCI(rdm::RDM_class *rdm) {
    rdm->compute1RDM();
    const Eigen::MatrixXd &RDMaa = rdm->getOneRDMaa();
    const Eigen::MatrixXd &RDMbb = rdm->getOneRDMbb();
    double mulliken_population = 0;
    for(size_t i = 0; i<this->K;i++) {
        mulliken_population += mulliken_matrix(i,i)*(RDMaa(i,i)+RDMbb(i,i));
    }
    return mulliken_population;



}
