#ifndef CI_CI_BASIS_MULLIKEN_HPP
#define CI_CI_BASIS_MULLIKEN_HPP

#include "CI_basis.hpp"
namespace doci {

class CI_basis_mulliken : public CI_basis {
private:
    double lagrange_multiplier = 0;

    Eigen::MatrixXd C; //Canonical mat hf
    Eigen::MatrixXd S; //Overlap matrix
    Eigen::MatrixXd S_inverse; //Overlap matrix
    Eigen::MatrixXd mulliken_matrix; //Matrix contains evaluation of the mulliken operator (one electron operator)
private:
    double evaluateMullikenOperator(size_t molecular_orbital1,size_t molecular_orbital2,size_t atomic_orbital);

public:
    /** Constructor based on a given RHF instance
     */
    CI_basis_mulliken(hf::rhf::RHF& rhf);
    /** Calculates the mulliken matrix for a set of AO's
     */
    void calculateMullikenMatrix(std::vector<size_t> set_of_AO);

    /** getOne_int now returns the one electron value + langrange multiplied corresponding value of the mulliken matrix.
     */
    double getOne_int(size_t index1, size_t index2) const override;
    /** Setters
     */
    void set_lagrange_multiplier(double lagrange_multiplier);
};
}

#endif //CI_CI_BASIS_MULLIKEN_HPP
