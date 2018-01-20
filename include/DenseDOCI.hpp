#ifndef DOCI_DOCI_DENSE_HPP
#define DOCI_DOCI_DENSE_HPP


#include "DOCI_Class.hpp"


class DenseDOCI : public DOCI_Class {
public:
    DenseDOCI(CI_basis basis);
private:
    Eigen::MatrixXd hamiltonian;
protected:
    void addToHamiltonian(double value, size_t index1, size_t index2) override;

public:
    Eigen::MatrixXd getHam();
    void print() override;
};


#endif // DOCI_DOCI_DENSE_HPP
