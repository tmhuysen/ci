#ifndef DOCI_DOCI_DENSE_HPP
#define DOCI_DOCI_DENSE_HPP


#include "DOCI_Class.hpp"


class DenseDOCI : public DOCI {  // FIXME: add documentation/comments

private:
    Eigen::MatrixXd hamiltonian;


// Protected methods
protected:
    void addToHamiltonian(double value, size_t index1, size_t index2) override;


public:
    DenseDOCI(doci::CI_basis basis);

    Eigen::MatrixXd getHam();
    void print() override;
};


#endif // DOCI_DOCI_DENSE_HPP
