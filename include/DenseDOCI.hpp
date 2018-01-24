#ifndef DOCI_DOCI_DENSE_HPP
#define DOCI_DOCI_DENSE_HPP


#include "DOCI_Class.hpp"


namespace doci {

class DenseDOCI : public doci::DOCI {  // FIXME: add documentation/comments

private:
    Eigen::MatrixXd hamiltonian;


// Protected methods
protected:
    /**
     * function that adds
     */
    void addToHamiltonian(double value, size_t index1, size_t index2) override;


public:
    DenseDOCI(doci::CI_basis basis);

    Eigen::MatrixXd getHam();

    void print() override;
};

} // namespace doci

#endif // DOCI_DOCI_DENSE_HPP
