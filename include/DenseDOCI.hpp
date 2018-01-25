#ifndef DOCI_DOCI_DENSE_HPP
#define DOCI_DOCI_DENSE_HPP


#include "DOCI_Class.hpp"


namespace doci {
/**
 * Dense DOCI for calculations where the hamiltonian is stored in a dense matrix from the eigen lib
 * Use only
 */
class DenseDOCI : public doci::DOCI {


private:
    Eigen::MatrixXd hamiltonian;


// Protected methods
protected:
    /**
     * function that stores a calculated value in the Hamiltonian
     * Us
     */
    void addToHamiltonian(double value, size_t index1, size_t index2) override;


public:
    DenseDOCI(doci::CI_basis basis);

    Eigen::MatrixXd getHam();

    void print() override;
};

} // namespace doci

#endif // DOCI_DOCI_DENSE_HPP
