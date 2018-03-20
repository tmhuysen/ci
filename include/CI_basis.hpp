#ifndef CI_CI_BASIS_HPP
#define CI_CI_BASIS_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <hf.hpp>
#include "utility.hpp"

namespace doci {

class CI_basis {

protected:
    Eigen::MatrixXd one_ints; // The one-electron integrals
    Eigen::Tensor<double, 4> two_ints;  // The two-electron integrals
    double internuclear_repulsion;  // The internuclear repulsion energy

    size_t K;  // The number of spatial orbitals
    size_t nelec;  // The number of electrons

public:
    /** Default constructor
     */
    CI_basis();

    /** Constructor based on a given RHF instance
     */
    CI_basis(hf::rhf::RHF& rhf);

    /**
     * Constructor based on a given FCIDUMP file
     */
    CI_basis(const std::string& filename);

    /**
     * Apply a Jacobi rotation on the CI_basis
     * @param rot in degrees
     * @param index1 row1 you want to apply
     * @param index2 row2 you want to apply
     */
    void rotate(double rot, size_t index1, size_t index2);

    /**
     * Getters
     */

    virtual double getOne_int(size_t index1, size_t index2) const;

    virtual double getTwo_int(size_t index1, size_t index2, size_t index3, size_t index4) const;

    double getInternuclear_repulsion() const;

    size_t getK() const;

    size_t getNelec() const;
};

} // namespace doci


#endif // CI_CI_BASIS_HPP
