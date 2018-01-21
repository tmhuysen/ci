#ifndef DOCI_CI_BASIS_HPP
#define DOCI_CI_BASIS_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <hf.hpp>


namespace doci {

class CI_basis {
private:
    Eigen::MatrixXd one_ints;  // The one-electron integrals
    Eigen::Tensor<double, 4> two_ints;  // The two-electron integrals
    double internuclear_repulsion;  // The internuclear repulsion energy
    size_t nbf;  // The number of basis functions
    size_t nelec;  // The number of electrons

public:

    /** Constructor based on a given RHF instance
     */
    CI_basis(hf::rhf::RHF& rhf);

    /** Constructor based on a given filename
     */
    CI_basis(const std::string& filename);

};

} // namespace doci


#endif // DOCI_CI_BASIS_HPP
