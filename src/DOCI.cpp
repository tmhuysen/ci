#include "DOCI.hpp"



#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace ci {


/*
 *  OVERRIDDEN PRIVATE METHODS
 */

/**
 *  Given a @param matrix_solver, construct the DOCI Hamiltonian matrix in the solver's matrix representation. An
 *  @param addressing_scheme is used to speed up the searches for the addresses of coupling spin strings.
 */
void DOCI::constructHamiltonian(ci::solver::BaseMatrixSolver* matrix_solver, const bmqc::AddressingScheme& addressing_scheme) {

}


/**
 *  Given a @param addressing_scheme, @return the action of the DOCI Hamiltonian of the coefficient vector @param x.
 */
Eigen::VectorXd DOCI::matrixVectorProduct(const bmqc::AddressingScheme& addressing_scheme, const Eigen::VectorXd& x) {

}


/**
 *  @return the diagonal of the matrix representation of the DOCI Hamiltonian.
 */
Eigen::VectorXd DOCI::calculateDiagonal() {

}



/*
 *  CONSTRUCTORS
 */
DOCI::DOCI(libwint::SOBasis& so_basis, const libwint::Molecule& molecule) :
    N_P (molecule.get_N() / 2),
    K (so_basis.get_K()),
    BaseCI(this->calculateDimension(this->K, this->N_P), so_basis)
{

    // Do some input checks
    if ((molecule.get_N() % 2) != 0) {
        throw std::invalid_argument("Your basis contains an odd amount of electrons and is not suitable for DOCI");
    }

    if (this->K < this->N_P) {
        throw std::invalid_argument("Too many electrons to place into the given number of spatial orbitals");
    }

}


/*
 *  STATIC PUBLIC METHODS
 */
/**
 *  Given a number of spatial orbitals @param K and a number of electron pairs @param N_P, @return the dimension of
 *  the DOCI space.
 */
size_t DOCI::calculateDimension(size_t K, size_t N_P) {

    // K and N_P are expected to be small, so static-casting them to unsigned (what boost needs) is permitted.
    auto dim_double = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_P));


    // Check if the resulting dimension is appropriate to be stored in size_t
    return boost::numeric::converter<double, size_t>::convert(dim_double);
}


}  // namespace ci
