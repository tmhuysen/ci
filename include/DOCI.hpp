#ifndef CI_DOCI_HPP
#define CI_DOCI_HPP



#include "BaseCI.hpp"



namespace ci {


class DOCI : public ci::BaseCI {
private:
    const size_t N_P;  // number of electron pairs
    const size_t K;  // number of spatial orbitals


    // OVERRIDDEN PRIVATE METHODS
    /**
     *  Given a @param matrix_solver, construct the DOCI Hamiltonian matrix in the solver's matrix representation. An
     *  @param addressing_scheme is used to speed up the searches for the addresses of coupling spin strings.
     */
    void constructHamiltonian(ci::solver::BaseMatrixSolver* matrix_solver, const bmqc::AddressingScheme& addressing_scheme) override;

    /**
     *  Given a @param addressing_scheme, @return the action of the DOCI Hamiltonian of the coefficient vector @param x.
     */
    Eigen::VectorXd matrixVectorProduct(const bmqc::AddressingScheme& addressing_scheme, const Eigen::VectorXd& x) override;

    /**
     *  @return the diagonal of the matrix representation of the DOCI Hamiltonian.
     */
    Eigen::VectorXd calculateDiagonal() override;



public:
    // DESTRUCTOR
    ~DOCI() override = default;


    // CONSTRUCTORS
    DOCI(libwint::SOBasis& so_basis, const libwint::Molecule& molecule);


    // STATIC PUBLIC METHODS
    /**
     *  Given a number of spatial orbitals @param K and a number of electron pairs @param N_P, @return the dimension of
     *  the DOCI space.
     */
    static size_t calculateDimension(size_t K, size_t N_P);
};


}  // namespace ci



#endif  // CI_DOCI_HPP
