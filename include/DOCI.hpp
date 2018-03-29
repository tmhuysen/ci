#ifndef CI_DOCI_HPP
#define CI_DOCI_HPP



#include "BaseCI.hpp"



namespace ci {


class DOCI : public ci::BaseCI {
private:
    const size_t K;  // number of spatial orbitals
    const size_t N_P;  // number of electron pairs
    const bmqc::AddressingScheme addressing_scheme;


    // OVERRIDDEN PRIVATE METHODS
    /**
     *  Given a @param matrix_solver, construct the DOCI Hamiltonian matrix in the solver's matrix representation.
     */
    void constructHamiltonian(numopt::eigenproblem::BaseMatrixSolver* matrix_solver) override;

    /**
     *  @return the action of the DOCI Hamiltonian of the coefficient vector @param x.
     */
    Eigen::VectorXd matrixVectorProduct(const Eigen::VectorXd& x) override;

    /**
     *  @return the diagonal of the matrix representation of the DOCI Hamiltonian.
     */
    Eigen::VectorXd calculateDiagonal() override;



public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param so_basis and a number of electrons @param N.
     */
    DOCI(libwint::SOMullikenBasis& so_basis, size_t N);


    // DESTRUCTOR
    ~DOCI() override = default;


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  Compute all of the one reduced density matrix.
     */
    void compute1RDM() override;

    /**
     *  Compute all of the two reduced density matrix.
     */
    void compute2RDM()override;

    /**
     * test
     */
    double computeE();
    // STATIC PUBLIC METHODS
    /**
     *  Given a number of spatial orbitals @param K and a number of electron pairs @param N_P, @return the dimension of
     *  the DOCI space.
     */
    static size_t calculateDimension(size_t K, size_t N_P);
};


}  // namespace ci



#endif  // CI_DOCI_HPP
