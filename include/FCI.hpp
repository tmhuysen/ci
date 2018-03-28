#ifndef CI_FCI_HPP
#define CI_FCI_HPP


#include "BaseCI.hpp"



namespace ci {


class FCI : public ci::BaseCI {
private:
    const size_t dim_alpha;  // the dimension of the alpha CI space
    const size_t dim_beta;  // the dimension of the beta CI space
    const size_t K;  // number of spatial orbitals
    const size_t N_A;  // number of alpha electrons
    const size_t N_B;  // number of beta electrons
    const bmqc::AddressingScheme addressing_scheme_alpha;
    const bmqc::AddressingScheme addressing_scheme_beta;


    // OVERRIDDEN PRIVATE METHODS
    /**
     *  Given a @param matrix_solver, construct the FCI Hamiltonian matrix in the solver's matrix representation.
     */
    void constructHamiltonian(numopt::eigenproblem::BaseMatrixSolver* matrix_solver) override;

    /**
     *  @return the action of the FCI Hamiltonian of the coefficient vector @param x.
     */
    Eigen::VectorXd matrixVectorProduct(const Eigen::VectorXd& x) override;

    /**
     *  @return the diagonal of the matrix representation of the FCI Hamiltonian.
     */
    Eigen::VectorXd calculateDiagonal() override;



public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param so_basis and a number of alpha electron and beta electrons @param N_A and N_B respectively.
     */
    FCI(libwint::SOBasis& so_basis, size_t N_A, size_t N_B);


    // DESTRUCTOR
    ~FCI() override = default;


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  Compute all of the one reduced density matrix.
     */
    void compute1RDM() override;

    /**
     *  Compute all of the two reduced density matrix.
     */
    void compute2RDM()override;


    // STATIC PUBLIC METHODS
    /**
     *  Given a number of spatial orbitals @param K and a number of alpha electrons @param N_A and beta electrons @param N_B, @return the dimension of
     *  the FCI space.
     */
    static size_t calculateDimension(size_t K, size_t N_A, size_t N_B);
};


}  // namespace ci



#endif //CI_FCI_HPP
