#ifndef CI_FCI_HPP
#define CI_FCI_HPP


#include "BaseCI.hpp"

typedef numopt::eigenproblem::BaseMatrixSolver MatrixSolver;

namespace ci {

struct spin_evaluation{
    int sign;
    size_t p;
    size_t q;
    size_t address;
};


class FCI : public ci::BaseCI {
private:
    const size_t dim_alpha;  // the dimension of the alpha CI space
    const size_t dim_beta;  // the dimension of the beta CI space
    const size_t K;  // number of spatial orbitals
    const size_t N_A;  // number of alpha electrons
    const size_t N_B;  // number of beta electrons
    const bmqc::AddressingScheme addressing_scheme_alpha;
    const bmqc::AddressingScheme addressing_scheme_beta;

    spin_evaluation** alpha_evaluation;
    spin_evaluation** beta_evaluation;


    // OVERRIDDEN PRIVATE METHODS
    /**
     *  Given a @param matrix_solver, construct the FCI Hamiltonian matrix in the solver's matrix representation.
     */
    void constructHamiltonian(MatrixSolver* matrix_solver) override;

    /**
     *  @return the action of the FCI Hamiltonian of the coefficient vector @param x.
     */
    Eigen::VectorXd matrixVectorProduct(const Eigen::VectorXd& x) override;

    /**
     *  @return the diagonal of the matrix representation of the FCI Hamiltonian.
     */
    Eigen::VectorXd calculateDiagonal() override;

    /**
     * Alpha en Beta parts or branches of the FCI calculation
     */
    void alpha_branch(MatrixSolver* matrix_solver);
    void beta_branch(MatrixSolver* matrix_solver);
    void mixed_branch(MatrixSolver* matrix_solver);



public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param so_basis and a number of alpha electron and beta electrons @param N_A and N_B respectively.
     */
    FCI(libwint::SOBasis& so_basis, size_t N_A, size_t N_B);


    // DESTRUCTOR
    ~FCI() override = default;


    // STATIC PUBLIC METHODS
    /**
     *  Given a number of spatial orbitals @param K and a number of alpha electrons @param N_A and beta electrons @param N_B, @return the dimension of
     *  the FCI space.
     */
    static size_t calculateDimension(size_t K, size_t N_A, size_t N_B);
};


}  // namespace ci



#endif //CI_FCI_HPP
