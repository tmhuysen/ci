#ifndef CI_BASECI_HPP
#define CI_BASECI_HPP



#include <libwint.hpp>
#include <bmqc.hpp>

#include "BaseMatrixSolver.hpp"



namespace ci {



class BaseCI {
protected:
    const size_t dim;  // the dimension of the Fock space (i.e. the dimension of the CI space)

    libwint::SOBasis& so_basis;


    // PURE VIRTUAL PROTECTED METHODS
    /**
     *  Given a @param matrix_solver, construct the Hamiltonian matrix in the solver's matrix representation. An
     *  @param addressing_scheme is used to speed up the searches for the addresses of coupling spin strings.
     */
    virtual void constructHamiltonian(ci::solver::BaseMatrixSolver* matrix_solver, const bmqc::AddressingScheme& addressing_scheme) = 0;

    /**
     *  Given a @param addressing_scheme, @return the action of the Hamiltonian of the coefficient vector @param x.
     */
    virtual Eigen::VectorXd matrixVectorProduct(const bmqc::AddressingScheme& addressing_scheme, const Eigen::VectorXd& x) = 0;

    /**
     *  @return the diagonal of the matrix representation of the Hamiltonian.
     */
    virtual Eigen::VectorXd calculateDiagonal() = 0;


public:
    // PUBLIC METHODS
    /**
     *  Find the lowest energy eigenpair of the Hamiltonian.
     */
    void solve(ci::solver::SolverType solver_type, const bmqc::AddressingScheme& addressing_scheme);


};



}  // namespace ci



#endif  // CI_BASECI_HPP
