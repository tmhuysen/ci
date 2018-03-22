#ifndef CI_BASECI_HPP
#define CI_BASECI_HPP



#include <libwint.hpp>
#include <bmqc.hpp>

#include "numopt.hpp"



namespace ci {



class BaseCI {
protected:
    libwint::SOBasis& so_basis;

    // Every derived class should have the following two members. However, until runtime we don't know anything about
    // their exact value.
    size_t dim = 0;  // a good value for an uninitialized state, since the real dimension is always greater than 0
    numopt::eigenproblem::BaseEigenproblemSolver* eigensolver_ptr = nullptr;

    bmqc::AddressingScheme* addressing_scheme_ptr = nullptr;

    bool is_solved = false;
    double eigenvalue;
    Eigen::VectorXd eigenvector;


    // PROTECTED CONSTRUCTORS
    /**
     *  Protected constructor to initialize the reference @member so_basis by @param so_basis and the @param
     *  eigensolver_ptr.
     */
    explicit BaseCI(libwint::SOBasis& so_basis, numopt::eigenproblem::SolverType solver_type);


    // PURE VIRTUAL PROTECTED METHODS
    /**
     *  Given a @param matrix_solver, construct the Hamiltonian matrix in the solver's matrix representation.
     */
    virtual void constructHamiltonian() = 0;

    /**
     *  @return the action of the Hamiltonian of the coefficient vector @param x.
     */
    virtual Eigen::VectorXd matrixVectorProduct(const Eigen::VectorXd& x) = 0;

    /**
     *  @return the diagonal of the matrix representation of the Hamiltonian.
     */
    virtual Eigen::VectorXd calculateDiagonal() = 0;


    // PROTECTED METHODS
    /**
     *  Initialize and subsequently solve the eigenvalue problem associated to the derived CI class, using a pointer
     *  to a @param matrix_solver.
     */
    void solveMatrixEigenvalueProblem(ci::solver::BaseMatrixSolver* matrix_solver);



public:
    // DESTRUCTOR
    virtual ~BaseCI();


    // GETTERS
    double get_eigenvalue() const { return this->eigensolver_ptr->get_eigenvalue(); }
    Eigen::VectorXd get_eigenvector() const { return this->eigensolver_ptr->get_eigenvector(); }


    // PUBLIC METHODS
    /**
     *  Find the lowest energy eigenpair of the Hamiltonian, using a @param solver_type.
     */
    void solve(ci::solver::SolverType solver_type);
};



}  // namespace ci



#endif  // CI_BASECI_HPP
