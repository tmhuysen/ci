#ifndef CI_BASECI_HPP
#define CI_BASECI_HPP



#include <libwint.hpp>
#include <bmqc.hpp>

#include "numopt.hpp"



namespace ci {



class BaseCI {
protected:

    libwint::SOBasis& so_basis;
    numopt::eigenproblem::BaseEigenproblemSolver* eigensolver_ptr = nullptr;

    const size_t dim;  // the dimension of the CI space

    bool are_computed_one_rdms = false;
    bool are_computed_two_rdms = false;

    Eigen::MatrixXd one_rdm_aa;  // 1-RDM for alpha-alpha (a-a)
    Eigen::MatrixXd one_rdm_bb;  // 1-RDM for beta-beta
    Eigen::Tensor<double, 4> two_rdm_aaaa;  // 2-RDM for a-a-a-a
    Eigen::Tensor<double, 4> two_rdm_abba;  // 2-RDM for a-b-b-a
    Eigen::Tensor<double, 4> two_rdm_baab;  // 2-RDM for b-a-a-b
    Eigen::Tensor<double, 4> two_rdm_bbbb;  // 2-RDM for b-b-b-b


    // PROTECTED CONSTRUCTORS
    /**
     *  Protected constructor given a @param so_basis and a dimension @dim.
     */
    explicit BaseCI(libwint::SOBasis& so_basis, size_t dim);


    // PURE VIRTUAL PROTECTED METHODS
    /**
     *  Given a @param matrix_solver, construct the Hamiltonian matrix in the solver's matrix representation.
     */
    virtual void constructHamiltonian(numopt::eigenproblem::BaseMatrixSolver* matrix_solver) = 0;

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
     *  Initialize and subsequently solve the eigenvalue problem associated to the derived CI class, using a @param
     *  matrix_solver_ptr.
     */
    void solveMatrixEigenvalueProblem(numopt::eigenproblem::BaseMatrixSolver* matrix_solver_ptr);



public:
    // DESTRUCTOR
    virtual ~BaseCI();


    // GETTERS
    double get_eigenvalue() const { return this->eigensolver_ptr->get_eigenvalue(); }
    Eigen::VectorXd get_eigenvector() const { return this->eigensolver_ptr->get_eigenvector(); }
    Eigen::MatrixXd get_one_rdm_aa() const;
    Eigen::MatrixXd get_one_rdm_bb() const;
    Eigen::Tensor<double, 4> get_two_rdm_aaaa() const;
    Eigen::Tensor<double, 4> get_two_rdm_abba() const;
    Eigen::Tensor<double, 4> get_two_rdm_baab() const;
    Eigen::Tensor<double, 4> get_two_rdm_bbbb() const;


    // PUBLIC METHODS
    /**
     *  Find the lowest energy eigenpair of the Hamiltonian, using a @param solver_type.
     */
    void solve(numopt::eigenproblem::SolverType solver_type);

    /**
     *  Compute all the 1-RDMS.
     */
    virtual void compute1RDMs() = 0;

    /**
     *  Compute all the 2-RDMS.
     */
    virtual void compute2RDMs() = 0;
};



}  // namespace ci



#endif  // CI_BASECI_HPP
