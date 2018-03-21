#ifndef CI_BASEMATRIXSOLVER_HPP
#define CI_BASEMATRIXSOLVER_HPP



#include <cstddef>

#include <Eigen/Dense>



namespace ci {
namespace solver {


class BaseMatrixSolver {
protected:
    const size_t dim;  // the dimension of the Hamiltonian matrix representation

    bool is_solved = false;
    double eigenvalue;
    Eigen::VectorXd eigenvector;


    // PROTECTED CONSTRUCTORS
    /**
     *  Protected constructor to initialize the const @member dim by @param dim.
     */
    explicit BaseMatrixSolver(size_t dim);



public:
    // DESTRUCTOR
    virtual ~BaseMatrixSolver() = default;


    // PURE VIRTUAL METHODS
    /**
     *  Solve the eigenvalue problem associated to the matrix solver.
     */
    virtual void solve() = 0;

    /**
     *  Add @param value to the matrix representation of the Hamiltonian at (@param index1, @param index2).
     */
    virtual void addToMatrix(double value, size_t index1, size_t index2) = 0;


    // GETTERS
    double get_eigenvalue() const;
    Eigen::VectorXd get_eigenvector() const;
};


}  // namespace solver
}  // ci



#endif  // CI_BASEMATRIXSOLVER_HPP
