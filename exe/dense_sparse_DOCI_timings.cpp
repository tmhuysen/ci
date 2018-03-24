/**
 *  Do some timings concerning for dense and sparse DOCI calculations
 */


#include <chrono>
#include <iostream>

#include <Eigen/Dense>
#include <cpputil.hpp>

#include "numopt.hpp"



/*
 *  HELPER FUNCTIONS
 */

/**
 *  Given a dimension @param dim, construct an example matrix. The example matrix can be found in liu1978
 */
Eigen::MatrixXd constructLiuExampleMatrix(size_t dim) {

    Eigen::MatrixXd A = Eigen::MatrixXd::Ones(dim, dim);
    for (size_t i = 0; i < dim; i++) {
        if (i < 5) {
            A(i, i) = 1 + 0.1 * i;
        } else {
            A(i, i) = 2 * (i + 1) - 1;
        }
    }

    return A;
}


/**
 *  Solve the full eigenvalue problem for a matrix @param A using the SelfAdjointEigenSolver.
 */
void solveEigenvalueProblemWithEigen(const Eigen::MatrixXd& A) {

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver (A);
}


/**
 *  Find the lowest eigenpair of a matrix @param A using Davidson's diagonalization method.
 */
void solveEigenvalueProblemWithDavidson(const Eigen::MatrixXd& A) {

    Eigen::VectorXd t_0 = Eigen::VectorXd::Zero(A.cols());
    t_0(0) = 1;
    numopt::eigenproblem::DavidsonSolver davidson_solver (A, t_0);
    davidson_solver.solve();
}


/**
 *  Print how long it takes to solve the full eigenvalue problem for a matrix @param A using Eigen.
 */
void printEigenTimings(const Eigen::MatrixXd& A) {

    // Since our timer wrapper expects a std::function<void ()>, we should 'bind' the argument A to the function call
    // We can use this easily using a lambda function
    std::function<void ()> callable_eigen_function = [A] { solveEigenvalueProblemWithEigen(A); };


    cpputil::printExecutionTime("Eigen", callable_eigen_function);
}


/**
 *  Print how long it takes to find the lowest eigenpair of a matrix @param A using Davidson's diagonalization method.
 */
void printDavidsonTimings(const Eigen::MatrixXd& A) {

    // Since our timer wrapper expects a std::function<void ()>, we should 'bind' the argument A to the function call
    // We can use this easily using a lambda function
    std::function<void ()> callable_eigen_function = [A] { solveEigenvalueProblemWithDavidson(A); };


    cpputil::printExecutionTime("Davidson", callable_eigen_function);
}



/*
 *  MAIN FUNCTION
 */

int main () {

    Eigen::MatrixXd A = constructLiuExampleMatrix(1000);


    printEigenTimings(A);
    printDavidsonTimings(A);
}