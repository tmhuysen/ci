/**
 *  Do some timings concerning for dense and sparse DOCI calculations
 */


#include <chrono>
#include <iostream>

#include <Eigen/Dense>
#include <cpputil.hpp>

#include <numopt.hpp>
#include <libwint.hpp>
#include <ci.hpp>



/*
 *  HELPER FUNCTIONS
 */


/**
 *  Solve a dense DOCI problem (full diagonalization of the DOCI Hamiltonian matrix).
 */
void solveDenseDoci() {

    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../../tests/reference_data/beh_cation_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons
    doci.solve(numopt::eigenproblem::SolverType::DENSE);


    // Calculate the total energy
    double internuclear_repulsion_energy = 1.5900757460937498e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;
}


/**
 *  Solve a sparse DOCI problem (find the lowest eigenpair of the DOCI Hamiltonian matrix).
 */
void solveSparseDoci() {

    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../../tests/reference_data/beh_cation_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons
    doci.solve(numopt::eigenproblem::SolverType::SPARSE);


    // Calculate the total energy
    double internuclear_repulsion_energy = 1.5900757460937498e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;
}


/**
 *  Print how long it takes to solve the dense DOCI problem.
 */
void printDenseTimings() {

    cpputil::printExecutionTime("Dense DOCI", solveDenseDoci);  // (void *)() is implicitly converted to std::function<void ()>
}


/**
 *  Print how long it takes to solve the sparse DOCI problem.
 */
void printSparseTimings() {

    cpputil::printExecutionTime("Sparse DOCI", solveSparseDoci);  // (void *)() is implicitly converted to std::function<void ()>
}



/*
 *  MAIN FUNCTION
 */

int main () {

    printDenseTimings();
    printSparseTimings();
}
