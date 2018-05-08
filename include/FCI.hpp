// This file is part of GQCG-ci.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-ci is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-ci is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-ci.  If not, see <http://www.gnu.org/licenses/>.
#ifndef CI_FCI_HPP
#define CI_FCI_HPP


#include "BaseCI.hpp"

typedef numopt::eigenproblem::BaseMatrixSolver MatrixSolver;

namespace ci {

struct SpinEvaluations{
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

    // Rectangular Matrix of Spin Evaluations
    SpinEvaluations** alpha_evaluation;  // Creation, annihilation pair evaluations of alpha spin strings.
    SpinEvaluations** beta_evaluation;  // Creation, annihilation pair evaluations of beta spin strings.


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
     *  @set the diagonal of the matrix representation of the FCI Hamiltonian.
     */
    void calculateDiagonal() override;

    /**
     * Performs all alpha related operations
     */
    void alpha_branch(MatrixSolver* matrix_solver);

    /**
     * Performs all beta related operations
     */
    void beta_branch(MatrixSolver* matrix_solver);

    /**
     * Recombines all alpha and beta single excitation related operations
     * That were stored in the previous branches resulting
     */
    void mixed_branch(MatrixSolver* matrix_solver);



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
     *  Calculate all the 1-RDMs for FCI.
     */
    void calculate1RDMs() override;

    /**
     *  Calculate all the 2-RDMS for FCI.
     */
    void calculate2RDMs() override;


    // STATIC PUBLIC METHODS
    /**
     *  Given a number of spatial orbitals @param K and a number of alpha electrons @param N_A and beta electrons @param N_B, @return the dimension of
     *  the FCI space.
     */
    static size_t calculateDimension(size_t K, size_t N_A, size_t N_B);


    // PUBLIC METHODS
    /**
     *  Calculate all the 1-RDMs.
     */
    void calculate1RDMs() override;

    /**
     *  Calculate all the 2-RDMS.
     */
    void calculate2RDMs() override;
};


}  // namespace ci



#endif //CI_FCI_HPP
