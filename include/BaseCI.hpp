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
#ifndef CI_BASECI_HPP
#define CI_BASECI_HPP



#include <libwint.hpp>
#include <bmqc.hpp>

#include "numopt.hpp"



namespace ci {



class BaseCI {
protected:

    libwint::SOMullikenBasis& so_basis;
    numopt::eigenproblem::BaseEigenproblemSolver* eigensolver_ptr = nullptr;
    double population_set;

    const size_t dim;  // the dimension of the CI space

    Eigen::VectorXd diagonal;  // the diagonal of the Hamiltonian matrix

    bool are_computed_one_rdms = false;
    bool are_computed_two_rdms = false;

    Eigen::MatrixXd one_rdm_aa;  // alpha-alpha (a-a) 1-RDM
    Eigen::MatrixXd one_rdm_bb;  // beta-beta (b-b) 1-RDM
    Eigen::MatrixXd one_rdm;  // spin-summed (total) 1-RDM

    Eigen::Tensor<double, 4> two_rdm_aaaa;  // a-a-a-a 2-RDM
    Eigen::Tensor<double, 4> two_rdm_aabb;  // a-a-b-b 2-RDM
    Eigen::Tensor<double, 4> two_rdm_bbaa;  // b-a-a-b 2-RDM
    Eigen::Tensor<double, 4> two_rdm_bbbb;  // b-b-b-b 2-RDM
    Eigen::Tensor<double, 4> two_rdm;  // spin-summed (total) 2-RDM


    // PROTECTED CONSTRUCTORS
    /**
     *  Protected constructor given a @param so_basis and a dimension @dim.
     */
    explicit BaseCI(libwint::SOMullikenBasis& so_basis, size_t dim);


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
     *  @set the diagonal of the matrix representation of the Hamiltonian.
     */
    virtual void calculateDiagonal() = 0;


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
    size_t get_dim() const { return this->dim; }
    double get_eigenvalue() const { return this->eigensolver_ptr->get_eigenvalue(); }
    Eigen::VectorXd get_eigenvector() const { return this->eigensolver_ptr->get_eigenvector(); }
    double get_population_set() const { return this->population_set; }
    libwint::SOMullikenBasis get_basis() { return this->so_basis;}

    Eigen::MatrixXd get_one_rdm_aa() const;
    Eigen::MatrixXd get_one_rdm_bb() const;
    Eigen::MatrixXd get_one_rdm() const;

    Eigen::Tensor<double, 4> get_two_rdm_aaaa() const;
    Eigen::Tensor<double, 4> get_two_rdm_aabb() const;
    Eigen::Tensor<double, 4> get_two_rdm_bbaa() const;
    Eigen::Tensor<double, 4> get_two_rdm_bbbb() const;
    Eigen::Tensor<double, 4> get_two_rdm() const;


    //SETTERS


    //void set_constraint_base(libwint::SOMullikenBasis& so_basis){this->so_basis=so_basis;};


    // PUBLIC METHODS
    /**
     *  Find the lowest energy eigenpair of the Hamiltonian, using a @param solver_type.
     */
    void solve(numopt::eigenproblem::SolverType solver_type);

    /**
     *  Solves the eigenvalue problem with a for a lagrange multiplier mulliken constraint and returns the energy.
     */
    double solveConstrained(numopt::eigenproblem::SolverType solver_type, std::vector<size_t> AO_set, double multiplier);


    /**
     *  Calculate all the 1-RDMs.
     */
    virtual void calculate1RDMs() = 0;

    /**
     *  Calculate all the 2-RDMS.
     */
    virtual void calculate2RDMs() = 0;
};



}  // namespace ci



#endif  // CI_BASECI_HPP
