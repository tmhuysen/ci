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
#include "FCI.hpp"



#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>



namespace ci {


/*
 *  OVERRIDDEN PRIVATE METHODS
 */

/**
 *  Given a @param matrix_solver, construct the FCI Hamiltonian matrix in the solver's matrix representation.
 */
void FCI::constructHamiltonian(numopt::eigenproblem::BaseMatrixSolver* matrix_solver) {

    // The construction of the FCI Hamiltonian is implemented in three parts: alpha-alpha, beta-beta, and alpha-beta


    // 1. ALPHA-ALPHA
    alpha_one_electron_couplings = new OneElectronCoupling*[dim_alpha];

    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<boost::dynamic_bitset<>> spin_string_alpha (0, this->addressing_scheme_alpha);  // alpha spin string with address 0


    for (size_t Ia = 0; Ia < this->dim_alpha; Ia++) {  // Ia loops over all the addresses of the alpha spin strings

        // (K+1-N_A) * N_A -> N_A annihilations followed by K+1-N_A creations (K-N_A available MO's +1 from the inplace excitation)
        alpha_one_electron_couplings[Ia] = new OneElectronCoupling[((K+1)-N_A)*N_A];  // stores all excitation operator evaluations for each spin string
        size_t eval_index = 0;

        for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            int sign_p = 1;  // sign for the annihilation operator (a_p)

            if (spin_string_alpha.annihilate(p, sign_p)) {
                for (size_t q = 0; q < this->K; q++) {  // q loops over SOs

                    // one-electron contributions for alpha, i.e. one electron excitation
                    int sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                    if (spin_string_alpha.create(q, sign_pq)) {

                        size_t Ja = spin_string_alpha.address(this->addressing_scheme_alpha);

                        // For the 'diagonal beta contributions', i.e. Ib = Jb, the one-electron alpha contributions
                        // are the same
                        // We are storing the alpha addresses as 'major', i.e. the total address IaIb = Ia * dim_b + I_b
                        for (size_t Ib = 0; Ib < this->dim_beta; Ib++) {
                            double value = sign_pq * this->so_basis.get_h_SO(p,q);
                            matrix_solver->addToMatrix(value, Ia * dim_beta + Ib, Ja * dim_beta + Ib);
                        }

                        alpha_one_electron_couplings[Ia][eval_index] = OneElectronCoupling{sign_pq,p,q,Ja};
                        eval_index++;
                        spin_string_alpha.annihilate(q);  // undo the previous creation on q
                    }  // create on q (alpha)


                    // two-electron contributions for beta-beta, i.e. two electron excitations
                    sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                                       // we have to reset this because we changed this in the previous if-statement

                    if (spin_string_alpha.annihilate(q, sign_pq)) {

                        for (size_t r = 0; r < this->K; r++) {
                            int sign_pqr = sign_pq;  // sign for total operator (a^\dagger_r a_q a_p)

                                if (spin_string_alpha.create(r, sign_pqr)) {
                                for (size_t s = 0; s < this->K; s++) {

                                    int sign_pqrs = sign_pqr;  // sign for total operator (a^dagger_s a^\dagger_r a_q a_p)
                                    if (spin_string_alpha.create(s, sign_pqrs)) {

                                        size_t Ja = spin_string_alpha.address(this->addressing_scheme_alpha);

                                        // For the 'diagonal beta contributions', i.e. Ib = Jb, the two-electron alpha
                                        // contributions are the same

                                        // We are storing the alpha addresses as 'major', i.e. the total address IaIb = Ia * dim_b + I_b
                                        for (size_t Ib = 0; Ib < this->dim_beta; Ib++) {
                                            double value = sign_pqrs * 0.5 * this->so_basis.get_g_SO(p,s,q,r);
                                            matrix_solver->addToMatrix(value, Ia * dim_beta + Ib, Ja * dim_beta + Ib);
                                        }

                                        spin_string_alpha.annihilate(s);  // undo the previous creation on s
                                    }  // create on s (alpha)
                                }  // loop over s

                                spin_string_alpha.annihilate(r);  // undo the previous creation on r
                            }  // create on r (alpha)
                        }  // loop over r

                        spin_string_alpha.create(q);  // undo the previous annihilation on q
                    }  // annihilate on q (alpha)
                }  // loop over q

                spin_string_alpha.create(p);  // undo the previous annihilation on p
            }  // annihilate p (alpha)
        }  // loop over p


        if (Ia < dim_alpha - 1) {  // prevent the last permutation to occur
            spin_string_alpha.nextPermutation();
        }
    }  // loop over alpha addresses (Ia)


    // 2. BETA-BETA

    beta_one_electron_couplings = new OneElectronCoupling*[dim_beta];

    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>

    bmqc::SpinString<boost::dynamic_bitset<>> spin_string_beta (0, this->addressing_scheme_beta);  // beta spin string with address 0

    for (size_t Ib = 0; Ib < this->dim_beta; Ib++) {  // Ib loops over addresses of all beta spin strings

        // (K+1-N_B) * N_B -> N_B annihilations followed by K+1-N_B creations (K-N_B available MO's +1 from the inplace excitation)
        beta_one_electron_couplings[Ib] = new OneElectronCoupling[((K+1)-N_B)*N_B];  // stores all excitation operator evaluations for each spin string
        size_t eval_index = 0;

        for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            int sign_p = 1;

            if (spin_string_beta.annihilate(p, sign_p)) {
                for (size_t q = 0; q < this->K; q++) {  // q loops over SOs

                    // one-electron contributions for beta, i.e. one electron excitation
                    int sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                    if (spin_string_beta.create(q, sign_pq)) {

                        size_t Jb = spin_string_beta.address(this->addressing_scheme_beta);

                        // For the 'diagonal alpha contributions', i.e. Ia = Ja, the one-electron beta contributions are
                        // the same
                        // We are storing the alpha addresses as 'major', i.e. the total address IaIb = Ia * dim_b + I_b

                        for (size_t Ia = 0; Ia < this->dim_alpha; Ia++) {
                            double value = sign_pq * this->so_basis.get_h_SO(p,q);
                            matrix_solver->addToMatrix(value, Ia * dim_beta + Ib, Ia * dim_beta + Jb);
                        }

                        beta_one_electron_couplings[Ib][eval_index] = OneElectronCoupling{sign_pq,p,q,Jb};

                        eval_index++;
                        spin_string_beta.annihilate(q);  // undo the previous creation on q
                    }  // create on q (beta)


                    // two-electron contributions for beta-beta, i.e. two electron excitations
                    sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                                       // we have to reset this because we changed this in the previous if-statement

                    if (spin_string_beta.annihilate(q, sign_pq)) {

                        for (size_t r = 0; r < this->K; r++) {  // r loops over SOs
                            int sign_pqr = sign_pq;  // sign for total operator (a^\dagger_r a_q a_p)

                            if (spin_string_beta.create(r, sign_pqr)) {
                                for (size_t s = 0; s < this->K; s++) {  // s loops over SOs

                                    int sign_pqrs = sign_pqr;  // sign for total operator (a^dagger_s a^\dagger_r a_q a_p)
                                    if (spin_string_beta.create(s, sign_pqrs)) {

                                        size_t Jb = spin_string_beta.address(this->addressing_scheme_beta);

                                        // For the 'diagonal alpha contributions', i.e. Ia = Ja, the two-electron beta
                                        // contributions are the same

                                        // We are storing the alpha addresses as 'major', i.e. the total address IaIb = Ia * dim_b + I_b
                                        for (size_t Ia = 0; Ia < this->dim_alpha; Ia++) {
                                            double value = sign_pqrs * 0.5 * this->so_basis.get_g_SO(p,s,q,r);
                                            matrix_solver->addToMatrix(value, Ia * dim_beta + Ib, Ia * dim_beta + Jb);
                                        }

                                        spin_string_beta.annihilate(s);  // undo the previous creation on s
                                    }  // create on s (beta)
                                }  // loop over s

                                spin_string_beta.annihilate(r);  // undo the previous creation on r
                            }  // create on r (beta)
                        }  // loop over r

                        spin_string_beta.create(q);  // undo the previous annihilation on q
                    }  // annihilate on q (beta)
                }  // loop over q

                spin_string_beta.create(p);  // undo the previous annihilation on p
            } // annihilate on p (beta)
        }  // loop over p

        if (Ib < dim_beta - 1) {  // prevent last permutation to occur
            spin_string_beta.nextPermutation();
        }
    }  // loop over beta addresses (Ib)


    // 3. ALPHA-BETA
    for (size_t Ia = 0; Ia < this->dim_alpha; Ia++) {  // loop over alpha addresses
        for (size_t Ib = 0; Ib < this->dim_beta; Ib++) {  // loop over beta addresses

            for(size_t a = 0; a<((K+1)-N_A)*N_A; a++){  // iterate over dimensions of the amount possible single excitations for alpha string
                OneElectronCoupling alpha = alpha_one_electron_couplings[Ia][a];

                for(size_t b = 0; b < ((K+1)-N_B)*N_B; b++){  // iterate over dimensions of the amount possible single excitations for beta string
                    OneElectronCoupling beta = beta_one_electron_couplings[Ib][b];

                    int sign = beta.sign*alpha.sign;
                    //  Construct combined evaluations from single excitations to retrieve two-electron integrals from alpha-beta-beta-alpha operations.
                    matrix_solver->addToMatrix(sign*this->so_basis.get_g_SO(alpha.p, alpha.q, beta.p, beta.q), Ia * dim_beta + Ib, alpha.address * dim_beta + beta.address);
                }
            }

        }  // loop over beta addresses (Ib)
    }  // loop over alpha addresses (Ia)
}


/**
 *  @return the action of the FCI Hamiltonian of the coefficient vector @param x.
 */
Eigen::VectorXd FCI::matrixVectorProduct(const Eigen::VectorXd& x) {
    throw std::logic_error("This function hasn't been implemented yet.");
}


/**
 *  @set the diagonal of the matrix representation of the FCI Hamiltonian.
 */
void FCI::calculateDiagonal() {


    // Initialize the diagonal
    this->diagonal = Eigen::VectorXd::Zero(this->dim);


    // Calculate the effective one-electron integrals
    // TODO: move this to libwint
    Eigen::MatrixXd k_SO = this->so_basis.get_h_SO();
    for (size_t p = 0; p < this->K; p++) {
        for (size_t q = 0; q < this->K; q++) {
            for (size_t r = 0; r < this->K; r++) {
                k_SO(p,q) -= 0.5 * this->so_basis.get_g_SO(p,r,r,q);
            }
        }
    }


    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string_alpha (0, this->addressing_scheme_alpha);
    for (size_t Ia = 0; Ia < this->dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        bmqc::SpinString<unsigned long> spin_string_beta (0, this->addressing_scheme_beta);
        for (size_t Ib = 0; Ib < this->dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

            for (size_t p = 0; p < this->K; p++) {  // p loops over SOs

                if (spin_string_alpha.isOccupied(p)) {  // p is in Ia
                    this->diagonal(Ia * this->dim_beta + Ib) += k_SO(p, p);

                    for (size_t q = 0; q < this->K; q++) {  // q loops over SOs
                        if (spin_string_alpha.isOccupied(q)) {  // q is in Ia
                            this->diagonal(Ia * this->dim_beta + Ib) += 0.5 * this->so_basis.get_g_SO(p, p, q, q);
                        } else {  // q is not in I_alpha
                            this->diagonal(Ia * this->dim_beta + Ib) += 0.5 * this->so_basis.get_g_SO(p, q, q, p);
                        }

                        if (spin_string_beta.isOccupied(q)) {  // q is in Ib
                            this->diagonal(Ia * this->dim_beta + Ib) += this->so_basis.get_g_SO(p, p, q, q);
                        }
                    }  // q loop
                }


                if (spin_string_beta.isOccupied(p)) {  // p is in Ib
                    this->diagonal(Ia * this->dim_beta + Ib) += k_SO(p, p);


                    for (size_t q = 0; q < this->K; q++) {  // q loops over SOs
                        if (spin_string_beta.isOccupied(q)) {  // q is in Ib
                            this->diagonal(Ia * this->dim_beta + Ib) += 0.5 * this->so_basis.get_g_SO(p, p, q, q);

                        } else {  // q is not in I_beta
                            this->diagonal(Ia * this->dim_beta + Ib) += 0.5 * this->so_basis.get_g_SO(p, q, q, p);
                        }
                    }  // q loop
                }

            }  // p loop

            if (Ib < this->dim_beta - 1) {  // prevent last permutation to occur
                spin_string_beta.nextPermutation();
            }
        }  // beta address (Ib) loop

        if (Ia < this->dim_alpha - 1) {  // prevent last permutation to occur
            spin_string_alpha.nextPermutation();
        }
    }  // alpha address (Ia) loop
}



/*
 *  CONSTRUCTORS
 */
/**
 *  Constructor based on a given @param so_basis, a number of alpha electrons @param N_A and a number of beta electrons
 *  @param N_B.
 */
FCI::FCI(libwint::SOBasis& so_basis, size_t N_A, size_t N_B) :
        BaseCI(so_basis, ci::FCI::calculateDimension(so_basis.get_K(), N_A, N_B)),
        K (so_basis.get_K()),
        N_A (N_A),
        N_B (N_B),
        addressing_scheme_alpha (bmqc::AddressingScheme(this->K, this->N_A)),
        addressing_scheme_beta (bmqc::AddressingScheme(this->K, this->N_B)),
        dim_alpha (ci::FCI::calculateDimension(so_basis.get_K(), N_A, 0)),
        dim_beta (ci::FCI::calculateDimension(so_basis.get_K(), 0, N_B))
{
    // Do some input checks
    if (this->K < this->N_A) {
        throw std::invalid_argument("Too many alpha-electrons to place in the number of spatial orbitals.");
    }

    if (this->K < this->N_B) {
        throw std::invalid_argument("Too many beta-electrons to place in the number of spatial orbitals.");
    }
}



/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  Given a number of spatial orbitals @param K, a number of alpha electrons @param N_A, and a number of beta electrons
 *  @param N_B, @return the dimension of the FCI space.
 */
size_t FCI::calculateDimension(size_t K, size_t N_A, size_t N_B) {

    // K, N_A, N_B are expected to be small, so static-casting them to unsigned (what boost needs) is permitted.
    auto dim_double_alpha = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_A));
    auto dim_double_beta = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_B));
    auto dim_double = dim_double_alpha * dim_double_beta;

    // Check if the resulting dimension is appropriate to be stored in size_t
    return boost::numeric::converter<double, size_t>::convert(dim_double);
}

void FCI::calculate2RDMs() {

}

void FCI::calculate1RDMs() {

}


/*
 *  PUBLIC METHODS
 */

/**
 *  Calculate all the 1-RDMs.
 */
void FCI::calculate1RDMs() {
    throw std::logic_error("This function hasn't been implemented yet.");
}


/**
 *  Calculate all the 2-RDMs.
 */
void FCI::calculate2RDMs() {
    throw std::logic_error("This function hasn't been implemented yet.");

}



}  // namespace ci
