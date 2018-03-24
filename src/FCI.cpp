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
    alpha_branch(matrix_solver);
    beta_branch(matrix_solver);
    mixed_branch(matrix_solver);
}

void FCI::alpha_branch(MatrixSolver *matrix_solver) {
    alpha_evaluation = new spin_evaluation*[dim_alpha];
    // Create the first spin string.
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<boost::dynamic_bitset<>> spin_string_alpha(0, this->addressing_scheme_alpha);  // spin string with address 0
    for (size_t Ia = 0; Ia < this->dim_alpha; Ia++) {  // Ib loops over all the addresses of the alpha spin strings
        alpha_evaluation[Ia] = new spin_evaluation[((K+1)-N_A)*N_A];
        size_t eval_index = 0;
        for (size_t p = 0; p < this->K; p++) {  // SO iteration 1
            int sign_p = 1;  // sign for the operator acting on the p-th SO
            if (spin_string_alpha.annihilate(p, sign_p)) {
                for (size_t q = 0; q < this->K; q++) {  // SO iteration 2
                    // alpha one-electron branch
                    int sign_q = sign_p;  // sign for the operator acting on the q-th SO
                    if (spin_string_alpha.create(q, sign_q)) {

                        // Retrieve address of the transformed SpinString
                        size_t Ja = spin_string_alpha.address(this->addressing_scheme_alpha);
                        // If alpha is major relative address in the total basis is multiplied by all beta combinations
                        alpha_evaluation[Ia][eval_index] = spin_evaluation{sign_q,p,q,Ja};
                        for (size_t Ib = 0; Ib < this->dim_beta; Ib++) {
                            matrix_solver->addToMatrix(sign_q * this->so_basis.get_h_SO(p, q), Ia * dim_beta + Ib, Ja * dim_beta + Ib);
                        }
                        eval_index++;
                        spin_string_alpha.annihilate(q);  // undo

                    }  // create q alpha

                    // alpha-alpha two-electron branch
                    sign_q = sign_p;  // sign for the operator acting on the q-th SO
                    if (spin_string_alpha.annihilate(q, sign_q)) {
                        for (size_t r = 0; r < this->K; r++) {
                            int sign_r = sign_q;  // sign for the operator acting on the r-th SO
                            if (spin_string_alpha.create(r, sign_r)) {
                                for (size_t s = 0; s < this->K; s++) {
                                    int sign_s = sign_r;  // sign for the operator acting on the s-th SO
                                    if (spin_string_alpha.create(s, sign_s)) {

                                        size_t Ja = spin_string_alpha.address(this->addressing_scheme_alpha);
                                        for (size_t Ib = 0; Ib < this->dim_beta; Ib++) {
                                            matrix_solver->addToMatrix(sign_s * this->so_basis.get_g_SO(p, s, q, r) / 2, Ia * dim_beta + Ib, Ja * dim_beta + Ib);
                                        }
                                        spin_string_alpha.annihilate(s);  // undo
                                    }  // create s alpha
                                }  // s
                                spin_string_alpha.annihilate(r);  // undo
                            }  // create r alpha
                        }  // r
                        spin_string_alpha.create(q);  // undo
                    }  // annihilate q alpha
                }
                spin_string_alpha.create(p);
            }  // annihilate p
        }
        if(Ia < dim_alpha-1){  // prevent last permutation to occur.
            spin_string_alpha.nextPermutation();
        }
    }
}

void FCI::beta_branch(MatrixSolver *matrix_solver) {
    beta_evaluation = new spin_evaluation*[dim_beta];
    // Create the first spin string.
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<boost::dynamic_bitset<>> spin_string_beta(0, this->addressing_scheme_beta);  // spin string with address 0

    for (size_t Ib = 0; Ib < this->dim_beta; Ib++) {
        beta_evaluation[Ib] = new spin_evaluation[((K+1)-N_B)*N_B];
        size_t eval_index = 0;
        for (size_t p = 0; p < this->K; p++) {  // SO iteration 1
            int sign_p = 1;
            // beta branch
            if (spin_string_beta.annihilate(p, sign_p)) {
                for (size_t q = 0; q < this->K; q++) {
                    // beta one-electron branch
                    int sign_q = sign_p;
                    if (spin_string_beta.create(q, sign_q)) {
                        size_t Jb = spin_string_beta.address(this->addressing_scheme_beta);

                        beta_evaluation[Ib][eval_index] = spin_evaluation{sign_q,p,q,Jb};

                        for (size_t Ia = 0; Ia < this->dim_alpha; Ia++) {
                            matrix_solver->addToMatrix(sign_q * this->so_basis.get_h_SO(p, q), Ia * dim_beta + Ib, Ia * dim_beta + Jb);
                        }
                        eval_index++;
                        spin_string_beta.annihilate(q);  // undo

                    }  // create q
                    // beta-beta two-electron branch
                    sign_q = sign_p;
                    if (spin_string_beta.annihilate(q, sign_q)) {
                        for (size_t r = 0; r < this->K; r++) {
                            int sign_r = sign_q;
                            if (spin_string_beta.create(r, sign_r)) {
                                for (size_t s = 0; s < this->K; s++) {
                                    int sign_s = sign_r;
                                    if (spin_string_beta.create(s, sign_s)) {

                                        size_t Jb = spin_string_beta.address(this->addressing_scheme_beta);
                                        for (size_t Ia = 0; Ia < this->dim_alpha; Ia++) {
                                            matrix_solver->addToMatrix(sign_s * this->so_basis.get_g_SO(p, s, q, r) / 2, Ia * dim_beta + Ib, Ia * dim_beta + Jb);
                                        }

                                        spin_string_beta.annihilate(s);  // undo
                                    }  // create s beta
                                }  // s
                                spin_string_beta.annihilate(r);  // undo
                            }  // create r beta
                        }  // r
                        spin_string_beta.create(q);  // undo
                    }  // annihilate q beta
                }  // q
                spin_string_beta.create(p);  // undo
            } // annihilate p beta
        }  // p
        if (Ib < dim_beta - 1) {  // prevent last permutation to occur.
            spin_string_beta.nextPermutation();
        }
    }



}

void FCI::mixed_branch(MatrixSolver *matrix_solver) {
    for (size_t Ia = 0; Ia < this->dim_alpha; Ia++) {
        for (size_t Ib = 0; Ib < this->dim_beta; Ib++) {
            for(size_t a = 0; a<((K+1)-N_A)*N_A; a++){
                spin_evaluation alpha = alpha_evaluation[Ia][a];
                for(size_t b = 0; b < ((K+1)-N_B)*N_B; b++){
                    spin_evaluation beta = beta_evaluation[Ib][b];
                    int sign = beta.sign*alpha.sign;
                    matrix_solver->addToMatrix(sign*this->so_basis.get_g_SO(alpha.p, alpha.q, beta.p, beta.q), Ia * dim_beta + Ib, alpha.address * dim_beta + beta.address);

                }

            }

        }
    }
}


/**
 *  @return the action of the FCI Hamiltonian of the coefficient vector @param x.
 */
Eigen::VectorXd FCI::matrixVectorProduct(const Eigen::VectorXd& x) {

}


/**
 *  @return the diagonal of the matrix representation of the FCI Hamiltonian.
 */
Eigen::VectorXd FCI::calculateDiagonal() {

}



/*
 *  CONSTRUCTORS
 */
/**
 *  Constructor based on a given @param so_basis and a number of alpha electron and beta electrons @param N_A and N_B respectively.
 */
FCI::FCI(libwint::SOBasis& so_basis, size_t N_A, size_t N_B) :
        BaseCI(so_basis, this->calculateDimension(so_basis.get_K(), N_A, N_B)),
        K (so_basis.get_K()),
        N_A (N_A),N_B (N_B),
        addressing_scheme_alpha (bmqc::AddressingScheme(this->K, this->N_A)),
        addressing_scheme_beta (bmqc::AddressingScheme(this->K, this->N_B)),
        dim_alpha (this->calculateDimension(so_basis.get_K(), N_A, 0)),
        dim_beta (this->calculateDimension(so_basis.get_K(), 0, N_B))
{
    // Do some input checks.
    if (this->K < this->N_A || this->K < this->N_B) {
        throw std::invalid_argument("Too many electrons of one spin to place into the given number of spatial orbitals.");
    }
}



/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  Given a number of spatial orbitals @param K and a number of alpha electrons @param N_A and beta electrons @param N_B, @return the dimension of
 *  the FCI space.
 */
size_t FCI::calculateDimension(size_t K, size_t N_A, size_t N_B) {

    // K N_A, N_B are expected to be small, so static-casting them to unsigned (what boost needs) is permitted.
    auto dim_double_alpha = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_A));
    auto dim_double_beta = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_B));
    auto dim_double_total = dim_double_alpha*dim_double_beta;

    // Check if the resulting dimension is appropriate to be stored in size_t
    return boost::numeric::converter<double, size_t>::convert(dim_double_total);
}


}  // namespace ci
