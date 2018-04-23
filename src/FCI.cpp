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
    auto start = std::chrono::high_resolution_clock::now();
    // Create the first spin string.
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string_alpha (0, this->addressing_scheme_alpha);  // spin string with address 0


    for (size_t Ia = 0; Ia < this->dim_alpha; Ia++) {  // Ib loops over all the addresses of the alpha spin strings
        bmqc::SpinString<boost::dynamic_bitset<>> spin_string_beta(0,
                                                         this->addressing_scheme_beta);  // spin string with address 0
        for (size_t Ib = 0; Ib < this->dim_beta; Ib++) {
            for (size_t p = 0; p < this->K; p++) {  // SO iteration 1

                // alpha branch
                int sign_p = 1;  // sign for the operator acting on the p-th SO
                if (spin_string_alpha.annihilate(p, sign_p)) {
                    for (size_t q = 0; q < this->K; q++) {  // SO iteration 2

                        // alpha one-electron branch
                        int sign_q = sign_p;  // sign for the operator acting on the q-th SO
                        if (spin_string_alpha.create(q, sign_q)) {
                            // Retrieve address of the transformed SpinString
                            size_t Ja = spin_string_alpha.address(this->addressing_scheme_alpha);
                            // If alpha is major relative address in the total basis is multiplied by all beta combinations
                            matrix_solver->addToMatrix(sign_q*this->so_basis.get_h_SO(p, q),
                                                       Ia * dim_beta + Ib, Ja * dim_beta + Ib);
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

                                            matrix_solver->addToMatrix(sign_s*this->so_basis.get_g_SO(p, s, q, r)/2,
                                                                       Ia * dim_beta + Ib, Ja * dim_beta + Ib);
                                            spin_string_alpha.annihilate(s);  // undo
                                        }  // create s alpha
                                    }  // s
                                    spin_string_alpha.annihilate(r);  // undo
                                }  // create r alpha
                            }  // r
                            spin_string_alpha.create(q);  // undo
                        }  // annihilate q alpha

                        // alpha-beta two-electron branch
                        sign_q = sign_p;
                        if (spin_string_beta.annihilate(q, sign_q)) {
                            for (size_t r = 0; r < this->K; r++) {
                                int sign_r = sign_q;
                                if (spin_string_beta.create(r, sign_r)) {
                                    for (size_t s = 0; s < this->K; s++) {
                                        int sign_s = sign_r;
                                        if (spin_string_alpha.create(s, sign_s)) {
                                            // Retrieve address of the transformed SpinString alpha
                                            size_t Ja = spin_string_alpha.address(this->addressing_scheme_alpha);
                                            // Retrieve address of the transformed SpinString beta
                                            size_t Jb = spin_string_beta.address(this->addressing_scheme_beta);

                                            matrix_solver->addToMatrix(sign_s*this->so_basis.get_g_SO(p, s, q, r),
                                                                       Ia * dim_beta + Ib, Ja * dim_beta + Jb);  // We do not divide by 2 because we will not evaluate this branch in the beta branch.
                                            spin_string_alpha.annihilate(s);  // undo
                                        }  // create s alpha
                                    }  // s
                                    spin_string_beta.annihilate(r);  // undo
                                }  // annihilate r beta
                            }  // r
                            spin_string_beta.create(q);  // undo
                        }  // annihilate q beta
                    }  // q
                    spin_string_alpha.create(p);  // undo
                }  // annihilate p alpha
                sign_p = 1;
                // beta branch 
                if (spin_string_beta.annihilate(p, sign_p)) {
                    for (size_t q = 0; q < this->K; q++) {
                        // beta one-electron branch
                        int sign_q = sign_p;
                        if (spin_string_beta.create(q, sign_q)) {
                            size_t Jb = spin_string_beta.address(this->addressing_scheme_beta);
                            matrix_solver->addToMatrix(sign_q*this->so_basis.get_h_SO(p, q),
                                                       Ia * dim_beta + Ib, Ia * dim_beta + Jb);
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

                                            matrix_solver->addToMatrix(sign_s*this->so_basis.get_g_SO(p, s, q, r)/2,
                                                                       Ia * dim_beta + Ib, Ia * dim_beta + Jb);
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
            if(Ib < dim_beta-1){  // prevent last permutation to occur.
                spin_string_beta.nextPermutation();
            }
        }
        if(Ia < dim_alpha-1){  // prevent last permutation to occur.
            spin_string_alpha.nextPermutation();
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();

    // Process the chrono time and output
    auto elapsed_time = stop - start;           // in nanoseconds
    auto seconds = elapsed_time.count() / 1e9;  // in seconds
    std::cout << dim << " : " << seconds << std::endl;
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


/**
 *  Computes all of the one reduced density matrix.
 */
void FCI::compute1RDM(){
    //initialize as zero matrix.
    this->one_rdm_aa = Eigen::MatrixXd::Zero(this->K,this->K);
    this->one_rdm_bb = Eigen::MatrixXd::Zero(this->K,this->K);

    bmqc::SpinString<boost::dynamic_bitset<>> spin_string_alpha (0, this->addressing_scheme_alpha);  // spin string with address 0
    for (size_t I = 0; I < this->dim_alpha; I++) {  // I loops over all the addresses of the spin strings
        if (I > 0) {
            spin_string_alpha.nextPermutation();
        }for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            int sign_p = 1;
            if(spin_string_alpha.annihilate(p,sign_p)){
                double coefficient =  0;
                for(size_t i = 0;i<dim_beta;i++){
                    coefficient += this->eigensolver_ptr->get_eigenvector()(I*dim_beta+i)*this->eigensolver_ptr->get_eigenvector()(I*dim_beta+i);
                }  // dim_beta loop
                this->one_rdm_aa(p,p) += coefficient;
                for(size_t q = 0; q < p;q++){
                    int sign_q = sign_p;
                    if(spin_string_alpha.create(q,sign_q)){
                        size_t J = spin_string_alpha.address(this->addressing_scheme_alpha);

                        double coefficient = 0;
                        for(size_t i = 0;i<dim_beta;i++){
                            coefficient += this->eigensolver_ptr->get_eigenvector()(I*dim_beta+i)
                                           *this->eigensolver_ptr->get_eigenvector()(J*dim_beta+i);
                        }  // dim_beta loop
                        one_rdm_aa(q,p) += sign_q*coefficient;
                        one_rdm_aa(p,q) += sign_q*coefficient;

                        spin_string_alpha.annihilate(q);

                    }  // create q
                }  // q
                spin_string_alpha.create(p);
            }  // anni p
        }  // p
    }  // I
    bmqc::SpinString<boost::dynamic_bitset<>> spin_string_beta (0, this->addressing_scheme_beta);  // spin string with address 0
    for (size_t I = 0; I < this->dim_beta; I++) {  // I loops over all the addresses of the spin strings
        if (I > 0) {
            spin_string_beta.nextPermutation();
        }for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            int sign_p = 1;
            if(spin_string_beta.annihilate(p,sign_p)){

                double coefficient =  0;
                for(size_t i = 0;i<dim_alpha;i++){
                    coefficient += this->eigensolver_ptr->get_eigenvector()(i*dim_beta+I)*this->eigensolver_ptr->get_eigenvector()(i*dim_beta+I);
                }  // dim_alpha
                this->one_rdm_bb(p,p) += coefficient;
                for(size_t q = 0; q < p;q++){
                    int sign_q = sign_p;
                    if(spin_string_beta.create(q,sign_q)){
                        size_t J = spin_string_beta.address(this->addressing_scheme_beta);

                        double coefficient = 0;
                        for(size_t i = 0;i<dim_alpha;i++){
                            coefficient += this->eigensolver_ptr->get_eigenvector()(i*dim_beta+I)
                                           *this->eigensolver_ptr->get_eigenvector()(i*dim_beta+J);
                        }
                        one_rdm_bb(q,p) += sign_q*coefficient;
                        one_rdm_bb(p,q) += sign_q*coefficient;
                        spin_string_beta.annihilate(q);

                    }  // create q


                }  // q
                spin_string_beta.create(p);
            }  // anni p
        }  // p
    }  // I
    this->are_computed_one_rdm = true;
}


/**
 *  Computes all of the two reduced density matrix.
 */
void FCI::compute2RDM(){

}



/*
 *  CONSTRUCTORS
 */
/**
 *  Constructor based on a given @param so_basis and a number of alpha electron and beta electrons @param N_A and N_B respectively.
 */
FCI::FCI(libwint::SOMullikenBasis& so_basis, size_t N_A, size_t N_B) :
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
