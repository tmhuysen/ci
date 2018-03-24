#include "DOCI.hpp"



#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>



namespace ci {


/*
 *  OVERRIDDEN PRIVATE METHODS
 */

/**
 *  Given a @param matrix_solver, construct the DOCI Hamiltonian matrix in the solver's matrix representation.
 */
void DOCI::constructHamiltonian(numopt::eigenproblem::BaseMatrixSolver* matrix_solver) {

    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one.
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string (0, this->addressing_scheme);  // spin string with address 0


    for (size_t I = 0; I < this->dim; I++) {  // I loops over all the addresses of the spin strings

        double diagonal_element = 0.0;
        for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            if (spin_string.annihilate(p)) {  // if p in I
                diagonal_element += 2 * this->so_basis.get_h_SO(p,p) + this->so_basis.get_g_SO(p,p,p,p);

                for (size_t q = 0; q < p; q++) {  // q loops over SOs
                    if (spin_string.create(q)) {  // if q not in I
                        size_t J = spin_string.address(this->addressing_scheme);  // J is the address of a string that couples to I

                        // The loops are p->K and q<p. So, we should normally multiply by a factor 2 (since the summand is symmetric).
                        // However, we are setting both of the symmetric indices of Hamiltonian, so no factor 2 is required.
                        matrix_solver->addToMatrix(this->so_basis.get_g_SO(p,q,p,q), I, J);
                        matrix_solver->addToMatrix(this->so_basis.get_g_SO(p,q,p,q), J, I);

                        spin_string.annihilate(q);  // reset the spin string after the previous creation
                    }

                    else {  // if q in I
                        diagonal_element += 2 * (2*this->so_basis.get_g_SO(p,p,q,q) - this->so_basis.get_g_SO(p,q,q,p));
                    }
                }  // q

                spin_string.create(p);  // reset the spin string after the previous annihilation
            }
        }  // p

        matrix_solver->addToMatrix(diagonal_element, I, I);
        spin_string.nextPermutation();
    }  // address (I) loop
}


/**
 *  @return the action of the DOCI Hamiltonian of the coefficient vector @param x.
 */
Eigen::VectorXd DOCI::matrixVectorProduct(const Eigen::VectorXd& x) {

}


/**
 *  @return the diagonal of the matrix representation of the DOCI Hamiltonian.
 */
Eigen::VectorXd DOCI::calculateDiagonal() {

}


/**
 *  Computes all of the one reduced density matrix.
 */
void DOCI::compute1RDM(){
    //initialize as zero matrix.
    this->one_rdm_aa = Eigen::MatrixXd::Zero(this->K,this->K);
    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one.
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string (0, this->addressing_scheme);  // spin string with address 0

    for (size_t I = 0; I < this->dim; I++) {  // I loops over all the addresses of the spin strings
        if (I > 0) {
            spin_string.nextPermutation();
        }for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            if(spin_string.isOccupied(p)){
                double coefficient = this->eigensolver_ptr->get_eigenvector()(I)*this->eigensolver_ptr->get_eigenvector()(I);
                this->one_rdm_aa(p,p) += coefficient;
            }
        }
    }
    this->one_rdm_bb = this->one_rdm_aa;  // For DOCI these are equal
    one_rdms_computed = true;
}


/**
 *  Computes all of the two reduced density matrix.
 */
void DOCI::compute2RDM(){
    //initialize as zero tensors.
    this->two_rdm_aaaa = Eigen::Tensor<double, 4>(this->K,this->K,this->K,this->K);
    this->two_rdm_aaaa.setZero();
    this->two_rdm_abba = Eigen::Tensor<double, 4>(this->K,this->K,this->K,this->K);
    this->two_rdm_abba.setZero();
    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one.
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string (0, this->addressing_scheme);  // spin string with address 0

    for (size_t I = 0; I < this->dim; I++) {  // I loops over all the addresses of the spin strings
        if (I > 0) {
            spin_string.nextPermutation();
        }for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            if(spin_string.annihilate(p)){

                double coefficient = this->eigensolver_ptr->get_eigenvector()(I)*this->eigensolver_ptr->get_eigenvector()(I);
                this->two_rdm_abba(p,p,p,p) += coefficient;

                for(size_t q = 0;q<p; q++){  // q loops over SOs with a smaller index than p
                    if(spin_string.create(q)){
                        size_t address = spin_string.address(this->addressing_scheme);
                        double coefficient = this->eigensolver_ptr->get_eigenvector()(address)*this->eigensolver_ptr->get_eigenvector()(I);
                        // Physical notation
                        this->two_rdm_abba(q,q,p,p) += coefficient;
                        // Symmetry
                        this->two_rdm_abba(p,p,q,q) += coefficient;

                        spin_string.annihilate(q);
                    }else {  // if we can't create we can annihilated
                        // Physical notation
                        // In-place and exchange excitation for 2 alpha electrons.
                        this->two_rdm_aaaa(p, q, p, q) += coefficient;
                        // Symmetry
                        this->two_rdm_aaaa(q, p, q, p) += coefficient;
                        // Anti-symmetry (exchange)
                        this->two_rdm_aaaa(p, q, q, p) -= coefficient;
                        this->two_rdm_aaaa(q, p, p, q) -= coefficient;

                        // In-place excitation for alpha beta pair (abba)
                        this->two_rdm_abba(p, q, p, q) += coefficient;
                        this->two_rdm_abba(q, p, q, p) += coefficient;
                    }
                }
                spin_string.create(p);
            }
        }
    }
    this->two_rdm_bbbb = this->two_rdm_aaaa;  // For DOCI these are equal
    this->two_rdm_baab = this->two_rdm_abba;  // For CI these are equal
    two_rdms_computed = true;
}



/*
 *  CONSTRUCTORS
 */
/**
 *  Constructor based on a given @param so_basis and a number of electrons @param N.
 */
DOCI::DOCI(libwint::SOBasis& so_basis, size_t N) :
    BaseCI(so_basis, this->calculateDimension(so_basis.get_K(), N / 2)),
    K (so_basis.get_K()),
    N_P (N / 2),
    addressing_scheme (bmqc::AddressingScheme(this->K, this->N_P)) // since in DOCI, alpha==beta, we should make an
                                                                   // addressing scheme with the number of PAIRS.
{
    // Do some input checks.
    if ((N % 2) != 0) {
        throw std::invalid_argument("You gave an odd amount of electrons, which is not suitable for DOCI.");
    }

    if (this->K < this->N_P) {
        throw std::invalid_argument("Too many electrons to place into the given number of spatial orbitals.");
    }
}



/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  Given a number of spatial orbitals @param K and a number of electron pairs @param N_P, @return the dimension of
 *  the DOCI space.
 */
size_t DOCI::calculateDimension(size_t K, size_t N_P) {

    // K and N_P are expected to be small, so static-casting them to unsigned (what boost needs) is permitted.
    auto dim_double = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_P));

    // Check if the resulting dimension is appropriate to be stored in size_t
    return boost::numeric::converter<double, size_t>::convert(dim_double);
}


}  // namespace ci
