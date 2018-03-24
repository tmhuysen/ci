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
        if (I > 0) {
            spin_string.nextPermutation();
        }
        for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            if (spin_string.annihilate(p)) { // single excitation
                //A single excitation in doci can only be done in place.
                //Exciting only one electron to a vacant SO, will break the double occupancy(not part of the basis).
                double one_int = this->so_basis.get_h_SO(p,p);
                matrix_solver->addToMatrix(2 * one_int, I, I); //Twice : alpha and beta.
                //There are also two in-place double excitations of abba en baab combination.
                double two_int = this->so_basis.get_g_SO(p,p,p,p);
                matrix_solver->addToMatrix(two_int, I, I);

                for(size_t q = 0; q < p; q++){// q loops over SOs creation l=j is covered in the first loop and since we can't annihilate twice this combination would be redundant.
                    if (spin_string.create(q)){ //we can never excite a single electron to a new site we have to do it in pairs.
                        size_t address = spin_string.address(this->addressing_scheme);
                        //integrals parameters are entered in chemical notation!
                        // This means that first 2 parameters are for the first electrons and subsequent ones are for the second
                        //Multiply by 2 getting rid of 1/2 two electron term because we have 2 equal combinations:
                        //abba and baab. We do not correct for the truncated SO iteration because we only calculate the lower triagonal.
                        //and then copy accordingly
                        double mix_spin_two_int = this->so_basis.get_g_SO(p,q,p,q);
                        matrix_solver->addToMatrix(mix_spin_two_int, I, address);
                        matrix_solver->addToMatrix(mix_spin_two_int, address, I);

                        spin_string.annihilate(q);//flip back (so we don't need to copy the set)
                    }else{ //if we can't create we can annihilated (but only on the diagonal, no transform required

                        // Integral parameters are entered in chemical notation!
                        double same_spin_two_int = this->so_basis.get_g_SO(p,p,q,q);//=mixed_spin_two_int (exciting a beta and an alpha in-place)
                        double same_spin_two_int_negative = -this->so_basis.get_g_SO(p,q,q,p);//mixed_spin does not have this because it would result in 0 term (integral of alpha-beta)

                        //We don't iterate over all the SO's the second time so multiply by 2 getting rid of 1/2 two electron term.
                        //multiply by 2 again because alpha,alpha is the same as beta,beta combinations.
                        //same_spin (positive) = mixed, so multiply that by 2 again.
                        matrix_solver->addToMatrix(4*same_spin_two_int + 2*same_spin_two_int_negative, I, I);
                    }

                }
                spin_string.create(p);
            }
        }
    }
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
