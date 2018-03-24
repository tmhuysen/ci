#include "FCI.hpp"



#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>



namespace ci {


/*
 *  OVERRIDDEN PRIVATE METHODS
 */

/**
 *  Given a @param matrix_solver, construct the DOCI Hamiltonian matrix in the solver's matrix representation.
 */
void FCI::constructHamiltonian(numopt::eigenproblem::BaseMatrixSolver* matrix_solver) {
    /*
    // Create the first spin string.
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string_alpha (0, this->addressing_scheme_alpha);  // spin string with address 0


    for (size_t I = 0; I < this->dim; I++) {  // I loops over all the addresses of the alpha spin strings

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
    */
}


/**
 *  @return the action of the DOCI Hamiltonian of the coefficient vector @param x.
 */
Eigen::VectorXd FCI::matrixVectorProduct(const Eigen::VectorXd& x) {

}


/**
 *  @return the diagonal of the matrix representation of the DOCI Hamiltonian.
 */
Eigen::VectorXd FCI::calculateDiagonal() {

}



/*
 *  CONSTRUCTORS
 */
/**
 *  Constructor based on a given @param so_basis and a number of electrons @param N.
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
 *  Given a number of spatial orbitals @param K and a number of electron pairs @param N_P, @return the dimension of
 *  the DOCI space.
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
