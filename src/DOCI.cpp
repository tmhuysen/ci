#include "DOCI.hpp"



#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>



namespace ci {


/*
 *  OVERRIDDEN PRIVATE METHODS
 */

/**
 *  Given a @param matrix_solver, construct the DOCI Hamiltonian matrix in the solver's matrix representation. An
 *  @param addressing_scheme is used to speed up the searches for the addresses of coupling spin strings.
 */
void DOCI::constructHamiltonian(ci::solver::BaseMatrixSolver* matrix_solver, const bmqc::AddressingScheme& addressing_scheme) {



    std::cout << "I got to constructHamiltonian." << std::endl;

    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one.
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string (0, addressing_scheme);  // spin string with address 0

    for (size_t I = 0; I < this->dim; I++) {  // I loops over all the addresses of the spin strings
        if (I > 0) {
            spin_string.nextPermutation();
        }
        std::cout << "I is " << I << std::endl;

        for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            std::cout << "p is " << p << std::endl;
            if (spin_string.annihilate(p)) { // single excitation
                std::cout << "p is occupied in I" << std::endl;
                //A single excitation in doci can only be done in place.
                //Exciting only one electron to a vacant SO, will break the double occupancy(not part of the basis).
                double one_int = this->so_basis.get_h_SO(p,p);
                matrix_solver->addToMatrix(2 * one_int, I, I); //Twice : alpha and beta.
                //There are also two in-place double excitations of abba en baab combination.
                double two_int = this->so_basis.get_g_SO(p,p,p,p);
                matrix_solver->addToMatrix(two_int, I, I);

                for(size_t q = 0; q < p; q++){// q loops over SOs creation l=j is covered in the first loop and since we can't annihilate twice this combination would be redundant.
                    std::cout << "q is " << q << std::endl;
                    if (spin_string.create(q)){ //we can never excite a single electron to a new site we have to do it in pairs.
                        size_t address = spin_string.address(addressing_scheme);
                        std::cout << "q is unoccupied in I and the resulting address is: " << address << std::endl;
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

    std::cout << "I got out of constructHamiltonian." << std::endl;
}


/**
 *  Given a @param addressing_scheme, @return the action of the DOCI Hamiltonian of the coefficient vector @param x.
 */
Eigen::VectorXd DOCI::matrixVectorProduct(const bmqc::AddressingScheme& addressing_scheme, const Eigen::VectorXd& x) {

}


/**
 *  @return the diagonal of the matrix representation of the DOCI Hamiltonian.
 */
Eigen::VectorXd DOCI::calculateDiagonal() {

}



/*
 *  CONSTRUCTORS
 */
/**
 *  Constructor based on a given @param so_basis and a number of electrons @param N.
 */
DOCI::DOCI(libwint::SOBasis& so_basis, size_t N) :
    BaseCI(this->calculateDimension(so_basis.get_K(), N / 2), so_basis),
    N_P (N / 2),
    K (so_basis.get_K())

{

    // Do some input checks
    if ((N % 2) != 0) {
        throw std::invalid_argument("Your basis contains an odd amount of electrons and is not suitable for DOCI");
    }

    if (this->K < this->N_P) {
        throw std::invalid_argument("Too many electrons to place into the given number of spatial orbitals");
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

    std::cout << "I got to calculateDimension"
              << " with K=" << K
              << " and N_P=" << N_P << std::endl;

    // K and N_P are expected to be small, so static-casting them to unsigned (what boost needs) is permitted.
    auto dim_double = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_P));
    std::cout << "dim_double: " << dim_double << std::endl;

    // Check if the resulting dimension is appropriate to be stored in size_t
    return boost::numeric::converter<double, size_t>::convert(dim_double);
}


}  // namespace ci
