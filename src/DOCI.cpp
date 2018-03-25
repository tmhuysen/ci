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

    Eigen::VectorXd matvec = Eigen::VectorXd::Zero(this->dim);

    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one.
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string (0, this->addressing_scheme);  // spin string with address 0


    for (size_t I = 0; I < this->dim; I++) {  // I loops over all the addresses of the spin strings

        for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            if (spin_string.annihilate(p)) {  // if p in I
                matvec(I) += 2 * this->so_basis.get_h_SO(p,p) * x(I);
                matvec(I) += this->so_basis.get_g_SO(p,p,p,p) * x(I);

                for (size_t q = 0; q < this->K; q++) {  // q loops over SOs
                    if (spin_string.create(q)) {  // if q not in I, so q!=p automatically
                        size_t J = spin_string.address(this->addressing_scheme);  // J is the address of a string that couples to I

                        matvec(I) += this->so_basis.get_g_SO(p,q,p,q) * x(J);

                        spin_string.annihilate(q);  // reset the spin string after the previous creation
                    }

                    else {  // if q in I
                        if (p != q) {
                            matvec(I) += (2*this->so_basis.get_g_SO(p,p,q,q) - this->so_basis.get_g_SO(p,q,q,p)) * x(I);
                        }
                    }
                }  // q

                spin_string.create(p);  // reset the spin string after the previous annihilation
            }
        }  // p

        spin_string.nextPermutation();
    }  // address (I) loop

    return matvec;
}


/**
 *  @return the diagonal of the matrix representation of the DOCI Hamiltonian.
 */
Eigen::VectorXd DOCI::calculateDiagonal() {

    Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(this->dim);

    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string (0, this->addressing_scheme);

    for (size_t I = 0; I < this->dim; I++) {  // I loops over addresses of spin strings

        for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            if (spin_string.isOccupied(p)) {  // p is in I
                diagonal(I) += 2 * this->so_basis.get_h_SO(p,p) + this->so_basis.get_g_SO(p,p,p,p);

                for (size_t q = 0; q < p; q++) {  // q loops over SOs
                    if (spin_string.isOccupied(q)) {  // q is in I

                        // Since we are doing a restricted summation q<p, we should multiply by 2 since the summand argument is symmetric.
                        diagonal(I) += 2 * (2*this->so_basis.get_g_SO(p,p,q,q) - this->so_basis.get_g_SO(p,q,q,p));
                    }
                }  // q loop
            }
        }  // p loop

        spin_string.nextPermutation();
    }  // address (I) loop

    return diagonal;
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

    std::cout << "DOCI instance was created." << std::endl;
}


/**
 *  Constructor based on a given @param so_basis and a @param molecule.
 */
DOCI::DOCI(libwint::SOBasis& so_basis, const libwint::Molecule& molecule) :
    DOCI (so_basis, molecule.get_N())
{}



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
