#ifndef CI_BASECI_ABSTRACT_CLASS_HPP
#define CI_BASECI_ABSTRACT_CLASS_HPP


#include <libwint.hpp>

#include "State.hpp"
#include "BaseHamiltonianSolver.hpp"




namespace ci {


class BaseCI {
protected:
    const size_t K;  // number of spatial orbitals
    const size_t dim;  // the dimension of the Fock space related to this problem

    libwint::SOBasis& so_basis;

    // We also need a protected member that refers to the chosen HamiltonianSolver.
    // We can't do this with a reference (since an abstract class cannot be instantiated) so we must use a pointer.
    // This pointer can point to any instance of any derived HamiltonianSolver.
    std::unique_ptr<ci::BaseHamiltonianSolver> hamiltonian_solver_ptr;


protected:
    /**
     * calculate hamiltonian elements.
     * @param start,end : indicates the bf you want to start with and where the iteration ends(excluded).
     */
    virtual void calculateCI(size_t start, size_t end) = 0;
//    /**
//     * Helper function for the constructors
//     */
//    virtual void construct()=0;


public:
    /** Constructor based on a given CI_basis
     * Initializes the base variables of the class, and calculates the number of basis functions (nbf)
     * Hamiltonian is picked by the program based on the nbf
     */
    CI(CI_basis *ciBasis);

    /**
     * Default constructor
     */
    CI();

    /**
     * Getters
     */
    const State &getLowestEigenState() const;
    Hamiltonian *getHamiltonian() const;
};


}  // namespace ci


#endif // CI_BASECI_ABSTRACT_CLASS_HPP
