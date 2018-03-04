#ifndef CI_CI_ABSTRACT_CLASS_HPP
#define CI_CI_ABSTRACT_CLASS_HPP

#include "State.hpp"
#include "Hamiltonian.hpp"
#include "CI_basis.hpp"


namespace doci {
class CI {

protected:
	size_t nelec; // number of electrons
	size_t K; // number of spatial orbitals
	size_t nbf; // number of basis functions
	doci::State lowestEigenState;
	doci::CI_basis* basis; //contains all information required to do a CI calculation
	doci::Hamiltonian* hamiltonian; //abstract hamiltonian

protected:

	/**
	* calculate hamiltonian elements.
	* @param start,end : indicates the bf you want to start with and where the iteration ends(excluded).
	*/
	virtual void calculateCI(size_t start, size_t end)=0;
	/**
	 * Helper function for the constructors
	 */
	virtual void construct()=0;

public:
	/** Constructor based on a given CI_basis
	 * Initializes the base variables of the class, and calculates the number of basis functions (nbf)
	 * Hamiltonian is picked by the program based on the nbf
	 */
	CI(CI_basis *ciBasis);

	/** Default Constructor
	 */
	CI();

    /** Destructor
     * clears the hamiltonian pointer.
     */
    virtual ~CI();

	/** Getters
	 */
	const State &getLowestEigenState() const;
	Hamiltonian *getHamiltonian() const;
};

}  // namespace doci


#endif // CI_CI_ABSTRACT_CLASS_HPP
