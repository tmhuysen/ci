#ifndef CI_FCI_HPP
#define CI_FCI_HPP


#include <bmqc.hpp>
#include <boost/math/special_functions.hpp>
#include <iostream>

#include "CI_Abstract_Class.hpp"


namespace doci {

class FCI: public doci::CI {
private:
	size_t nelec_a; // number of alpha electrons
	size_t nelec_b; // number of beta electrons
	size_t nbf_a; // number of alpha basis functions
	size_t nbf_b; // number of beta basis functions
	bmqc::AddressingScheme ad_mat_a; // AddressingScheme for alpha electrons
	bmqc::AddressingScheme ad_mat_b; // AddressingScheme for beta electrons

	/**
	* calculate hamiltonian elements.
	* @param start,end : indicates the bf you want to start with and where the iteration ends(excluded).
	*/
	void calculateCI(size_t start, size_t end) override ;

	/**
	 * Helper function for the constructors
	 */
	void construct() override;


public:
	FCI(CI_basis *ciBasis);

};

} // namespace doci


#endif // CI_FCI_HPP
