#ifndef CI_DOCI_HPP
#define CI_DOCI_HPP

#include "CI_Abstract_Class.hpp"

#include <bmqc.hpp>
#include <iostream>


namespace doci {

class DOCI : public CI {
private:
	size_t npairs; // number of electron pairs
	bmqc::AddressingScheme ad_mat;

	/**
	* calculate hamiltonian elements.
	* @param start,end : indicates the bf you want to start with and where the iteration ends(excluded).
	*/
	void calculateCI(size_t start, size_t end) override;

	/**
	 * Helper function for the constructors
	 */
	void construct() override;

public:
	DOCI(CI_basis *ciBasis);
	DOCI(CI_basis_mulliken *ciBasis, double constraint, std::vector<size_t> set_of_AO);
};

}  // namespace doci

#endif // CI_DOCI_HPP
