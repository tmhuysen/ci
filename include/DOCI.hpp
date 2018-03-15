#ifndef CI_DOCI_HPP
#define CI_DOCI_HPP

#include "CI_Abstract_Class.hpp"
#include "bitlong.hpp"
#include <bmqc.hpp>
#include <iostream>


namespace doci {

class DOCI : public CI {
public:
    size_t npairs; // number of electron pairs
    bmqc::AddressingScheme ad_mat;

    /** Calculates the diagonal elements
     */
    void calculateDiagonal();

    /** Calculates the off diagonal elements
     */
    void calculateOffDiagonal();

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
    DOCI(CI_basis *ciBasis, StorageType type);
};

}  // namespace doci

#endif // CI_DOCI_HPP
