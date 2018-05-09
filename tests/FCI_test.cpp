#define BOOST_TEST_MODULE "FCI_test"



#include "FCI.hpp"
#include "RHF.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( FCI_dimension ) {

    BOOST_CHECK_EQUAL(ci::FCI::calculateDimension(10, 1, 1), 100);
    BOOST_CHECK_EQUAL(ci::FCI::calculateDimension(6, 2, 2), 225);
    BOOST_CHECK_EQUAL(ci::FCI::calculateDimension(8, 3, 3), 3136);
}


BOOST_AUTO_TEST_CASE ( FCI_constructor ) {


    // Prepare an SO basis to test the constructor
    libwint::Molecule water ("../tests/reference_data/h2o_Psi4_GAMESS.xyz");
    libwint::AOBasis ao_basis (water, "STO-3G");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (water, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());
    BOOST_CHECK_NO_THROW(ci::FCI (so_basis, 5, 5));  // N_alpha = 5, N_beta = 5
}
