#define BOOST_TEST_MODULE "RDM_test"



#include "DOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( compute_test ) {

    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../tests/reference_data/beh_cation_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons
    doci.solve(numopt::eigenproblem::SolverType::DENSE);
    doci.compute1RDMs();
    doci.compute2RDMs();

    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE ( reference_tests ) {

    // TODO: add reference
    BOOST_CHECK(false);
}
