#define BOOST_TEST_MODULE "DOCI_test"



#include "DOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( DOCI_dimension ) {

    BOOST_CHECK_EQUAL(ci::DOCI::calculateDimension(10, 1), 10);
    BOOST_CHECK_EQUAL(ci::DOCI::calculateDimension(6, 2), 15);
    BOOST_CHECK_EQUAL(ci::DOCI::calculateDimension(8, 3), 56);
}


//BOOST_AUTO_TEST_CASE ( DOCI_constructor ) {
//
//    libwint::SOMullikenBasis so_basis ("../tests/reference_data/beh_cation_631g_caitlin.FCIDUMP", 16);  // 16 SOs
//    BOOST_CHECK_NO_THROW(ci::DOCI (so_basis, 4));  // non-faulty constructor
//}
