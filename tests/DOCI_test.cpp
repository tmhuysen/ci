#define BOOST_TEST_MODULE "DOCI_test"



#include "DOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( dimension ) {

    BOOST_CHECK_EQUAL(ci::DOCI::calculateDimension(10, 1), 10);
    BOOST_CHECK_EQUAL(ci::DOCI::calculateDimension(6, 2), 15);
    BOOST_CHECK_EQUAL(ci::DOCI::calculateDimension(8, 3), 56);
}
