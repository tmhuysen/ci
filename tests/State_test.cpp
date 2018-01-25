#define BOOST_TEST_MODULE "State_tests"

#include "State.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( State_test) {
    doci::State test_1 (100.00101, Eigen::VectorXd());
    doci::State test_2 (100.001010000000001, Eigen::VectorXd());
    BOOST_CHECK(test_1==test_2);
    doci::State test_3 (200, Eigen::VectorXd());
    BOOST_CHECK((test_3>test_1));

}


