#define BOOST_TEST_MODULE "State"

#include "State.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( isDegenerate ) {

    double eigenvalue_1 = 100.00101;
    double eigenvalue_2 = 100.001010000000001;
    Eigen::VectorXd eigenvector = Eigen::VectorXd::Zero(2);  // an unimportant vector for this test

    // We expect the following States to have equal eigenvalues
    ci::State state_1 (eigenvalue_1, eigenvector);
    ci::State state_2 (eigenvalue_2, eigenvector);
    BOOST_CHECK(state_1.isDegenerate(state_2));
}


BOOST_AUTO_TEST_CASE ( hasLowerEnergy ) {

    double eigenvalue_1 = 100.0;
    double eigenvalue_2 = 100.5;
    Eigen::VectorXd eigenvector = Eigen::VectorXd::Zero(2);  // an unimportant vector for this test

    // We expect the state_1 to have a lower energy than state_2
    ci::State state_1 (eigenvalue_1, eigenvector);
    ci::State state_2 (eigenvalue_2, eigenvector);
    BOOST_CHECK(state_1.hasLowerEnergy(state_2));
}
