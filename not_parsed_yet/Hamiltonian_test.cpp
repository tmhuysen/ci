#define BOOST_TEST_MODULE "Hamiltonian_test"
#include <Hamiltonian.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( Hamiltonian ) {
    //make a 3by3
    doci::Hamiltonian *hamPtr = doci::Hamiltonian::make_hamiltonian(3);
    hamPtr->add(1,0,0); //fill diagonal
    hamPtr->add(2,1,1);
    hamPtr->add(3,2,2);
    hamPtr->solve(); //"solve" the eigenvalue problem
    doci::State ground = hamPtr->getGroundstates().at(0);
    BOOST_CHECK(std::abs(1- ground.get_eigenvalue()) < 1.0e-04);
}
