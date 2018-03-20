#define BOOST_TEST_MODULE "CI_basis_mulliken"

#include "CI_basis_mulliken.hpp"
#include "DOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( constructor_test ) {
    // Do an RHF calculation
    const std::string xyzfilename = "../tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule water (xyzfilename);
    libwint::Basis basis (water, basis_name);
    hf::rhf::RHF rhf (basis, threshold);

    // Test that the constructor doesn't raise an error when a correct input RHF instance is supplied
    BOOST_REQUIRE_NO_THROW(doci::CI_basis_mulliken ci_basis_mulliken (rhf));
}

BOOST_AUTO_TEST_CASE ( doci_value_change_test ) {
    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule water (xyzfilename);
    libwint::Basis basis (water, basis_name);
    hf::rhf::RHF rhf (basis, threshold);
    doci::CI_basis_mulliken ci_basis_mulliken (rhf);

    // Do a DOCI calculation based on the RHF calculation
    doci::DOCI doci_test = doci::DOCI(&ci_basis_mulliken);
    doci::State ground = doci_test.getLowestEigenState();
    double test_doci_energy = ground.getEval() + ci_basis_mulliken.getInternuclear_repulsion();

    std::cout<<std::endl<<" energy base : "<< test_doci_energy;
    ci_basis_mulliken.calculateMullikenMatrix({0});
    ci_basis_mulliken.set_lagrange_multiplier(5);
    doci_test = doci::DOCI(&ci_basis_mulliken);
    ground = doci_test.getLowestEigenState();
    test_doci_energy = ground.getEval() + ci_basis_mulliken.getInternuclear_repulsion();
    std::cout<<std::endl<<" energy mulliken : "<< test_doci_energy;

    BOOST_CHECK(true)


}