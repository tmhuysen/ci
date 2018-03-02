#define BOOST_TEST_MODULE "CI_basis_test"

#include "CI_basis.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( constructor_rhf ) {

    // Do an RHF calculation
    const std::string xyzfilename = "../tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule water (xyzfilename);
    libwint::Basis basis (water, basis_name);
    hf::rhf::RHF rhf (basis, threshold);

    // Test that the constructor doesn't raise an error when a correct input RHF instance is supplied
    BOOST_REQUIRE_NO_THROW(doci::CI_basis ciBasis (rhf));
}


BOOST_AUTO_TEST_CASE ( constructor_filename ) {

    // Test that the constructor raises an error when a file can't be read
    const std::string doci_faulty_path = "this_is_a_non_existent_path";
    BOOST_REQUIRE_THROW(doci::CI_basis ciBasis (doci_faulty_path), std::runtime_error);


    // Test that the constructor doesn't raise an error when a correct FCIDUMP file is provided
    const std::string co_631g_fcidump = "../tests/reference_data/co_631g.FCIDUMP";
    BOOST_REQUIRE_NO_THROW(doci::CI_basis ciBasis (co_631g_fcidump));
}