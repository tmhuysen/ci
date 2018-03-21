#define BOOST_TEST_MODULE "DOCI_test"



#include "DOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( dimension ) {

    BOOST_CHECK_EQUAL(ci::DOCI::calculateDimension(10, 1), 10);
    BOOST_CHECK_EQUAL(ci::DOCI::calculateDimension(6, 2), 15);
    BOOST_CHECK_EQUAL(ci::DOCI::calculateDimension(8, 3), 56);
}


BOOST_AUTO_TEST_CASE ( DOCI_beh_cation_klaas_dense ) {

    // Klaas' reference DOCI energy for BeH+ (obtained through Caitlin)
    double reference_doci_energy = -14.8782216937;


    // Do a DOCI calculation based on a given FCIDUMP file
    size_t K = 16;
    size_t N = 4;
    libwint::SOBasis so_basis ("../tests/reference_data/beh_cation_631g_caitlin.FCIDUMP", K);
    ci::DOCI doci (so_basis, N);

    bmqc::AddressingScheme addressing_scheme (K, N);
    doci.solve(ci::solver::SolverType::DENSE, addressing_scheme);


    // Calculate the total energy
    double internuclear_repulsion_energy = 1.5900757460937498e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-06);
}
