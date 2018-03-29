#define BOOST_TEST_MODULE "DOCI_dense_test"



#include "DOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

/*

BOOST_AUTO_TEST_CASE ( DOCI_beh_cation_klaas_dense ) {

    // Klaas' reference DOCI energy for BeH+ (obtained through Caitlin)
    double reference_doci_energy = -14.8782216937;


    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOMullikenBasis so_basis ("../tests/reference_data/beh_cation_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons
    doci.solve(numopt::eigenproblem::SolverType::DENSE);


    // Calculate the total energy
    double internuclear_repulsion_energy = 1.5900757460937498e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_lih_klaas_dense ) {


    // Klaas' reference DOCI energy for LiH (obtained through Caitlin)
    double reference_doci_energy = -8.0029560313;


    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOMullikenBasis so_basis ("../tests/reference_data/lih_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons
    doci.solve(numopt::eigenproblem::SolverType::DENSE);


    // Calculate the total energy
    double internuclear_repulsion_energy = 9.6074293445896852e-01;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}
*/