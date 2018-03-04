#define BOOST_TEST_MODULE "DOCI_test"

#include "DOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( DOCI_beh_klaas ) {

	// Klaas' reference DOCI energy for BeH+ (obtained through Caitlin)
	double reference_doci_energy = -14.8782216937;

	// Do a DOCI calculation based on a given FCIDUMP file
	const std::string beh_cation_631g_fcidump = "../tests/reference_data/beh_cation_631g_caitlin.FCIDUMP";
    doci::CI_basis ciBasis (beh_cation_631g_fcidump);
    doci::DOCI doci_test (&ciBasis);

	// Calculate the total energy as the sum of the lowest energy eigenstate + the internuclear repulsion
	const doci::State& ground = doci_test.get_lowest_eigenstate();
    double test_doci_energy = ground.getEval() + ciBasis.getInternuclear_repulsion();

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( DOCI_lih_klaas ) {

	// Klaas' reference DOCI energy for LiH (obtained through Caitlin)
	double reference_doci_energy = -8.0029560313;

	// Do a DOCI calculation based on a given FCIDUMP file
	const std::string lih_631g_fcidump = "../tests/reference_data/lih_631g_caitlin.FCIDUMP";
    doci::CI_basis ciBasis (lih_631g_fcidump);
    doci::DOCI doci (&ciBasis);

	// Calculate the total energy as the sum of the lowest energy eigenstate + the internuclear repulsion
    const doci::State &ground = doci.get_lowest_eigenstate();
    double en = ground.getEval() + ciBasis.getInternuclear_repulsion();

    BOOST_CHECK(std::abs(en - (reference_doci_energy)) < 1.0e-06);
}
