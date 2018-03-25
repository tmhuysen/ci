#define BOOST_TEST_MODULE "DOCI_Davidson_test"


#include <hf.hpp>

#include "DOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( DOCI_h2_sto3g_Davidson ) {

    // Prepare an SOBasis from an RHF calculation
    libwint::Molecule h2 ("../tests/reference_data/h2.xyz");
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();

    libwint::AOBasis ao_basis (h2, "STO-3G");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-08);
    rhf.solve();

    Eigen::MatrixXd coefficient_matrix = rhf.get_C_canonical();
    libwint::SOBasis so_basis (ao_basis, coefficient_matrix);


    // Calculate the DOCI energy using full diagonalization of the dense Hamiltonian
    ci::DOCI doci_dense (so_basis, h2);
    doci_dense.solve(numopt::eigenproblem::SolverType::DENSE);
    double doci_energy_dense = doci_dense.get_eigenvalue() + internuclear_repulsion_energy;


    // Calculate the DOCI energy using the Davidson algorithm
    ci::DOCI doci_davidson (so_basis, h2);
    doci_davidson.solve(numopt::eigenproblem::SolverType::DAVIDSON);
    double doci_energy_davidson = doci_davidson.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(doci_energy_dense - doci_energy_davidson) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( DOCI_beh_cation_klaas_Davidson ) {

    // Klaas' reference DOCI energy for BeH+ (obtained through Caitlin)
    double reference_doci_energy = -14.8782216937;


    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../tests/reference_data/beh_cation_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons
    doci.solve(numopt::eigenproblem::SolverType::DAVIDSON);


    // Calculate the total energy
    double internuclear_repulsion_energy = 1.5900757460937498e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_lih_klaas_Davidson ) {


    // Klaas' reference DOCI energy for LiH (obtained through Caitlin)
    double reference_doci_energy = -8.0029560313;


    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../tests/reference_data/lih_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons
    doci.solve(numopt::eigenproblem::SolverType::DAVIDSON);


    // Calculate the total energy
    double internuclear_repulsion_energy = 9.6074293445896852e-01;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}
