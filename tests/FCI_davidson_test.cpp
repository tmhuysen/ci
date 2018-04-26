#define BOOST_TEST_MODULE "FCI_dense_test"

#include "FCI.hpp"
#include <RHF.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( FCI_NO_Psi4_Games ) {
    // Psi4 and Games's NO+ FCI energy
    double reference_fci_energy = -127.4123904082415;

    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/NO+_Psi4_Games.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule NO (xyzfilename,1);
    libwint::AOBasis ao_basis (NO, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (NO, ao_basis, threshold);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis,rhf.get_C_canonical());

    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,7,7);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);

    // Calculate the total energy
    double internuclear_repulsion_energy = NO.calculateInternuclearRepulsionEnergy();  // this comes straight out of the FCIDUMP file
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}
