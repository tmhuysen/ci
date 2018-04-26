#define BOOST_TEST_MODULE "FCI_dense_test"

#include "FCI.hpp"
#include <RHF.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( FCI_He_Cristina ) {
    // Cristina's He FCI energy
    double reference_fci_energy = -2.902533599;

    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/he_cristina.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "aug-cc-pVQZ";
    libwint::Molecule helium (xyzfilename);
    libwint::AOBasis ao_basis (helium, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (helium, ao_basis, threshold);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis,rhf.get_C_canonical());

    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,1,1);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);

    // Calculate the total energy
    double internuclear_repulsion_energy = helium.calculateInternuclearRepulsionEnergy();  // this comes straight out of the FCIDUMP file
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}

BOOST_AUTO_TEST_CASE ( FCI_H2_Cristina ) {
    // Cristina's H2 FCI energy/OO-DOCI energy
    double reference_fci_energy = -1.1651486697;

    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/h2_cristina.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "6-31g**";
    libwint::Molecule hydrogen_gas (xyzfilename);
    libwint::AOBasis ao_basis (hydrogen_gas, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (hydrogen_gas, ao_basis, threshold);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis,rhf.get_C_canonical());

    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,1,1);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);

    // Calculate the total energy
    double internuclear_repulsion_energy = hydrogen_gas.calculateInternuclearRepulsionEnergy();  // this comes straight out of the FCIDUMP file
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;
    std::cout<<std::endl<<test_fci_energy<<std::endl;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}

BOOST_AUTO_TEST_CASE ( FCI_H2O_Psi4_Games ) {
    // Psi4 and Games's FCI energy
    double reference_fci_energy = -75.0129803939602;

    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/h2o_Psi4_Games.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule water (xyzfilename);
    libwint::AOBasis ao_basis (water, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (water, ao_basis, threshold);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis,rhf.get_C_canonical());

    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,5,5);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);

    // Calculate the total energy
    double internuclear_repulsion_energy = water.calculateInternuclearRepulsionEnergy();  // this comes straight out of the FCIDUMP file
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;
    std::cout<<std::endl<<test_fci_energy<<std::endl;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}

BOOST_AUTO_TEST_CASE ( FCI_H2O_Psi4_Games ) {
    // Cristina's H2 FCI energy/OO-DOCI energy
    double reference_fci_energy = -75.0129803939602;

    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/h2o_Psi4_Games.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule water (xyzfilename);
    libwint::AOBasis ao_basis (water, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (water, ao_basis, threshold);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis,rhf.get_C_canonical());

    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,1,1);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);

    // Calculate the total energy
    double internuclear_repulsion_energy = water.calculateInternuclearRepulsionEnergy();  // this comes straight out of the FCIDUMP file
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}