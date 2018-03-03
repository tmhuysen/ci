#define BOOST_TEST_MODULE "FCI_test"

#include "FCI.hpp"

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
    libwint::Basis basis (helium, basis_name);
    hf::rhf::RHF rhf (basis, threshold);

    // Do a DOCI calculation based on the RHF calculation
    doci::CI_basis ci_basis (rhf);
    doci::FCI fci_test(&ci_basis);

    // Calculate the total energy as the sum of the lowest energy eigenstate + the internuclear repulsion
    const doci::State& ground = fci_test.getLowestEigenState();
    double test_fci_energy = ground.getEval() + ci_basis.getInternuclear_repulsion();

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
    libwint::Basis basis (hydrogen_gas, basis_name);
    hf::rhf::RHF rhf (basis, threshold);

    // Do a DOCI calculation based on the RHF calculation
    doci::CI_basis ci_basis (rhf);
    doci::FCI fci_test(&ci_basis);

    // Calculate the total energy as the sum of the lowest energy eigenstate + the internuclear repulsion
    const doci::State& ground = fci_test.getLowestEigenState();
    double test_fci_energy = ground.getEval() + ci_basis.getInternuclear_repulsion();

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);

}

