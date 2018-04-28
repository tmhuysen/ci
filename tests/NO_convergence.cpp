#define BOOST_TEST_MODULE "NO_convict"

#include "ci.hpp"
#include <RHF.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( NO_convergence_lemmens ) {

    std::cout<<std::setprecision(16);
    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/1NO.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule NO (xyzfilename,1);
    libwint::AOBasis ao_basis (NO, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (NO, ao_basis, threshold);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());

    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,7,7);
    fci.solve(numopt::eigenproblem::SolverType::DAVIDSONLEMMENS);

    //Calculate the total energy
    double internuclear_repulsion_energy = NO.calculateInternuclearRepulsionEnergy();  // this comes straight out of the FCIDUMP file
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;

    std::cout<<test_fci_energy;
}

BOOST_AUTO_TEST_CASE ( NO_convergence_lemmens_plain ) {

    std::cout<<std::setprecision(16);
    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/1NO.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule NO (xyzfilename,1);
    libwint::AOBasis ao_basis (NO, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (NO, ao_basis, threshold);
    rhf.solve(hf::rhf::solver::SCFSolverType::PLAIN);
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());

    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,7,7);
    fci.solve(numopt::eigenproblem::SolverType::DAVIDSONLEMMENS);

    //Calculate the total energy
    double internuclear_repulsion_energy = NO.calculateInternuclearRepulsionEnergy();  // this comes straight out of the FCIDUMP file
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;

    std::cout<<test_fci_energy;
}

BOOST_AUTO_TEST_CASE ( NO_convergence_plain ) {

    std::cout<<std::setprecision(16);
    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/1NO.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule NO (xyzfilename,1);
    libwint::AOBasis ao_basis (NO, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (NO, ao_basis, threshold);
    rhf.solve(hf::rhf::solver::SCFSolverType::PLAIN);
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());

    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,7,7);
    fci.solve(numopt::eigenproblem::SolverType::DAVIDSON);

    //Calculate the total energy
    double internuclear_repulsion_energy = NO.calculateInternuclearRepulsionEnergy();  // this comes straight out of the FCIDUMP file
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;

    std::cout<<test_fci_energy;
}

BOOST_AUTO_TEST_CASE ( NO_convergence_DIIS ) {

    std::cout<<std::setprecision(16);
    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/1NO.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule NO (xyzfilename,1);
    libwint::AOBasis ao_basis (NO, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (NO, ao_basis, threshold);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());

    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,7,7);
    fci.solve(numopt::eigenproblem::SolverType::DAVIDSON);

    //Calculate the total energy
    double internuclear_repulsion_energy = NO.calculateInternuclearRepulsionEnergy();  // this comes straight out of the FCIDUMP file
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;

    std::cout<<test_fci_energy;
}