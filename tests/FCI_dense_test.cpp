#define BOOST_TEST_MODULE "FCI_dense_test"

#include "FCI.hpp"
#include <RHF.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

/*
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
    rhf.solve(hf::rhf::solver::SCFSolverType::PLAIN);
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());

    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,1,1);
    fci.solve(numopt::eigenproblem::SolverType::DAVIDSON);

    // Calculate the total energy
    double internuclear_repulsion_energy = helium.calculateInternuclearRepulsionEnergy();  // this comes straight out of the FCIDUMP file
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}
*/
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
    rhf.solve(hf::rhf::solver::SCFSolverType::PLAIN);
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());


    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,1,1);
    fci.solve(numopt::eigenproblem::SolverType::DAVIDSON);

    // Calculate the total energy
    double internuclear_repulsion_energy = hydrogen_gas.calculateInternuclearRepulsionEnergy();  // this comes straight out of the FCIDUMP file
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;



    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}

BOOST_AUTO_TEST_CASE ( FCI_NO) {

    std::cout<<std::setprecision(16);
    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/1NO.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule hydrogen_gas (xyzfilename,1);
    libwint::AOBasis ao_basis (hydrogen_gas, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (hydrogen_gas, ao_basis, threshold);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());

    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,7,7);

    std::cout<<" DO I GET HERE??";

    fci.solve(numopt::eigenproblem::SolverType::DAVIDSON);

    //Calculate the total energy
    double internuclear_repulsion_energy = hydrogen_gas.calculateInternuclearRepulsionEnergy();  // this comes straight out of the FCIDUMP file
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;





    std::cout<<test_fci_energy<< " HERE ";

}

BOOST_AUTO_TEST_CASE ( FCI_H2O) {


    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule hydrogen_gas (xyzfilename);
    libwint::AOBasis ao_basis (hydrogen_gas, basis_name);
    ao_basis.calculateIntegrals();
    //hf::rhf::RHF rhf (hydrogen_gas, ao_basis, threshold);
    //rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);
    Eigen::MatrixXd F = ao_basis.get_V()+ao_basis.get_T();
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gsaes(F,ao_basis.get_S());
    Eigen::MatrixXd C = gsaes.eigenvectors();
    libwint::SOMullikenBasis so_basis (ao_basis,C);
    // Do a FCI calculation based on a given SObasis
    so_basis.rotateJacobi(3,6,14);
    so_basis.rotateJacobi(3,6,14);
    ci::FCI fci(so_basis,5,5);

    fci.solveConstrained(numopt::eigenproblem::SolverType::DENSE,{1},0);
    std::cout<<std::setprecision(16);
    double coeff = 0;
    size_t aa = 0;
    for (size_t a= 0;a<441;a++){
        if(coeff<fci.get_eigenvector()(a)){
            coeff = fci.get_eigenvector()(a);
            aa = a;

        }
        a++;
    }
    std::cout<<"add and coef "<<aa<<" "<<coeff<<std::endl;
    // Calculate the total energy
    double internuclear_repulsion_energy = hydrogen_gas.calculateInternuclearRepulsionEnergy();  // this comes straight out of the FCIDUMP file
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;



    std::cout<<test_fci_energy<< " HERE ";
}

BOOST_AUTO_TEST_CASE ( FCI_H2O_2) {


    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule hydrogen_gas (xyzfilename);
    libwint::AOBasis ao_basis (hydrogen_gas, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (hydrogen_gas, ao_basis, threshold);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);

    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());
    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,5,5);
    fci.solveConstrained(numopt::eigenproblem::SolverType::DAVIDSON,{1},0);
    std::cout<<std::setprecision(16);
    // Calculate the total energy
    double internuclear_repulsion_energy = hydrogen_gas.calculateInternuclearRepulsionEnergy();  // this comes straight out of the FCIDUMP file
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;
    std::cout<<fci.get_eigenvector();


    std::cout<<test_fci_energy<< " HERE ";
}