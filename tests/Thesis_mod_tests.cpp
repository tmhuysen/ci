#define BOOST_TEST_MODULE "Thesis_mod_tests"



#include "ci.hpp"
#include "RHF.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>
#include <SOptimizer.hpp>


BOOST_AUTO_TEST_CASE ( Constructor_mulliken_type_basis_test ) {

    // Prepare an SO basis to test the constructor
    libwint::Molecule water ("../tests/reference_data/h2o_Psi4_GAMESS.xyz");
    libwint::AOBasis ao_basis (water, "STO-3G");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (water, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOMullikenBasis so_basis (ao_basis, rhf.get_C_canonical());
    libwint::SOMullikenBasis so_basis2 (ao_basis, rhf.get_C_canonical());
    ci::FCI fci (so_basis, 5, 5);
    ci::DOCI doci (so_basis, water);

    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE ( optimize_doci ) {

    // Prepare an SO basis to test the constructor
    libwint::Molecule water ("../tests/reference_data/h2o_Psi4_GAMESS.xyz");
    libwint::AOBasis ao_basis (water, "STO-3G");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (water, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOMullikenBasis so_basis (ao_basis, rhf.get_C_canonical());
    ci::DOCI doci (so_basis, water);
    ci::FCI fci (so_basis, 5,5);
    doci.solve(numopt::eigenproblem::SolverType::DENSE);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);
    ci::SOptimizer sop (so_basis, 10, numopt::eigenproblem::SolverType::DENSE);
    sop.optimize(20000, 1000);
    std::cout<<sop.get_eigenvalue()+water.calculateInternuclearRepulsionEnergy()<<std::endl;
    std::cout<<doci.get_eigenvalue()+water.calculateInternuclearRepulsionEnergy()<<std::endl;
    std::cout<<fci.get_eigenvalue()+water.calculateInternuclearRepulsionEnergy()<<std::endl;
    BOOST_CHECK(true);
}