#define BOOST_TEST_MODULE "RDM_test"


/*
#include "DOCI.hpp"
#include "FCI.hpp"
#include <RHF.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( compute_doci_test ) {

    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOMullikenBasis so_basis ("../tests/reference_data/beh_cation_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons
    doci.solve(numopt::eigenproblem::SolverType::DENSE);
    doci.compute1RDM();
    doci.compute2RDM();

    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE ( compute_fci_test ) {

    // Cristina's He FCI energy
    double reference_fci_energy = -2.902533599;

    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule helium (xyzfilename);
    libwint::AOBasis ao_basis (helium, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (helium, ao_basis, threshold);
    rhf.solve();
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());

    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,5,5);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);
    fci.compute1RDM();
    Eigen::MatrixXd rdm1 = fci.get_one_rdm_aa();
    Eigen::MatrixXd rdm2 = fci.get_one_rdm_bb();
    // Do da doci
    ci::DOCI doci(so_basis,10);
    doci.solve(numopt::eigenproblem::SolverType::DENSE);
    doci.compute1RDM();
    Eigen::MatrixXd rdm3 = doci.get_one_rdm_aa();
    Eigen::MatrixXd rdm4 = doci.get_one_rdm_bb();

    // fake RDM
    Eigen::MatrixXd rdm5 = Eigen::MatrixXd::Zero(7,7);
    rdm5.topLeftCorner(5,5) = Eigen::MatrixXd::Identity(5,5);

    // Mullikenbasis
    libwint::SOMullikenBasis soo(ao_basis,rhf.get_C_canonical());
    std::cout<<std::setprecision(12);
    double lol = 0;
    double lol2 = 0;
    for(size_t i = 0; i<7;i++){
        soo.calculateMullikenMatrix({i});
        lol += soo.mullikenPopulationFCI(rdm1, rdm2);
        lol2 += soo.mullikenPopulationCI(rdm1, rdm2);
        std::cout<<i<<" fci : "<<soo.mullikenPopulationCI(rdm1, rdm2)<<std::endl;
        std::cout<<i<<" Ffci : "<<soo.mullikenPopulationFCI(rdm1, rdm2)<<std::endl;
        std::cout<<i<<" doci : "<<soo.mullikenPopulationCI(rdm3, rdm4)<<std::endl;
        std::cout<<i<<" HF : "<<soo.mullikenPopulationCI(rdm5, rdm5)<<std::endl;
    }
    std::cout<<" LOL : "<<lol<<std::endl;
    std::cout<<" LOL2 : "<<lol2<<std::endl;



    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE ( reference_tests ) {

    // TODO: add reference
    BOOST_CHECK(false);
}

BOOST_AUTO_TEST_CASE ( test_rdm_trace ) {

    // Cristina's He FCI energy
    double reference_fci_energy = -2.902533599;

    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule helium (xyzfilename);
    libwint::AOBasis ao_basis (helium, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (helium, ao_basis, threshold);
    rhf.solve();
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());

    // Do da doci
    ci::DOCI doci(so_basis,10);
    doci.solve(numopt::eigenproblem::SolverType::DENSE);
    double test_doci_energy = doci.get_eigenvalue();
    std::cout<<std::setprecision(16);
    std::cout<<" trace : "<<doci.computeE();
    std::cout<<" safe : "<<test_doci_energy;




    BOOST_CHECK(true);
}*/