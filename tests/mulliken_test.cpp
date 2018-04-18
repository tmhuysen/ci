#define BOOST_TEST_MODULE "mulliken_test"



#include "DOCI.hpp"
#include "FCI.hpp"
#include <RHF.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( mulliken_fci ) {

    // Do a RHF calculation
    const std::string xyzfilename = "../tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule helium (xyzfilename);
    libwint::AOBasis ao_basis (helium, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (helium, ao_basis, threshold,1000);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());
    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,5,5);
    // Do da doci
    ci::DOCI doci(so_basis,10);
    double con = 8;
    std::vector<size_t> AO_set{0,1,2,3,4};
    doci.solve(numopt::eigenproblem::SolverType::DENSE);
    doci.compute1RDM();
    fci.solve(numopt::eigenproblem::SolverType::DENSE);
    fci.compute1RDM();
    so_basis.calculateMullikenMatrix(AO_set);
    Eigen::MatrixXd aa = fci.get_one_rdm_aa();
    Eigen::MatrixXd bb = fci.get_one_rdm_bb();
    std::cout<<std::endl;
    std::cout<<aa;
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<bb;
    std::cout<<std::endl;
    double pop = so_basis.mullikenPopulationCI(aa,bb);
    std::cout<<pop;


    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE ( test_exe ) {

    const std::string xyzfilename = "../exe/NO_2.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule helium (xyzfilename,+1);
    libwint::AOBasis ao_basis (helium, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (helium, ao_basis, 1e-3,1000000);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());
    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,7,7);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);
    BOOST_CHECK(true);

}
