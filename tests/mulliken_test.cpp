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
    rhf.solve();
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());
    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,5,5);
    // Do da doci
    ci::DOCI doci(so_basis,10);
    double con = 10;
    std::vector<size_t> AO_set{0,1,2,3,4};
    doci.solve(numopt::eigenproblem::SolverType::DENSE);
    //std::cout<<std::endl<<" DOCI : "<<doci.solveConstrained(numopt::eigenproblem::SolverType::DENSE,AO_set,con);
    std::cout<<std::endl<<" FCI : "<<fci.solveConstrained(numopt::eigenproblem::SolverType::DENSE,AO_set,con);


    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE ( test_exe ) {

    const std::string xyzfilename = "../exe/NO.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule helium (xyzfilename,+1);
    libwint::AOBasis ao_basis (helium, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (helium, ao_basis, 1e-3,1000000);
    rhf.solve();
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());
    // Do a FCI calculation based on a given SObasis
    ci::FCI fci(so_basis,7,7);
    ci::DOCI doci(so_basis,14);
    std::vector<size_t> AO_set{0,1,2,3,4};
    for(double con = 5; con<8.1; con += 0.1){
        std::cout<<std::endl<<" Con : "<<con;
        std::cout<<std::endl<<" FCI : "<<fci.solveConstrained(numopt::eigenproblem::SolverType::DENSE,AO_set,con);
        std::cout<<std::endl<<" DOCI : "<<doci.solveConstrained(numopt::eigenproblem::SolverType::DENSE,AO_set,con);
    }
}
