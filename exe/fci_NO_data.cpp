#include <chrono>
#include <iostream>

#include <Eigen/Dense>
#include <cpputil.hpp>

#include <numopt.hpp>
#include <libwint.hpp>
#include <ci.hpp>
#include <RHF.hpp>

int main () {
    // Do a RHF calculation
    const std::string xyzfilename = "../exe/NO_2.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule helium (xyzfilename,+1);
    libwint::AOBasis ao_basis (helium, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (helium, ao_basis, threshold,1000);
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