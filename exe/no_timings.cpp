#include <chrono>
#include <iostream>

#include <Eigen/Dense>
#include <cpputil.hpp>

#include <numopt.hpp>
#include <libwint.hpp>
#include <ci.hpp>
#include <RHF.hpp>
#include <iomanip>


int main () {
    std::string st = "../exe/";
    std::string no = "NO.xyz";
    std::vector<size_t> AO_set{0,1,2,3,4};
    std::cout<<std::setprecision(16);
    std::string xyzfilename = st + std::to_string(1) + no;
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule NO (xyzfilename,+1);
    libwint::AOBasis ao_basis (NO, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (NO, ao_basis, threshold,1000);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());
    ci::FCI fci(so_basis,7,7);
    double energy = fci.solveConstrained(numopt::eigenproblem::SolverType::DENSE,AO_set,0);
    double total_energy = energy + NO.calculateInternuclearRepulsionEnergy();
    double pop = fci.get_population();
    std::cout <<  0 << " & " << -pop + 7 << " & " << pop  << " & " << total_energy << " \\\\" << std::endl;





}