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
    std::cout<<std::setprecision(16);
    std::vector<size_t> AO_set{0,1,2,3,4};
    std::string xyzfilename = st + std::to_string(1) + no;
    double threshold = 1.0e-012;
    std::string basis_name = "STO-3G";
    libwint::Molecule NO (xyzfilename,+1);
    libwint::AOBasis ao_basis (NO, basis_name);
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (NO, ao_basis, threshold,500000);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);
    libwint::SOMullikenBasis so_basis (ao_basis,rhf.get_C_canonical());
    int j = 0;
    ci::FCI fci(so_basis,7,7);
    fci.solveConstrained(numopt::eigenproblem::SolverType::DAVIDSON,AO_set,0);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(fci.get_one_rdm_aa());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solverS(ao_basis.get_S());
    Eigen::MatrixXd C = solverS.operatorInverseSqrt()*solver.eigenvectors();

    for(int i = 2;i<24;i++){
        std::string xyzfilename = st + std::to_string(i) + no;
        std::string basis_name = "STO-3G";
        libwint::Molecule NO (xyzfilename,+1);
        libwint::AOBasis ao_basis (NO, basis_name);
        ao_basis.calculateIntegrals();
        libwint::SOMullikenBasis so_basis (ao_basis,C);
        int j = 0;
        ci::FCI fci(so_basis,7,7);
        fci.solveConstrained(numopt::eigenproblem::SolverType::DAVIDSON,AO_set,0);
        double total_energy = fci.get_eigenvalue() + NO.calculateInternuclearRepulsionEnergy();
        double pop = fci.get_population();
        std::cout << " file "<<i<<" " << j << " & " << -pop + 7 << " & " << pop  << " & " << total_energy << " \\\\" << std::endl;
        std::cout << " file "<<i<<" " << j << " & " << -14+pop + 8 << " & " << 14-pop  << " & " << total_energy << " \\\\" << std::endl;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(fci.get_one_rdm_aa());
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solverS(ao_basis.get_S());
        C = solverS.operatorInverseSqrt()*solver.eigenvectors();


    }




}