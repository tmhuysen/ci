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
    for(int i = 4;i<7;i++){
        std::string xyzfilename = st + std::to_string(i) + no;
        double threshold = 1.0e-02;
        std::string basis_name = "STO-3G";
        libwint::Molecule NO (xyzfilename,+1);
        libwint::AOBasis ao_basis (NO, basis_name);
        ao_basis.calculateIntegrals();
        Eigen::MatrixXd F = ao_basis.get_V()+ao_basis.get_T();
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gsaes(F,ao_basis.get_S());
        Eigen::MatrixXd C = gsaes.eigenvectors();
        libwint::SOMullikenBasis so_basis (ao_basis,C);
        int j = 0;
        ci::FCI fci(so_basis,7,7);
        fci.solve(numopt::eigenproblem::SolverType::DENSE);
        double total_energy = fci.get_eigenvalue() + NO.calculateInternuclearRepulsionEnergy();
        double pop = fci.get_population();
        std::cout <<  j << " & " << -pop + 7 << " & " << pop  << " & " << total_energy << " \\\\" << std::endl;




    }




}