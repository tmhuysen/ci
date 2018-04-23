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
    Eigen::Matrix<double,10,10> C;
    C <<  0.0 ,      0.9933096 ,     0.0 ,       0.0 ,   0.0 , 0.0 ,         0.2677014,    0.0120400 ,   0.0,        0.0,
         -0.0 ,      0.0269998 ,    -0.0 ,       0.0 ,   0.0 , -0.0000001,   -1.0273904,  -0.0470429 ,   0.0,        0.0,
          0.0 ,      -0.0002174 ,    0.0 ,       0.0 ,   0.0 , 0.0000004  ,  0.0457351,   -0.9989536  ,  0.0,        0.0,
          0.0 ,      0.0 ,           0.0 ,       0.0 ,   0.0 , 0.0  ,        0.0 ,         0.0   ,       0.0,        1.0,
          0.0 ,      0.0 ,           0.0 ,       0.0 ,   0.0 , 0.0,          0.0,          0.0,          1.0,        0.0,
          0.9966037, 0.0  ,          0.2571079 , 0.0 ,   0.0 , 0.0055501,    -0.0,         0.0,          0.0,        0.0,
          0.0139596, -0.0 ,          -1.0289185, 0.0 ,   0.0 , -0.0220485,   0.0000001,    -0.0000002,   0.0,        0.0,
          0.0000406, -0.0 ,          0.0214244,  0.0 ,   0.0 , -0.9997705,   0.0000003,    -0.0000007,   0.0,        0.0,
          0.0   ,    0.0 ,           0.0 ,       0.0 ,   1.0 , 0.0,          0.0,          0.0,          0.0,       -0.0,
          0.0   ,    0.0 ,           0.0 ,       1.0 ,   0.0 , 0.0,          0.0,          0.0,         -0.0,        0.0;
    Eigen::MatrixXd CC = C;

    std::string fileout3 = "../STO-3G_NO+_INTEGRALS.out";
    std::ofstream outfile3 (fileout3);
    std::string st = "../exe/";
    std::string no = "NO.xyz";
    std ::cout<<std::setprecision(16);
    std::vector<size_t> AO_set{0,1,2,3,4};
    std::string xyzfilename = "2NO.xyz";
    double threshold = 1.0e-012;
    std::string basis_name = "STO-3G";
    libwint::Molecule NO (xyzfilename,+1);
    libwint::AOBasis ao_basis (NO, basis_name);
    ao_basis.calculateIntegrals();
    outfile3<<std::setprecision(15);
    outfile3<<std::endl<<" S MATRIX "<<std::endl<<ao_basis.get_S()<<std::endl;
    outfile3<<std::endl<<" T MATRIX "<<std::endl<<ao_basis.get_T()<<std::endl;
    outfile3<<std::endl<<" G MATRIX "<<std::endl<<ao_basis.get_g()<<std::endl;


    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(ao_basis.get_S());
    Eigen::MatrixXd S_noise = solver.operatorInverseSqrt();
    int sign = 1;
    for(int i = 0; i<CC.cols();i++){
        for(int j = 0; j<CC.cols();j++) {

            S_noise(i,j) += sign*0.00000005;
            sign *= -1;
        }

    }
    libwint::SOMullikenBasis so_basis (ao_basis,(S_noise*CC));
    int j = 0;
    /*
    ci::FCI fci(so_basis,7,7);
    fci.solveConstrained(numopt::eigenproblem::SolverType::ARMADENSE,AO_set,0);
    double total_energy = fci.get_eigenvalue() + NO.calculateInternuclearRepulsionEnergy();
    double pop = fci.get_population();
    outfile3 << "file 8" << std::endl;
    outfile3 << CC << std::endl;
    outfile3 <<" " << j << " & " << -pop + 7 << " & " << pop  << " & " << total_energy << " \\\\" << std::endl;
    outfile3 <<" " << j << " & " << -14+pop + 8 << " & " << 14-pop  << " & " << total_energy << " \\\\" << std::endl;
    outfile3 <<std::endl;
    */
    outfile3.close();

}