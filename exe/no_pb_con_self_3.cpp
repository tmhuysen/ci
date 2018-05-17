/**
 *  Calculate constrained NO+ from PB's basis
 */


#include <libwint.hpp>
#include <numopt.hpp>
#include <hf.hpp>

#include <ci.hpp>
#include <iomanip>


int main() {
    std::ofstream outfile ("NO_FCI_Constrained_self4.data");

    outfile << std::setprecision(12);
    std::vector<size_t> ao_set = {0,1,2,3,4};
    std::vector<double> files = {2.3,2.5,3.0,4.0,5.0,7.0,10.0};
    double d = 1;
    for(double i : files){
        outfile << "BOHR :"<<std::to_string(i)<< std::endl;
        outfile << "multiplier" << " \t FCI energy" << " \t electronic FCI energy"<< " \t NO population" << std::endl;
        try{
            std::stringstream stream;
            stream << std::fixed << std::setprecision(1) << i;
            std::string s = stream.str();
            std::string filename = "../../exe/PB/no_" + s +"_PB";
            libwint::SOMullikenBasis so_basis_temp(filename,10);
            so_basis_temp.calculateMullikenMatrix({0});
            so_basis_temp.set_lagrange_multiplier(0);
            libwint::Molecule no (filename+".xyz",1);
            libwint::AOBasis ao_basis(no,"STO-3G");
            ao_basis.calculateIntegrals();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigs(so_basis_temp.get_S());
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigs2(ao_basis.get_S());
            Eigen::MatrixXd c_corrected2 = so_basis_temp.get_C()*eigs.operatorSqrt();
            Eigen::MatrixXd dd = so_basis_temp.get_C().transpose();
            std::cout<<std::endl<<ao_basis.get_T()+ao_basis.get_V()<<std::endl;
            std::cout<<std::endl<<libwint::transformations::transform_SO_to_AO(so_basis_temp.get_h_SO(),dd)<<std::endl;
            std::cout<<"libint :"<<std::endl<<c_corrected2*c_corrected2.transpose()<<std::endl;
            Eigen::MatrixXd c_corrected = c_corrected2*eigs2.operatorInverseSqrt();
            Eigen::MatrixXd cano = c_corrected.transpose();
            libwint::SOMullikenBasis so_basis(ao_basis,cano);
            std::cout<<"libint :"<<std::endl<<so_basis.get_h_SO()<<std::endl;
            std::cout<<"pb :" <<  std::endl<<so_basis_temp.get_h_SO()<<std::endl;

            so_basis.set_P(false);
            /*
            ci::FCI fci(so_basis,7,7);
            double energy = fci.solveConstrained(numopt::eigenproblem::SolverType::DAVIDSONLEMMENS,ao_set,0);
            double repulsion = no.calculateInternuclearRepulsionEnergy();
            double total_energy = energy+repulsion;
            double population = fci.get_population_set();
            outfile << std::to_string(0) << " \t " << total_energy << " \t "  << energy << " \t " << population<<std::endl;
            */


        }catch (const std::exception& e) {
            outfile << std::to_string(0) <<e.what();

        }
        d+=0.33;

    }

    outfile.close();
}
