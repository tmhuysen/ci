/**
 *  Calculate constrained NO+ from PB's basis
 */


#include <libwint.hpp>
#include <numopt.hpp>
#include <hf.hpp>

#include <ci.hpp>
#include <iomanip>


int main() {
    std::ofstream outfile ("NO_FCI_Constrained3.data");

    outfile << std::setprecision(12);
    std::vector<size_t> ao_set = {0,1,2,3,4};
    std::vector<double> files = {2.3,2.5,3.0,4.0,5.0,7.0,10.0};
    double d = 1;
    for(double i : files){
        outfile << "BOHR :"<<std::to_string(i)<< std::endl;
        outfile << "multiplier" << " \t FCI energy" << " \t electronic FCI energy"<< " \t NO population" << std::endl;
        for(double x = -d;x<d+0.1;x += 0.1){
            try{
                std::stringstream stream;
                stream << std::fixed << std::setprecision(1) << i;
                std::string s = stream.str();
                std::string filename = "../../exe/PB/no_" + s +"_PB";
                libwint::SOMullikenBasis so_basis(filename,10);
                so_basis.set_C(so_basis.get_C().transpose());
                so_basis.set_P(false);
                libwint::Molecule no (filename+".xyz",1);
                ci::FCI fci(so_basis,7,7);
                double energy = fci.solveConstrained(numopt::eigenproblem::SolverType::DAVIDSONLEMMENS,ao_set,x);
                double repulsion = no.calculateInternuclearRepulsionEnergy();
                double total_energy = energy+repulsion;
                double population = fci.get_population_set();
                outfile << std::to_string(x) << " \t " << total_energy << " \t "  << energy << " \t " << population<<std::endl;


            }catch (const std::exception& e) {
                outfile<<e.what()<<std::to_string(x)<<std::endl;

            }


        }
        d+=0.33;

    }

    outfile.close();
}
