/**
 *  Create FCI reference data for a given h2 geometry @STO-3G.
 *  For a 2-electron system DOCI is exact, so we can use the current DOCI implementation instead of a FCI one.
 */


#include <libwint.hpp>
#include <numopt.hpp>
#include <hf.hpp>

#include <ci.hpp>
#include <iomanip>


int main() {
    std::ofstream outfile ("NO_FCI_02569_5.data");
    outfile << "bohr          " << " \t FCI energy   " << " \t electronic FCI energy   "<< " \t NO population   " << std::endl;
    outfile << std::setprecision(12);
    std::vector<size_t> ao_set = {0,1,2,3,4};
    for(double i = 0.9;i<10.1;i += 0.1){
        std::stringstream stream;
        stream << std::fixed << std::setprecision(1) << i;
        std::string s = stream.str();
        std::string filename = "../../exe/PB/no_" + s +"_PB";
        libwint::SOMullikenBasis so_basis(filename,10);
        so_basis.set_C(so_basis.get_C().transpose());
        so_basis.set_P(false);
        libwint::Molecule no (filename+".xyz",1);
        ci::FCI fci(so_basis,7,7);
        double energy = fci.solveConstrained(numopt::eigenproblem::SolverType::DAVIDSONLEMMENS,ao_set,0);
        double repulsion = no.calculateInternuclearRepulsionEnergy();
        double total_energy = energy+repulsion;
        double population = fci.get_population_set();
        outfile << std::to_string(i) << " \t " << total_energy << " \t "  << energy << " \t " << population<<std::endl;

    }
    outfile.close();
}
