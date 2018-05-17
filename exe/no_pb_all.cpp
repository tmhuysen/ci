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
    std::ofstream outfile ("no_fci_pb_test6.data");
    outfile << "bohr          " << " \t FCI energy   " << " \t electronic FCI energy   "<< " \t NO population   " << std::endl;
    outfile << std::setprecision(12);
    for(size_t i = 0;i<10;i++){
        std::vector<size_t> ao_set = {i};
        std::string filename = "../../exe/PB/no_10.0_PB";
        libwint::SOMullikenBasis so_basis(filename,10);
        libwint::Molecule no (filename+".xyz",1);
        ci::FCI fci(so_basis,7,7);
        double energy = fci.solveConstrained(numopt::eigenproblem::SolverType::DAVIDSONLEMMENS,ao_set,0);
        double repulsion = no.calculateInternuclearRepulsionEnergy();
        double total_energy = energy+repulsion;
        double population = fci.get_population_set();
        outfile <<population<<", ";

    }
    outfile.close();
}
