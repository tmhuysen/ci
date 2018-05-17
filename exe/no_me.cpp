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
    std::ofstream outfile ("no_fci_pb_test2.data");
    outfile << "bohr          " << " \t FCI energy   " << " \t electronic FCI energy   "<< " \t NO population   " << std::endl;
    outfile << std::setprecision(12);
    for(size_t i = 0;i<10;i++){
        std::vector<size_t> ao_set = {i};
        std::string filename = "../../exe/PB/no_2.3_PB";
        libwint::Molecule no (filename+".xyz",1);
        libwint::AOBasis ao_basis (no, "STO-3G");
        ao_basis.calculateIntegrals();
        hf::rhf::RHF rhf (no, ao_basis, 1.0e-06);
        rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);
        libwint::SOMullikenBasis so_basis (ao_basis, rhf.get_C_canonical());
        ci::FCI fci(so_basis,7,7);
        double energy = fci.solveConstrained(numopt::eigenproblem::SolverType::DAVIDSONLEMMENS,ao_set,0);
        double repulsion = no.calculateInternuclearRepulsionEnergy();
        double total_energy = energy+repulsion;
        double population = fci.get_population_set();
        outfile << std::to_string(i) << " \t " << total_energy << " \t "  << energy << " \t " << population<<std::endl;

    }
    outfile.close();
}
