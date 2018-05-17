/**
 *  Test constrained DOCI
 */


#include <libwint.hpp>
#include <numopt.hpp>
#include <hf.hpp>

#include <ci.hpp>
#include <iomanip>
#include <RHFC.hpp>


int main() {
    std::ofstream outfile ("co_base.cdoci");
    outfile << "doci multiplier" << " \t RHF energy   " << " \t C population RHF    "<< " \t C DOCI energy  " << " C population DOCI " <<std::endl;
    outfile << std::setprecision(12);
    libwint::Molecule CO ("../../exe/Cdoci/CO.xyz");
    libwint::AOBasis ao_basis (CO, "STO-3G");
    ao_basis.calculateIntegrals();
    double repulsion = CO.calculateInternuclearRepulsionEnergy();
    std::vector<size_t> ao_set = {0,1,2,3,4};
    for(double i = -2.0;i<2.0;i += 0.1){
        hf::rhf::RHFC rhfc (CO, ao_basis, 1.0e-12, 10000);
        rhfc.solve( hf::rhf::solver::SCFSolverType::DIIS,ao_set,0);
        double total_energy_rhfc = rhfc.get_electronic_energy() + repulsion;
        double population_rhfc = rhfc.get_population_set();
        libwint::SOMullikenBasis so_basis(ao_basis,rhfc.get_C_canonical());
        ci::DOCI doci(so_basis,CO);
        double energy = doci.solveConstrained(numopt::eigenproblem::SolverType::ARMADENSE,ao_set,i);
        double total_energy_doci = energy+repulsion;
        double population_doci = doci.get_population_set();
        outfile <<std::endl<< std::to_string(i) <<"\t" << total_energy_rhfc <<"\t"  << population_rhfc <<"\t"  << total_energy_doci <<"\t"  << population_doci;

    }
    outfile <<std::endl<<" repulsion " << CO.calculateInternuclearRepulsionEnergy()<< " rhf multiplier : "<<0;
    outfile.close();
}
