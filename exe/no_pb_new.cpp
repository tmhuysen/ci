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
    std::ofstream outfile ("no_fci_pb_test_tep2.data");
    outfile << "bohr          " << " \t FCI energy   " << " \t electronic FCI energy   "<< " \t NO population   " << std::endl;
    outfile << std::setprecision(12);
    std::vector<size_t> ao_set = {0};
    std::string filename = "../../exe/PB/no_10.0_PB";
    libwint::SOMullikenBasis so_basis(filename,10);

    libwint::Molecule no (filename+".xyz",1);

    so_basis.set_P(false);
    so_basis.calculateMullikenMatrix(ao_set);
    ci::FCI fci(so_basis,7,7);

    fci.solve(numopt::eigenproblem::SolverType::DAVIDSONLEMMENS);


    fci.calculate1RDMs();
    Eigen::MatrixXd aa = fci.get_one_rdm_aa();
    Eigen::MatrixXd bb = fci.get_one_rdm_bb();
    outfile << std::endl<< aa << std::endl;
    outfile << std::endl<< bb << std::endl;
    outfile << std::endl<< "0  population  " << so_basis.mullikenPopulationCI(aa,bb) << std::endl;
    double total = so_basis.mullikenPopulationCI(aa,bb);
    for(size_t i = 1;i<10;i++){
        std::vector<size_t> ao_set = {i};
        so_basis.calculateMullikenMatrix(ao_set);
        outfile <<", " <<so_basis.mullikenPopulationCI(aa,bb);
        total += so_basis.mullikenPopulationCI(aa,bb);








    }
    outfile << std::endl<< "population" << total << std::endl;
    outfile.close();
}
