#ifndef CI_SOPTIMIZER_HPP
#define CI_SOPTIMIZER_HPP

#include "ci.hpp"

namespace ci {


class SOptimizer {
private:
    libwint::SOMullikenBasis so_basis;
    libwint::SOMullikenBasis so_basis_new;
    numopt::eigenproblem::SolverType solver_type;
    size_t N;
    size_t K;
    double eigenvalue = 0;
    double temp_mod = 2;
    Eigen::VectorXd eigenvector;

    double population = 0;
    bool population_counts = false;
    double population_interval;
    double population_target;
public:

    double l_m = 0;
    std::vector<size_t> AO_set = {};

    /**
     * Constructors
     */

    SOptimizer(libwint::SOMullikenBasis so_basis, size_t N, numopt::eigenproblem::SolverType solver_type) : so_basis(
            so_basis), so_basis_new(so_basis), N(N), solver_type(solver_type), K(so_basis.get_K()) {

    };

    /**
     * Optimize based on parameters
     */

    void optimize(size_t max_iterations = 20000, size_t max_fails = 1000,
                  double theta_interval = 80.214, double narrowing_factor_temperature = 0.99,
                  double narrowing_factor_intervals = 0.9999, double std_dev = 25);

    void setNewToCurrent(){
        this->so_basis.copy(so_basis_new);
        this->temp_mod = 200;

    }
    void setPopulationTarget(double interval, double target){
        this->population_counts = true;
        this->population_interval = interval;
        this->population_target = target;

    }

    // GETTERS
    const libwint::SOMullikenBasis &get_so_basis() const { return so_basis; }
    const libwint::SOMullikenBasis &get_so_basis_new() const { return so_basis_new; }

    numopt::eigenproblem::SolverType get_solver_type() const { return solver_type; }

    double get_eigenvalue() const { return eigenvalue; }
    double get_population() const { return population; }

    const Eigen::VectorXd &get_eigenvector() const { return eigenvector; }


};

}   // ci
#endif //CI_SOPTIMIZER_HPP
