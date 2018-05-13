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
    Eigen::VectorXd eigenvector;

    double l_m = 0;
    std::vector<size_t> AO_set = {};
    double population = 0;

public:

    /**
     * Constructors
     */

    SOptimizer(libwint::SOMullikenBasis so_basis, size_t N, numopt::eigenproblem::SolverType solver_type) : so_basis(
            so_basis), so_basis_new(libwint::SOMullikenBasis(so_basis.get_K())), N(N), solver_type(solver_type), K(so_basis.get_K()) {
        so_basis_new.copy(so_basis);
    };

    /**
     * Optimize based on parameters
     */

    void optimize(size_t max_iterations = 20000, size_t max_fails = 1000,
                  double theta_interval = 80.214, double narrowing_factor_temperature = 0.99,
                  double narrowing_factor_intervals = 0.9999, double std_dev = 25);

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
