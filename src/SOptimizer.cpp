#include "SOptimizer.hpp"
namespace ci {

void SOptimizer::optimize(size_t max_iterations, size_t max_fails, double theta_interval,
                          double narrowing_factor_temperature, double narrowing_factor_intervals, double std_dev) {

    if (max_iterations < max_fails) {
        std::__throw_logic_error(
                "maximum amount of allowed consecutive failure should be lower than the maximum total iterations");
    }
    ci::DOCI doci(this->so_basis,this->N);
    double current_energy = doci.solveConstrained(this->solver_type,this->AO_set,this->l_m);
    // We need random generation to pick which two spatial orbitals we are gonna rotate.
    // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device random_device;
    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 random_generator(random_device());
    // Define range for SO's random selection.
    std::uniform_int_distribution<> distribution(0, this->K - 1);
    // Define range higher energy leniency
    std::uniform_real_distribution<> double_distribution(0.0, 1);

    size_t failures = 0; //counts the amount of consecutive non accepted rotations
    size_t iterations = 0; //counts the amount of iterations
    NormalGenerator normal_generator(0, std_dev, -theta_interval,
                                     theta_interval); //random generator around 0 with standard deviation of 25 and max and min values of 80.214.

    double temperature = std::abs(current_energy / this->temp_mod); //indicator of how lenient we accept a new energy that is to high.

    libwint::SOMullikenBasis basis_best = this->so_basis; //copy the old basis
    libwint::SOMullikenBasis basis_prev = this->so_basis; //copy the old basis

    while (iterations < max_iterations and failures < max_fails) {

        size_t so_index1 = 0; //spatial orbital index
        size_t so_index2 = 0; //spatial orbital index
        while (so_index1 == so_index2) { //if they are the same we need to re-pick
            so_index1 = distribution(random_generator);
            so_index2 = distribution(random_generator);
            if (so_index1 > so_index2) {
                size_t temp = so_index1;
                so_index1 = so_index2;
                so_index2 = temp;
            }
        }

        double angle = normal_generator.generate(); //generate random angle normal distributed
        basis_best.rotateJacobi(so_index1, so_index2, angle);
        ci::DOCI doci(basis_best, this->N);
        double new_energy = doci.solveConstrained(this->solver_type,this->AO_set,this->l_m);
        double energy_diff = current_energy - new_energy;

        if (energy_diff < 0) {
            //Algorithm in thesis of mario : annealing (pg 54)
            double random_number = double_distribution(random_generator);
            double exp_value = exp(energy_diff / temperature);
            double compare_value = exp_value / (exp_value + 1);


            if (random_number > compare_value) {
                basis_best.copy(basis_prev); // return to the old state.
                failures++;
            } else {
                current_energy = new_energy;
                failures = 0;
            }
        } else {
            current_energy = new_energy;
        }

        temperature *= narrowing_factor_temperature;
        normal_generator.transform(narrowing_factor_intervals);
        basis_prev.copy(basis_best);  //copy the old basis
        iterations++;
    }

    if(this->eigenvalue > current_energy){
        this->so_basis_new.copy(basis_best);
        this->eigenvalue = current_energy;
        this->eigenvector = doci.get_eigenvector();
        this->population = doci.get_population_set();
    }



}

}  // ci