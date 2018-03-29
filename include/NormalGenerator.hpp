#ifndef CI_NORMALGENERATOR_HPP
#define CI_NORMALGENERATOR_HPP

#include <random>


/**
 * Random number generator with normal distribution and set bounds.
 */
class NormalGenerator {
private:
    // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device random_device;

    // Mersenne twister prime number generator, initialized with seed from previous random device instance
    std::mt19937 generator;

    // a normal distribution
    std::normal_distribution<double> distribution;

public:
    void setMin_bound(double min_bound);

    void setMax_bound(double max_bound);

    void setMean(double mean);

    void setStd_dev(double std_dev);

private:
    // maximum and minimum values that may be generated
    double min_bound;
    double max_bound;

    //distribution parameters
    double mean;
    double std_dev;

public:
    /**
     * constructor for the generator
     */
    NormalGenerator(double mean, double std_dev, double min_bound,double max_bound);

    /**
     * generates a random number within the bounds with the normal distribution
     */
    double generate();

    /**
     * multiplies the bounds and deviation with @magnitude. To pinch or stretch your distribution.
     */
    void transform(double magnitude);


};


#endif // CI_NORMALGENERATOR_HPP
