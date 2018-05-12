#include <iostream>
#include "NormalGenerator.hpp"

NormalGenerator::NormalGenerator(double mean, double std_dev, double min_bound, double max_bound):
        distribution(mean,std_dev), mean(mean), std_dev(std_dev), min_bound(min_bound), max_bound(max_bound), generator(random_device())
{}

double NormalGenerator::generate() {
    size_t safety = 0; //prevent loop
    while (safety < 10000) {
        double number = this->distribution(generator);
        if (number >= this->min_bound && number <= this->max_bound)
            return number;
        safety++;
    }
    std::cout<<std::endl<< " max: "<<max_bound<<" min: "<<min_bound<<" mean : "<<mean<<" std : "<<std_dev<<std::endl;
    throw std::overflow_error("could not generate a random number for given bounds within a reasonable amount of attempts");
}

void NormalGenerator::transform(double magnitude) {
    double interval = this->mean - this->min_bound;
    interval *= magnitude;
    this->min_bound = mean - interval;

    double interval_max = this->max_bound - this->mean;
    interval_max *= magnitude;
    this->max_bound = mean + interval_max;

    this->std_dev *= magnitude;
    this->distribution = std::normal_distribution<double>(this->mean,this->std_dev);
}

void NormalGenerator::setMin_bound(double min_bound) {
    this->min_bound = min_bound;
}

void NormalGenerator::setMax_bound(double max_bound) {
    this->max_bound = max_bound;
}

void NormalGenerator::setMean(double mean) {
    this->mean = mean;
    this->distribution = std::normal_distribution<double>(this->mean,this->std_dev);
}

void NormalGenerator::setStd_dev(double std_dev) {
    this->std_dev = std_dev;
    this->distribution = std::normal_distribution<double>(this->mean,this->std_dev);
}
