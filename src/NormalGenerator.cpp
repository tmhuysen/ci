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
	std::__throw_overflow_error("could not generate a random number for given bounds within a reasonable amount of attempts");
}

void NormalGenerator::transform(double magnitude) {
	this->min_bound *= magnitude;
	this->max_bound *= magnitude;
	this->std_dev *= magnitude;
	this->distribution = std::normal_distribution<double>(this->mean,this->std_dev);
}
