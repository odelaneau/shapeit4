#ifndef _RANDOM_NUMBER_H
#define _RANDOM_NUMBER_H

#include <cfloat>
#include <cstdint>
#include <random>
#include <vector>


class random_number_generator {
protected:
	unsigned int seed;
	std::mt19937 randomEngine;
	std::uniform_int_distribution < unsigned int > uniformDistributionInt;
	std::uniform_real_distribution < double > uniformDistributionDouble;

public:

	random_number_generator(unsigned int seed = 15052011) : randomEngine(seed), uniformDistributionInt(0, 32768), uniformDistributionDouble(0, 1.0) {
	}

	~random_number_generator(){
	}

	void setSeed(unsigned int _seed) {
		seed = _seed;
		randomEngine.seed(seed);
	}

	unsigned int getSeed() {
		return seed;
	}

	std::mt19937 & getEngine() {
		return randomEngine;
	}

	unsigned int getInt(unsigned int imin, unsigned int imax) {
		return uniformDistributionInt(randomEngine, std::uniform_int_distribution < unsigned int > {imin, imax}.param());
	}

	unsigned int getInt(unsigned int isize) {
		return getInt(0, isize - 1);
	}

	double getDouble(double fmin, double fmax) {
		return uniformDistributionDouble(randomEngine, std::uniform_real_distribution < double > {fmin, fmax}.param());
	}

	double getDouble() {
		return getDouble(0.0, 1.0);
	}

	bool flipCoin() {
		return (getDouble() < 0.5);
	}

	int sample(std::vector < double > & vec, double sum) {
		double csum = vec[0];
		double u = getDouble() * sum;
		for (int i = 0; i < vec.size() - 1; ++i) {
			if ( u < csum ) return i;
			csum += vec[i+1];
		}
		return vec.size() - 1;
	}

	int sample4(const double * vec, double sum) {
		double csum = vec[0];
		double u = getDouble() * sum;
		for (int i = 0; i < 3; ++i) {
			if ( u < csum ) return i;
			csum += vec[i+1];
		}
		return 3;
	}
};

#endif
