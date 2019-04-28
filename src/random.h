//
// Created by xgallom on 4/16/19.
//

#ifndef QM_MONTE_CARLO_RANDOM_H
#define QM_MONTE_CARLO_RANDOM_H

#include <random>

class Random {
public:
	Random();

	float get();
	float norm();
	float dist();

private:
	std::random_device m_device;
	std::mt19937 m_engine;
	std::uniform_real_distribution<float> m_dist = std::uniform_real_distribution<float>(-1., 1.);
	std::uniform_real_distribution<float> m_norm = std::uniform_real_distribution<float>();
};

#endif //QM_MONTE_CARLO_RANDOM_H
