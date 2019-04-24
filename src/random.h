//
// Created by xgallom on 4/16/19.
//

#ifndef QM_MONTE_CARLO_RANDOM_H
#define QM_MONTE_CARLO_RANDOM_H

#include "vector.h"

namespace Random
{
	void init();

	double get();
	double norm();
	double dist();
	Vector3D vector();
	Vector3D vectorOffset();
}

#endif //QM_MONTE_CARLO_RANDOM_H
