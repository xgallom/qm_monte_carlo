//
// Created by xgallom on 4/29/19.
//

#ifndef QM_MONTE_CARLO_MEAN_H
#define QM_MONTE_CARLO_MEAN_H

#include "context.h"
#include <iostream>

static float accumulate(const Mm *e, size_t n)
{
	Mm sum = _mm256_setzero_ps();

	while(n--)
		sum = _mm256_add_ps(sum, *e++);

	return sum[0] + sum[1] + sum[2] + sum[3] +
		   sum[4] + sum[5] + sum[6] + sum[7];
}

static float mean(const Mm *e, size_t n, size_t count)
{
	return accumulate(e, n) / count;
}

#endif //QM_MONTE_CARLO_MEAN_H
