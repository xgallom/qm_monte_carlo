//
// Created by xgallom on 4/29/19.
//

#ifndef QM_MONTE_CARLO_PRECOMPUTE_H
#define QM_MONTE_CARLO_PRECOMPUTE_H

#include "arithmetic.h"

static void computeWaveFunc(Mm *to, const Mm *r1, const Mm *r2, size_t n, const float a, const float b)
{
	const Mm
			an = _mm256_set1_ps(-a),
			bn = _mm256_set1_ps(-b);
	while(n--)
		*to++ = exp256_ps(_mm256_add_ps(_mm256_mul_ps(an, *r1++), _mm256_mul_ps(bn, *r2++)));
}

static void precomputeValues(Context &c, const float a, const float b)
{
	abs(mm(c.r1), c.x1, Batches);
	abs(mm(c.r2), c.x2, Batches);
	computeWaveFunc(mm(c.wv1), mm(c.r1), mm(c.r2), Batches, a, b);
	computeWaveFunc(mm(c.wv2), mm(c.r1), mm(c.r2), Batches, b, a);
}

#endif //QM_MONTE_CARLO_PRECOMPUTE_H
