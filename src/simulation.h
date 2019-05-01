//
// Created by xgallom on 4/29/19.
//

#ifndef QM_MONTE_CARLO_SIMULATION_H
#define QM_MONTE_CARLO_SIMULATION_H

#include "context.h"

static void computeProbabilities(Mm *to, const Mm *wv1, const Mm *wv2, size_t n)
{
	while(n--) {
		const auto wv = _mm256_add_ps(*wv1++, *wv2++);
		*to++ = _mm256_mul_ps(wv, wv);
	}
}

static size_t *compareProbabilities(size_t *updatesBuffer, const Mm *p, const Mm *pn, size_t n, Random &random)
{
	size_t i = 0;
	while(n--) {
		const auto
				norm = _mm256_set_ps(random.norm(), random.norm(), random.norm(), random.norm(),
									 random.norm(), random.norm(), random.norm(), random.norm()),

				compare = _mm256_cmp_ps(_mm256_div_ps(*pn++, *p++), norm, 13 /* GE */);

		for(size_t m = 0; m < PerBatch; ++m) {
			if(compare[m] != 0)
				*updatesBuffer++ = i;

			++i;
		}
	}

	return updatesBuffer;
}

static void replaceUpdated(Context &to, const Context &from, size_t *updatesBuffer, size_t n)
{
	while(n--) {
		const auto i = *updatesBuffer++;
		to.x1.x[i] = from.x1.x[i];
		to.x1.y[i] = from.x1.y[i];
		to.x1.z[i] = from.x1.z[i];
		to.x2.x[i] = from.x2.x[i];
		to.x2.y[i] = from.x2.y[i];
		to.x2.z[i] = from.x2.z[i];
		to.r1[i] = from.r1[i];
		to.r2[i] = from.r2[i];
		to.wv1[i] = from.wv1[i];
		to.wv2[i] = from.wv2[i];
		to.prob[i] = from.prob[i];
	}
}

static void simulateSteps(Context &c, Context &cn, size_t n, const float a, const float b, Random &random)
{
	while(n--) {
		stepPositions(cn, c, random);
		precomputeValues(cn, a, b);

		computeProbabilities(mm(cn.prob), mm(cn.wv1), mm(cn.wv2), Batches);

		size_t updatesBuffer[Config::WalkerCount];
		const auto *const updatesEnd = compareProbabilities(updatesBuffer, mm(c.prob), mm(cn.prob), Batches, random);

		replaceUpdated(c, cn, updatesBuffer, static_cast<size_t>(updatesEnd - updatesBuffer));
	}
}

#endif //QM_MONTE_CARLO_SIMULATION_H
