//
// Created by xgallom on 4/24/19.
//

#include "energyForParameters.h"
#include "mean.h"
#include "translations.h"
#include "precompute.h"
#include "simulation.h"
#include "localEnergy.h"

float energyForParameters(const float a, const float b, Random &random)
{
	WArrayD
			x1x, x1y, x1z,
			x2x, x2y, x2z,
			r1, r2,
			wv1, wv2,
			probability, energy,

			xn1x, xn1y, xn1z,
			xn2x, xn2y, xn2z,
			rn1, rn2,
			wvn1, wvn2,
			probn;

	Context
			context = {
			{x1x.data, x1y.data, x1z.data},
			{x2x.data, x2y.data, x2z.data},
			r1.data,
			r2.data,
			wv1.data,
			wv2.data,
			probability.data
	},
			newContext = {
			{xn1x.data, xn1y.data, xn1z.data},
			{xn2x.data, xn2y.data, xn2z.data},
			rn1.data,
			rn2.data,
			wvn1.data,
			wvn2.data,
			probn.data
	};

	// Prepare
	auto *e = mm(energy.data);
	clear(e, Batches);

	// Initialize
	generateRandomPositions(context, random);
	precomputeValues(context, a, b);
	computeProbabilities(mm(context.prob), mm(context.wv1), mm(context.wv2), Batches);

	// Thermalize
	simulateSteps(context, newContext, Config::Therm, a, b, random);

	// Simulate
	for(size_t n = 0; n < Config::PointsCount; ++n) {
		computeLocalEnergies(e, context, a, b);

		// Skip
		simulateSteps(context, newContext, Config::SkipSteps, a, b, random);
	}

	// Average
	return mean(e, Batches, Config::EnergiesCount);
}
