//
// Created by xgallom on 4/21/20.
//

#include "Context.h"
#include "transformations.h"
#include "constants.h"

DistanceContext generateContext(
		const Configuration &configuration,
		const ProtonConfiguration &protonConfiguration)
{
	return {
			.distances = distances(configuration, protonConfiguration),
			.electronDistance = distance(configuration.e1, configuration.e2)
	};
}

WaveContext generateContext(
		const DistanceContext &distanceContext,
		const Parameters &parameters)
{
	WaveContext result = {
			.electronProtonWaves = waveFunction(distanceContext.distances, parameters.a),
			.electronWaves = {},
			.wave = {},
			.waveSquared = {}
	};

	result.electronWaves = waveFunction(
			result.electronProtonWaves,
			distanceContext.electronDistance,
			parameters.alpha,
			parameters.beta
	);

	result.wave = sum(result.electronWaves);
	result.waveSquared = square(result.wave);

	return result;
}

Context generateContext(
		std::mt19937 &generator,
		std::uniform_real_distribution<double> &distribution,
		const Parameters &parameters,
		double S)
{
	Configuration configuration = generateConfiguration(generator, distribution, S);
	ProtonConfiguration protonConfiguration = generateConfiguration(S / 2.);

	DistanceContext distanceContext = generateContext(configuration, protonConfiguration);

	return {
			.electronConfiguration = configuration,
			.protonConfiguration = protonConfiguration,
			.distanceContext = distanceContext,
			.waveContext = generateContext(distanceContext, parameters)
	};
}

Context generateContext(
		const Context &oldContext,
		std::mt19937 &generator,
		std::uniform_real_distribution<double> &distribution,
		const Parameters &parameters,
		double S)
{
	Configuration configuration = generateConfiguration(
			oldContext.electronConfiguration,
			generator,
			distribution,
			S
	);
	ProtonConfiguration protonConfiguration = generateConfiguration(S / 2.);

	DistanceContext distanceContext = generateContext(configuration, protonConfiguration);

	return {
			.electronConfiguration = configuration,
			.protonConfiguration = protonConfiguration,
			.distanceContext = distanceContext,
			.waveContext = generateContext(distanceContext, parameters)
	};
}
