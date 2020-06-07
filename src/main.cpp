//
// Created by xgallom on 4/21/20.
//

#include "constants.h"
#include "Context.h"
#include "Parameters.h"
#include "localEnergy.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <random>

constexpr double
		SpaceDistMax = 2. * BohrRadius,
		DeltaDistMax = 0.05,

		SeparationMin = 0.3 * BohrRadius,
		SeparationMax = 6.3 * BohrRadius,
		SeparationStep = 0.05 * BohrRadius;

constexpr size_t
		Tries = 1u << 10u,
		Steps = 1u << 13u,
		Skips = 1u << 8u,
		SkipsMask = Skips - 1u,
		Computations = Tries * Steps / Skips,
		SeparationCount = (SeparationMax - SeparationMin) / SeparationStep;

static inline double simulateTry(
		std::fstream &debug,
		std::uniform_real_distribution<double> &normalDist,
		std::uniform_real_distribution<double> &spaceDist,
		std::uniform_real_distribution<double> &deltaDist,
		std::mt19937 &generator,
		const Parameters &parameters,
		double S)
{
	Context context = generateContext(generator, spaceDist, parameters, S);

	double totalEnergy = 0.;

	for(size_t step = 0; step < Steps; ++step) {
		Context newContext = generateContext(context, generator, deltaDist, parameters, S);

		if(newContext.waveContext.waveSquared / context.waveContext.waveSquared >= normalDist(generator))
			context = newContext;

		if(!(step & SkipsMask)) {
			double energy = localEnergy(context, parameters);

			totalEnergy += energy;

			debug
					<< (context.electronConfiguration.e1.x / BohrRadius) << " "
					<< (context.electronConfiguration.e1.y / BohrRadius) << " "
					<< (context.electronConfiguration.e1.z / BohrRadius) << " "

					<< (context.electronConfiguration.e2.x / BohrRadius) << " "
					<< (context.electronConfiguration.e2.y / BohrRadius) << " "
					<< (context.electronConfiguration.e2.z / BohrRadius) << " "

					<< context.waveContext.electronWaves.e1 << " "
					<< context.waveContext.electronWaves.e2 << " "
					<< context.waveContext.electronWaves.f << " "

					<< context.waveContext.wave << " "
					<< context.waveContext.waveSquared << " "
					<< ((2. * Q2 / S + energy) / ElectronCharge) << " "

					<< context.distanceContext.electronDistance << "\n";
		}
	}

	return totalEnergy;
}

static inline void simulateSeparation(
		std::fstream &file,
		std::uniform_real_distribution<double> &normalDist,
		std::uniform_real_distribution<double> &spaceDist,
		std::uniform_real_distribution<double> &deltaDist,
		std::mt19937 &generator,
		size_t separation)
{
	std::stringstream ss;
	ss << "debug_" << separation << ".txt";
	std::fstream debug(ss.str(), std::ios::out);

	const auto S = SeparationMin + separation * SeparationStep;

	Parameters parameters = generateParameters(S);

	double totalEnergy = 0.;

	for(size_t tries = 0; tries < Tries; ++tries)
		totalEnergy += simulateTry(
				debug,
				normalDist,
				spaceDist,
				deltaDist,
				generator,
				parameters,
				S
		);

	const auto
			averageEnergy = (2. * Q2 / S + totalEnergy / Computations) / ElectronCharge,
			S_a0 = S / BohrRadius;

	std::cout
			<< "Separation: " << S_a0 << " a0\n"
			<< "Result: " << averageEnergy << " eV\n\n";

	file
			<< S_a0 << " "
			<< averageEnergy << "\n";

	file.flush();
}

int main()
{
	std::random_device randomDevice;
	std::mt19937 generator(randomDevice());

	std::uniform_real_distribution<double>
			normalDist(0.0, 1.0),
			spaceDist(-SpaceDistMax, SpaceDistMax),
			deltaDist(-DeltaDistMax, DeltaDistMax);

	std::fstream file("out.txt", std::ios::out);

	if(!file.is_open())
		return 1;

	for(size_t separation = 0; separation <= SeparationCount; ++separation)
		simulateSeparation(
				file,
				normalDist,
				spaceDist,
				deltaDist,
				generator,
				separation
		);

	file.close();

	return 0;
}
