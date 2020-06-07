//
// Created by xgallom on 4/21/20.
//

#include "constants.h"
#include "Configuration.h"
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
	SkipsMask = Skips - 1u;

int main()
{
	std::random_device randomDevice;
	std::mt19937 generator(randomDevice());

	std::uniform_real_distribution<double>
			normalDist(0.0, 1.0),
			spaceDist(-SpaceDistMax, SpaceDistMax),
			deltaDist(-DeltaDistMax, DeltaDistMax);

	std::fstream file("out.txt", std::ios::out);

	if(!file.is_open())// || !debug.is_open())
		return 1;


	size_t n = 0;
	for(double S = SeparationMin; S <= SeparationMax; S += SeparationStep) {
		std::stringstream ss;
		ss << "debug_" << n << ".txt";
		std::fstream debug(ss.str(), std::ios::out);

		Parameters parameters = generateParameters(S);

		double averageEnergy = 0.;
		for(size_t tries = 0; tries < Tries; ++tries) {
			Context context = generateContext(generator, spaceDist, parameters, S);

			double totalEnergy = 0.;
			size_t count = 0;

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
						<< ((Q2 / ProtonOffest + energy) / ElectronCharge) << " "

						<< context.distanceContext.electronDistance << "\n";

					++count;
				}
			}

			const double result = 2. * Q2 / S + totalEnergy / count;
			averageEnergy += result;
		}

		std::cout 
			<< "Separation: " << (S / BohrRadius) << "\n"
			<< "Result: " << (averageEnergy / Tries / ElectronCharge) << "\n\n";

		file
			<< (S / BohrRadius) << " " 
			<< (averageEnergy / Tries / ElectronCharge) << "\n";

		file.flush();

		++n;
	}

	file.close();

	return 0;
}
