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
#include <random>

int main()
{
	std::random_device randomDevice;
	std::mt19937 generator(randomDevice());

	std::uniform_real_distribution<double>
			normalDist(0, 1.0),
			spaceDist(-2 * BohrRadius, 2 * BohrRadius),
			deltaDist(-0.05 * BohrRadius, 0.05 * BohrRadius);

	std::fstream
			file("out.txt", std::ios::out),
			debug("debug.txt", std::ios::out);

	if(!file.is_open() || !debug.is_open())
		return 1;


	for(double S = 0.3 * BohrRadius; S < 6.3 * BohrRadius; S += 0.05 * BohrRadius) {

		Parameters parameters = generateParameters(S);

		debug << "Starting simulation\n";
		debug << "parameters: " << parameters << "\n\n";

		double averageEnergy = 0.;

		for(size_t tries = 0; tries < 1u << 10u; ++tries) {
			Context context = generateContext(generator, spaceDist, parameters, S);

			double totalEnergy = 0.;
			size_t count = 0;

			double minimum = 0.;
			Context minimumContext = {};


			for(size_t step = 0; step < (1u << 16u); ++step) {
				Context newContext = generateContext(context, generator, deltaDist, parameters, S);

				if(/*newContext.waveContext.waveSquared < 1.0 &&*/
						newContext.waveContext.waveSquared / context.waveContext.waveSquared >= normalDist(generator)
						)
					context = newContext;

				if(!(step & ((1u << 8u) - 1))) {
					double energy = localEnergy(context, parameters);

					if(energy < minimum) {
						minimum = energy;
						minimumContext = context;
					}

					totalEnergy += energy;

//					debug << "[" << step << "]\n";
//					debug << "context: " << context << "\n";
//					debug << "energy: " << energy << " J = " << (energy / ElectronCharge) << " eV\n";
//					debug << "totalEnergy: " << totalEnergy << " J = " << (totalEnergy / ElectronCharge) << " eV\n\n";
/*
				file
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

						<< context.distanceContext.electronDistance << "\n";*/

					++count;
				}
			}

//			std::cout << "Separation: " << (S / BohrRadius) << "\n";
//			std::cout << "Minimum energy: " << (minimum / ElectronCharge) << "\n";
//		std::cout << "Minimum context: " << minimumContext << "\n\n";

			const double result = 2. * Q2 / S + totalEnergy / count;
//			std::cout << "\nResult: " << result << " J = " << (result / ElectronCharge) << " eV\n\n";

			averageEnergy += result;
		}

		std::cout << "Separation: " << (S / BohrRadius) << "\nResult: " << (averageEnergy / double(1u << 10u) / ElectronCharge) << "\n\n";
		file << (S / BohrRadius) << " " << (averageEnergy / double(1u << 10u) / ElectronCharge) << "\n";
		file.flush();
	}

	file.close();

	return 0;
}
