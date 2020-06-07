//
// Created by xgallom on 4/21/20.
//

#include "Configuration.h"

Configuration generateConfiguration(
		std::mt19937 &generator,
		std::uniform_real_distribution<double> &distribution,
		double protonDistance)
{
	return {{
					.e1 = {
							.data = {
									distribution(generator), distribution(generator),
									distribution(generator) - protonDistance / 2.
							}},
					.e2 = {
							.data = {
									distribution(generator), distribution(generator),
									distribution(generator) + protonDistance / 2.
							}}
			}};
}

Configuration generateConfiguration(
		const Configuration &oldConfiguration,
		std::mt19937 &generator,
		std::uniform_real_distribution<double> &distribution,
		double protonDistance)
{
	Configuration result;

	for(size_t e = 0; e < Proton::Count; ++e) {
		result.data[e][Dimension::X] = oldConfiguration.data[e][Dimension::X] +
									   protonDistance * distribution(generator);
		result.data[e][Dimension::Y] = oldConfiguration.data[e][Dimension::Y] +
									   protonDistance * distribution(generator);
		result.data[e][Dimension::Z] = oldConfiguration.data[e][Dimension::Z] +
									   protonDistance * distribution(generator);
	}

	return result;
}

ProtonConfiguration generateConfiguration(double protonOffset)
{
	return {{
					.p1 = {{.x = 0., .y = 0., .z = -protonOffset}},
					.p2 = {{.x = 0., .y = 0., .z = protonOffset}}
			}};
}
