//
// Created by xgallom on 4/21/20.
//

#include "Configuration.h"
//#include "constants.h"

//static constexpr double Range = 2., Delta = Range / 100.;

Configuration generateConfiguration(
		std::mt19937 &generator,
		std::uniform_real_distribution<double> &distribution,
		double protonDistance)
{
//	return {{
//					.e1 = {.data = {0., 0., -Range * BohrRadius}},
//					.e2 = {.data = {0., 0., -Range * BohrRadius}},
//	}};
	return {{
					.e1 = {.data = {distribution(generator), distribution(generator),
					 distribution(generator) - protonDistance / 2.}},
					.e2 = {.data = {distribution(generator), distribution(generator),
					 distribution(generator) + protonDistance / 2.}}
			}};
}

Configuration generateConfiguration(
		const Configuration &oldConfiguration,
		std::mt19937 &generator,
		std::uniform_real_distribution<double> &distribution,
		double protonDistance)
{
//	Configuration result = {};
//
//	result.e1.z = oldConfiguration.e1.z + Delta * BohrRadius;
//	result.e2.z = oldConfiguration.e2.z;
//
//	if(result.e1.z >= Range * BohrRadius) {
//		result.e1.z = -Range * BohrRadius;
//		result.e2.z += Delta * BohrRadius;
//	}
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
