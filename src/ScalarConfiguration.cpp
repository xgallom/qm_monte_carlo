//
// Created by xgallom on 4/21/20.
//

#include "ScalarConfiguration.h"

double distance(const Vector &vector)
{
	return sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

double distance(const Vector &vector, const Vector &origin)
{
	double square = 0.;

	for(size_t n = 0; n < Dimension::Count; ++n) {
		const double dimension = vector.data[n] - origin.data[n];
		square += dimension * dimension;
	}

	return sqrt(square);
}

ScalarConfiguration distance(const Configuration &configuration, const Vector &origin)
{
	return {
			{
					.e1 = distance(configuration.e1, origin),
					.e2 = distance(configuration.e2, origin)
			}
	};
}

ElectronProtonConfiguration distances(const Configuration &configuration, const ProtonConfiguration &origin)
{
	return {
			{
					.p1 = distance(configuration, origin.p1),
					.p2 = distance(configuration, origin.p2)
			}
	};
}

double waveFunction(const double &distance, const double &a)
{
	return exp(-distance / a);
}

ScalarConfiguration waveFunction(const ScalarConfiguration &distanceConfiguration, const double &a)
{
	return {{
					.e1 = waveFunction(distanceConfiguration.e1, a),
					.e2 = waveFunction(distanceConfiguration.e2, a)
			}};
}

ElectronProtonConfiguration waveFunction(const ElectronProtonConfiguration &distanceConfiguration, const double &a)
{
	return {{
					.p1 = waveFunction(distanceConfiguration.p1, a),
					.p2 = waveFunction(distanceConfiguration.p2, a)
			}};
}

ExtendedScalarConfiguration waveFunction(
		const ElectronProtonConfiguration &electronProtonConfiguration,
		const double &electronDistance, const double &alpha, const double &beta)
{
	return {{
					.e1 = electronProtonConfiguration.p1.e1 + electronProtonConfiguration.p2.e1,
					.e2 = electronProtonConfiguration.p1.e2 + electronProtonConfiguration.p2.e2,
					.f = exp(electronDistance / (alpha * (1. + beta * electronDistance))) / 10.
			}};
}

