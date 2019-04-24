//
// Created by xgallom on 4/17/19.
//

#ifndef QM_MONTE_CARLO_CONFIG_H
#define QM_MONTE_CARLO_CONFIG_H

namespace Config
{
	constexpr size_t
			WalkerCount = 300,
			Steps = 60000,
			ThermStep = 1000,
			Therm = 10000,

			TrialsAlpha = 40,

			PointsSize = (Steps - Therm) / ThermStep;

	constexpr double
			dR = 10.,
			Alpha1 = 0.1,
			dAlpha = 0.05,

			Dist = 0.5;
}

#endif //QM_MONTE_CARLO_CONFIG_H
