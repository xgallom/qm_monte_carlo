//
// Created by xgallom on 4/17/19.
//

#ifndef QM_MONTE_CARLO_CONFIG_H
#define QM_MONTE_CARLO_CONFIG_H

namespace Config
{
	constexpr size_t
			WalkerCount = 150,
			Steps = 6000,
			ThermStep = 1,
			Therm = 1000,

			TrialsAlpha = 40,
			TrialsBeta = 1,

			PointsSize = (Steps - Therm) / ThermStep;

	constexpr double
			dR = 10,
			Alpha1 = 0.1,
			dAlpha = 0.05,
			Beta1 = 1.9,
			dBeta = 0.05,

			Dist = 2.;
}

#endif //QM_MONTE_CARLO_CONFIG_H
