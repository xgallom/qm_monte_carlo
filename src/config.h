//
// Created by xgallom on 4/17/19.
//

#ifndef QM_MONTE_CARLO_CONFIG_H
#define QM_MONTE_CARLO_CONFIG_H

namespace Config
{
	constexpr size_t
			WalkerCount = 1,
			Steps = 10000000,
			SkipSteps = 500,
			Therm = 100000,
			TotalSteps = Steps + Therm,

			TrialsAlpha = 10,
			TrialsBeta = 12,

			PointsCount = Steps / SkipSteps,

			EnergiesCount = PointsCount * WalkerCount,

			ThreadCount = 1;

	constexpr double
			dR = 1.,
			Alpha1 = 0.9,
			dAlpha = 0.05,
			Beta1 = 1.9,
			dBeta = 0.05,

			Dist = 2.;
}

#endif //QM_MONTE_CARLO_CONFIG_H
