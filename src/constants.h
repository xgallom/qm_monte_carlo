//
// Created by xgallom on 4/21/20.
//

#ifndef QM_MONTE_CARLO_SRC_CONSTANTS_H
#define QM_MONTE_CARLO_SRC_CONSTANTS_H

static constexpr double
		PlanckConstant = 6.62607004e-34, // m2 kg s-1
		HBar = 1.0545718e-34, // m2 kg s-1
		ElectronMass = 9.10938356e-31,   // kg
		ElectronCharge = 1.60217662e-19, // C
		BohrRadius = 5.29177210903e-11, // m
		ProtonDistance = 1.2 * BohrRadius,
		ProtonOffest = ProtonDistance / 2.,

		ParameterAlpha = 2. * BohrRadius, // m
		ParameterBeta = 0.25 / BohrRadius, // m
		ParameterA = 3.71294e-11, // m
		CoulombConstant = 8.9875517923e9, // kg m3 s-2 C-2

		HBM = HBar * HBar / (2. * ElectronMass), // kg m4 s-2
		Q2 = CoulombConstant * ElectronCharge * ElectronCharge // kg m3 s-2
;

#endif //QM_MONTE_CARLO_SRC_CONSTANTS_H
