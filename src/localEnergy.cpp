//
// Created by xgallom on 4/21/20.
//

#include <iostream>
#include "localEnergy.h"
#include "constants.h"
#include "transformations.h"

static inline double momentumLaplacianPhi(const double &a)
{
	return 2. / (a * a);
}

static inline double momentumGradPhiGradF(const Context &context, Electron::Enum electron)
{
	const auto &electronConfiguration = context.electronConfiguration.data[electron];

	const double
			dotLeft = projection({ // r_nL
										 {
												 .x = electronConfiguration[Dimension::X] -
													  context.protonConfiguration.p1.x,
												 .y = electronConfiguration[Dimension::Y] -
													  context.protonConfiguration.p1.y,
												 .z = electronConfiguration[Dimension::Z] -
													  context.protonConfiguration.p1.z
										 }
								 }, { // r_12
										 {
												 .x = context.electronConfiguration.e1.x -
													  context.electronConfiguration.e2.x,
												 .y = context.electronConfiguration.e1.y -
													  context.electronConfiguration.e2.y,
												 .z = context.electronConfiguration.e1.z -
													  context.electronConfiguration.e2.z
										 }
								 },
								 context.distanceContext.distances.p1.data[electron]
	),
			dotRight = projection({ // r_nR
										  {
												  .x = electronConfiguration[Dimension::X] -
													   context.protonConfiguration.p2.x,
												  .y = electronConfiguration[Dimension::Y] -
													   context.protonConfiguration.p2.y,
												  .z = electronConfiguration[Dimension::Z] -
													   context.protonConfiguration.p2.z
										  }
								  }, { // r_12
										  {
												  .x = context.electronConfiguration.e1.x -
													   context.electronConfiguration.e2.x,
												  .y = context.electronConfiguration.e1.y -
													   context.electronConfiguration.e2.y,
												  .z = context.electronConfiguration.e1.z -
													   context.electronConfiguration.e2.z
										  }
								  },
								  context.distanceContext.distances.p2.data[electron]
	),

			left = dotLeft * context.waveContext.electronProtonWaves.p1.data[electron],
			right = dotRight * context.waveContext.electronProtonWaves.p2.data[electron];

	return (left + right) / context.waveContext.electronWaves.data[electron];
}

static inline double momentumGradPhiGradF(const Context &context, const double &a, const double &alpha,
										  const double &beta, const double &gamma, const double &gamma2)
{
	const double
			c1 = 2. / (a * alpha),
			c2 = 1. / (context.distanceContext.electronDistance * gamma) - beta / gamma2,

			electron1 = momentumGradPhiGradF(context, Electron::E1),
			electron2 = momentumGradPhiGradF(context, Electron::E2);

	return c1 * c2 * (electron2 - electron1);
}

static inline double momentumLaplacianF(const double &electronDistance, const double &alpha,
										const double &beta, const double &beta2,
										const double &gamma, const double &gamma2,
										const double &gamma3, const double &gamma4)
{
	const double
			c1 = 4. / alpha,
			v1 = 1. / (electronDistance * gamma)
				 - 2. * beta / gamma2
				 + beta2 * electronDistance / gamma3,

			c2 = 2. / (alpha * alpha),
			v2 = 1. / gamma2
				 - 2. * beta * electronDistance / gamma3
				 + beta2 * electronDistance * electronDistance / gamma4;

	return c1 * v1 + c2 * v2;
}

static inline double momentum(const Context &context, const Parameters &parameters)
{
	const double
			electronDistance = context.distanceContext.electronDistance,

			a = parameters.a,
			alpha = parameters.alpha,

			beta = parameters.beta,
			beta2 = beta * beta,

			gamma = 1 + beta * electronDistance,
			gamma2 = gamma * gamma,
			gamma3 = gamma * gamma2,
			gamma4 = gamma * gamma3;

	return -HBM * (
			momentumLaplacianPhi(a)
			+ momentumGradPhiGradF(context, a, alpha, beta, gamma, gamma2)
			+ momentumLaplacianF(electronDistance, alpha, beta, beta2, gamma, gamma2, gamma3, gamma4)
	);
}

static inline double potential(const DistanceContext &distanceContext)
{
	return Q2 * (1. / distanceContext.electronDistance - sumOfReciprocals(distanceContext.distances));
}

double localEnergy(const Context &context, const Parameters &parameters)
{
	const double
			m = momentum(context, parameters),
			p = potential(context.distanceContext),
			energy = m + p;

	return energy;
}
