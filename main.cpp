#include <iostream>
#include <cmath>
#include <random>
#include <cstring>
#include <fstream>

static constexpr double
		PlanckConstant = 6.62607004e-34, // m2 kg s-1
		HBar = 1.0545718e-34, // m2 kg s-1
		ElectronMass = 9.10938356e-31,   // kg
		ElectronCharge = 1.60217662e-19, // C
		BohrRadius = 5.29177210903e-11, // m
		ProtonDistance = 1.2 * BohrRadius,
		ProtonOffest = ProtonDistance / 2.,

		ParameterAlpha = 2. * BohrRadius, // m
		ParameterA = 3.71294e-11, // m
		CoulombConstant = 8.9875517923e9, // kg m3 s-2 C-2

		HBM = HBar * HBar / (2. * ElectronMass), // kg m4 s-2
		e2 = CoulombConstant * ElectronCharge * ElectronCharge // kg m3 s-2
;

enum Dim {
	X = 0,
	Y,
	Z,

	Dims
};

enum Ele {
	E1 = 0,
	E2,

	Eles
};

enum Prot {
	P1 = 0,
	P2,

	Prots
};

double evaluateA(size_t iterations, double S)
{
	double a = 0.;

	while(iterations--) {
		const double
				ex = exp(-S / a),
				paren = 1. + ex,
				bohrParen = BohrRadius / paren;

		a = a - (bohrParen - a) / (bohrParen * S * ex / paren - 1.);
	}

	return a;
}

double waveFunction(double r, double a)
{
	return exp(-r / a);
}

double cross(double r12, double alpha, double beta)
{
	return exp(r12 / (alpha * (1 + beta * r12)));
}

double distance(double x, double y, double z)
{
	return sqrt(x * x + y * y + z * z);
}

double momentum(double a, double r, double rs)
{
	return r * r / a / rs / rs / rs - 3. / a / rs + r * r / a / a / rs / rs;
}

double momentumCross(double alpha, double beta, double r12)
{
	const double
			gamma = 1. + beta * r12,
			gamma2 = gamma * gamma,
			gamma3 = gamma2 * gamma,
			gamma4 = gamma3 * gamma,
			beta2 = beta * beta;

	return 2. / alpha / alpha *
		   (
				   1. / gamma2
				   - 2. * beta * r12 / gamma3
				   + beta2 * r12 * r12 / gamma4
		   ) +
		   4. / alpha *
		   (
				   beta2 * r12 / gamma3
				   - beta / gamma2
		   );
}

double localEnergy(
		const double (&dist)[Eles][Prots],
		const double (&posDist)[Eles],
		const double (&wv)[Eles][Prots],
		const double (&wvs)[Eles + 1],
		double eleDist,
		double f,
		double a,
		double alpha,
		double beta)
{
	const double
			p11 = momentum(a, posDist[E1], dist[E1][P1]) * wv[E1][P1],
			p12 = momentum(a, posDist[E1], dist[E1][P2]) * wv[E1][P2],
			p21 = momentum(a, posDist[E2], dist[E2][P1]) * wv[E2][P1],
			p22 = momentum(a, posDist[E2], dist[E2][P2]) * wv[E2][P2],

			p1 = (p11 + p12) / wvs[E1],
			p2 = (p21 + p22) / wvs[E2],
			pf = momentumCross(alpha, beta, eleDist),

			p = -HBM * (p1 + p2 + pf),
			v = e2 * (1. / eleDist - 1. / dist[E1][P1] - 1. / dist[E1][P2] - 1. / dist[E2][P1] - 1. / dist[E2][P2]);

	return p + v;
}

int main()
{
	constexpr double
			a0 = BohrRadius,
			alpha = ParameterAlpha;

	std::random_device device;
	std::mt19937 rd(device());
	std::uniform_real_distribution<double> nd(0, 1.0);

	std::fstream file("out.txt", std::ios::out);
	std::uniform_real_distribution<double> dd(-a0 / 2., a0 / 2.);

	if(!file.is_open())
		return 1;

	for(size_t w = 0; w < 30; ++w) {
		const double
				S = ProtonOffest * (double(w + 5) / 5.0),
				a = evaluateA(64, S);
		double lowestEnergy = INFINITY, lowestBeta = 0.;
		std::cout << "[" << (w + 1) << "/30] S = " << S << "\n";

		for(size_t b = 0; b < 100; ++b) {
			const double beta = 0.025 * double(b + 1);
			std::cout << "[" << (w + 1) << "/30][" << (b + 1) << "/100] beta = " << beta << "\n";

			double pos[Eles][Dims] = {
					{dd(rd), dd(rd), dd(rd)},
					{dd(rd), dd(rd), dd(rd)}
			};

			double dist[Eles][Prots] = {
					{
							distance(pos[E1][X], pos[E1][Y], pos[E1][Z] + S),
							distance(pos[E1][X], pos[E1][Y], pos[E1][Z] - S)
					},
					{
							distance(pos[E2][X], pos[E2][Y], pos[E2][Z] + S),
							distance(pos[E2][X], pos[E2][Y], pos[E2][Z] - S)
					}
			};

			double posDist[Eles] = {
					distance(pos[E1][X], pos[E1][Y], pos[E1][Z]),
					distance(pos[E2][X], pos[E2][Y], pos[E2][Z])
			};

			double eleDist = distance(
					pos[E1][X] - pos[E2][X],
					pos[E1][Y] - pos[E2][Y],
					pos[E1][Z] - pos[E2][Z]
			);

			double wv[Eles][Prots] = {
					{waveFunction(dist[E1][P1], a), waveFunction(dist[E1][P2], a)},
					{waveFunction(dist[E2][P1], a), waveFunction(dist[E2][P2], a)}
			};

			double wvs[Eles + 1] = {
					wv[E1][P1] + wv[E1][P2],
					wv[E2][P1] + wv[E2][P2]
			};

			double f = cross(eleDist, alpha, beta);

			wvs[Eles] = wvs[E1] * wvs[E2] * f;
			double wvsqr = wvs[Eles] * wvs[Eles];

			double energy = 0;
			size_t count = 0;

			double minX1 = std::numeric_limits<double>::max(), maxX1 = std::numeric_limits<double>::min(),
					minX2 = std::numeric_limits<double>::max(), maxX2 = std::numeric_limits<double>::min();

			bool wasEvaluated = false;
			double evaluatedEnergy = 0.;

			for(size_t step = 0; step < 1u << 21u; ++step) {
				double newPos[Eles][Dims] = {
						{dd(rd) + pos[E1][X], dd(rd) + pos[E1][Y], dd(rd) + pos[E1][Z]},
						{dd(rd) + pos[E2][X], dd(rd) + pos[E2][Y], dd(rd) + pos[E2][Z]}
				};

				double newDist[Eles][Prots] = {
						{
								distance(newPos[E1][X], newPos[E1][Y], newPos[E1][Z] + S),
								distance(newPos[E1][X], newPos[E1][Y], newPos[E1][Z] - S)
						},
						{
								distance(newPos[E2][X], newPos[E2][Y], newPos[E2][Z] + S),
								distance(newPos[E2][X], newPos[E2][Y], newPos[E2][Z] - S)
						}
				};

				double newEleDist = distance(
						newPos[E1][X] - newPos[E2][X],
						newPos[E1][Y] - newPos[E2][Y],
						newPos[E1][Z] - newPos[E2][Z]
				);

				double newWv[Eles][Prots] = {
						{waveFunction(newDist[E1][P1], a), waveFunction(newDist[E1][P2], a)},
						{waveFunction(newDist[E2][P1], a), waveFunction(newDist[E2][P2], a)}
				};

				double newWvs[Eles + 1] = {
						newWv[E1][P1] + newWv[E1][P2],
						newWv[E2][P1] + newWv[E2][P2]
				};

				double newF = cross(newEleDist, alpha, beta);

				newWvs[Eles] = newWvs[E1] * newWvs[E2] * newF;
				double newWvsqr = newWvs[Eles] * newWvs[Eles];

				if(newWvsqr / wvsqr >= nd(rd)) {
					wasEvaluated = false;

					double newPosDist[Eles] = {
							distance(newPos[E1][X], newPos[E1][Y], newPos[E1][Z]),
							distance(newPos[E2][X], newPos[E2][Y], newPos[E2][Z])
					};

					memcpy(pos, newPos, sizeof(double) * Eles * Dims);
					memcpy(dist, newDist, sizeof(double) * Eles * Prots);
					memcpy(wv, newWv, sizeof(double) * Eles * Prots);
					memcpy(wvs, newWvs, sizeof(double) * (Eles + 1));
					memcpy(posDist, newPosDist, sizeof(double) * Eles);
					eleDist = newEleDist;
					f = newF;

					wvsqr = newWvsqr;

					minX1 = std::min(minX1, pos[E1][X]);
					maxX1 = std::max(maxX1, pos[E1][X]);
					minX2 = std::min(minX2, pos[E2][X]);
					maxX2 = std::max(maxX2, pos[E2][X]);
				}

				if(step >= 1u << 16u && !(step & 0x3ffu)) {
					if(!wasEvaluated) {
						evaluatedEnergy = localEnergy(dist, posDist, wv, wvs, eleDist, f, a, alpha, beta);
						wasEvaluated = true;
					}

					energy += evaluatedEnergy;
					++count;
				}
			}

			const double thisEnergy = energy / count;
			if(lowestEnergy > thisEnergy) {
				lowestEnergy = thisEnergy;
				lowestBeta = beta;
			}

			std::cout
					<< "  "
					<< (minX1 / a0)
					<< ":"
					<< (maxX1 / a0)
					<< " "
					<< (minX2 / a0)
					<< ":"
					<< (maxX2 / a0)
					<< "\n";
		}

		file << (2. * S / a0) << " " << lowestBeta << " " << (e2 / S + lowestEnergy) << "\n";
	}

	file.close();

	std::cout << "[-] Done \a\a\n";

	return 0;
}
