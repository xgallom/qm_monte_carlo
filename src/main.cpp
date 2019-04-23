#include <fstream>
#include "vector.h"
#include "random.h"
#include "config.h"
#include "matrix.h"

namespace
{
	void init()
	{
		//std::cout << std::scientific;
		Random::init();
	}

	double prob(const Vector<Vector3D> &r, double a, double b)
	{
		const auto
				r1 = abs(r[0]),
				r2 = abs(r[1]);

		return sqr(exp(-a * r1 - b * r2) + exp(-b * r1 - a * r2));
	}

	double energy(const Vector<Vector3D> &r, double a, double b)
	{
		const auto
				ra = abs(r[1] - r[0]),
				r1 = abs(r[0]),
				r2 = abs(r[1]);

		const auto V = 2. / r1 + 2. / r2 - 1. / ra;
		const auto wave = exp(-a * r1 - b * r2) + exp(-b * r1 - a * r2);

		return -0.5 * (sqr(a) + sqr(b)) + (
												  (a * exp(-a * r1 - b * r2) + b * exp(-b * r1 - a * r2)) / r1 +
												  (b * exp(-a * r1 - b * r2) + a * exp(-b * r1 - a * r2)) / r2
										  ) / wave - V;
	}

	Vector<Vector3D> step(const Vector<Vector3D> &r)
	{
		Vector<Vector3D> result = {};
		result.reserve(r.size());

		for(const auto &rn : r)
			result.push_back(rn + Random::vectorOffset() * Config::dR);

		return result;
	}
}

int main()
{
	init();

	Matrix<Config::TrialsBeta, Config::TrialsAlpha> A, B, E = {}, R0 = {}, R1 = {};

	for(size_t n = 0; n < A.size(); ++n) {
		const auto c = A(n);

		B[n] = Config::Beta1 + Config::dBeta * c.row;
		A[n] = Config::Alpha1 + Config::dAlpha * c.col;
	}

	std::cout << A << "\n\n";
	std::cout << B << "\n\n";

	double minEnergy = 0., minA, minB;

	for(size_t n = 0; n < A.size(); ++n) {
		const auto
				&a = A[n],
				&b = B[n];

		auto allEnergy = 0.;
		auto
				allR0 = 0.,
				allR1 = 0.;

		std::cout << "A: " << a << " B: " << b << "\n";
		for(size_t w = 0; w < Config::WalkerCount; ++w) {
			Vector<Vector3D> r = {Random::vector(), Random::vector()};

			auto walkerEnergy = 0.;
			auto
					walkerR0 = 0.,
					walkerR1 = 0.;

			for(size_t s = 0; s < Config::Steps; ++s) {
				auto rNew = step(r);

				const auto p = prob(rNew, a, b) / prob(r, a, b);

				if(p >= Random::get())
					r = rNew;

				if(s > Config::Therm && (s - Config::Therm) % Config::ThermStep == 0) {
					walkerEnergy += energy(r, a, b);
					walkerR0 += accumulate(r[0]);
					walkerR1 += accumulate(r[1]);
				}
			}

			allEnergy += walkerEnergy / Config::PointsCount;
			allR0 += walkerR0 / Config::PointsCount;
			allR1 += walkerR1 / Config::PointsCount;
		}

		const auto
				avgEnergy = allEnergy / Config::WalkerCount,
				avgR0 = allR0 / Config::WalkerCount,
				avgR1 = allR1 / Config::WalkerCount;
		std::cout << "Avg: " << avgEnergy << "\n";

		std::cout << "Avg R1: " << avgR0 << "\n";
		std::cout << "Avg R2: " << avgR1 << "\n";

		if(avgEnergy < minEnergy) {
			minEnergy = avgEnergy;
			minA = a;
			minB = b;
		}

		E[n] = avgEnergy;
		R0[n] = avgR0;
		R1[n] = avgR1;
	}

	std::ofstream f("out.txt");
	for(int n = 0; n < A.size(); ++n)
		f << A[n] << " " << B[n] << " " << E[n] << " " << R0[n] << " " << R1[n] << "\n";
	f.close();

	std::cout << "Minimums: " << minEnergy << " " << minA << " " << minB << "\n";

	return 0;
}