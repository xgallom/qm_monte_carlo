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

	double prob(double r, double a)
	{
		return exp(-2. * a * sqr(r));
	}

	double energy(double r, double a)
	{
		return a + sqr(r) * (0.5 - 2. * sqr(a));
	}

	double step(double r)
	{
		return std::min(10., std::max(-10., r + Random::get() * Config::dR));
	}
}

int main()
{
	init();

	VectorD A = {}, E = {};
	A.reserve(Config::TrialsAlpha);
	E.reserve(Config::TrialsAlpha);

	for(size_t n = 0; n < Config::TrialsAlpha; ++n)
		A.push_back(Config::Alpha1 + Config::dAlpha * n);

	std::cout << A << "\n\n";

	double minEnergy = 0., minA = 0.;

	for(const auto &a : A) {
		auto allEnergy = 0.;

		std::cout << "a: " << a << "\n";

		for(size_t w = 0; w < Config::WalkerCount; ++w) {
			auto r = Random::dist();

			VectorD points = {};
			points.reserve(Config::PointsSize);

			for(size_t s = 0; s < Config::Steps; ++s) {
				const auto rNew = step(r);

				const auto p = prob(rNew, a) / prob(r, a);

				if(p >= Random::norm())
					r = rNew;

				if(s > Config::Therm && (s - Config::Therm) % Config::ThermStep == 0)
					points.push_back(r);
			}

			allEnergy += accumulate(transform(points, [&a](const double &r) { return energy(r, a); })) / points.size();
		}

		const auto avgEnergy = allEnergy / Config::WalkerCount;
		std::cout << "Avg: " << avgEnergy << "\n";

		if(avgEnergy < minEnergy) {
			minEnergy = avgEnergy;
			minA = a;
		}

		E.push_back(avgEnergy);
	}

	std::ofstream f("out.txt");
	for(int n = 0; n < A.size(); ++n)
		f << A[n] << " " << E[n] << "\n";
	f.close();

	std::cout << "Minimums: " << minEnergy << " " << minA << "\n";

	return 0;
}