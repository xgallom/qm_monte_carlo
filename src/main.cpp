#include <fstream>
#include <mutex>
#include <thread>
#include "random.h"
#include "config.h"
#include "matrix.h"
#include "energyForParameters.h"

namespace
{
	void init()
	{
		Random::init();
	}

	void threadHandler(double a, double b, size_t t, std::mutex &mutex, VectorD &energies)
	{
		const auto energy = energyForParameters(a, b);

		std::lock_guard<std::mutex> lock(mutex);
		energies[t] = energy;
	}
}

int main()
{
	init();

	Matrix<Config::TrialsBeta, Config::TrialsAlpha> A, B, E = {};

	for(size_t n = 0; n < A.size(); ++n) {
		const auto c = A(n);

		B[n] = Config::Beta1 + Config::dBeta * c.row;
		A[n] = Config::Alpha1 + Config::dAlpha * c.col;
	}

	auto minEnergy = 0., minA = 0., minB = 0.;

	VectorD energies(Config::ThreadCount);
	std::mutex mutex;

	for(size_t n = 0; n < A.size(); n += Config::ThreadCount) {
		std::cerr << "A: " << A[n] << " B: " << B[n] << "\n";
		std::cerr.flush();

		Vector<std::thread> threads;

		for(size_t t = 0; t < Config::ThreadCount; ++t)
			threads.emplace_back(threadHandler, A[n + t], B[n + t], t, std::ref(mutex), std::ref(energies));

		for(auto &thread : threads)
			thread.join();

		for(size_t t = 0; t < Config::ThreadCount; ++t) {
			const auto avgEnergy = energies[t];

			if(avgEnergy < minEnergy) {
				minEnergy = avgEnergy;
				minA = A[n + t];
				minB = B[n + t];
			}

			E[n + t] = avgEnergy;
		}
	}

	std::ofstream f("out.txt");
	for(size_t n = 0; n < A.size(); ++n)
		f << A[n] << " " << B[n] << " " << E[n] << "\n";
	f.close();

	std::cerr << "\nMinimums: " << minEnergy << " " << minA << " " << minB << "\n";

	return 0;
}