#include "config.h"
#include "energyForParameters.h"

#include <fstream>
#include <iostream>

#include <mutex>
#include <thread>

#include <vector>

namespace
{
	void init()
	{
	}

	void prepare(float *a, float *b)
	{
		for(size_t n = 0; n < Config::Trials; ++n) {
			*b++ = Config::Beta1 + Config::dBeta * (n / Config::TrialsAlpha);
			*a++ = Config::Alpha1 + Config::dAlpha * (n % Config::TrialsAlpha);
		}
	}

	void threadHandler(const float *a, const float *b, float *e, size_t t,
					   std::mutex &mutex, size_t &progress)
	{
		float
				aBuffer[Config::TasksPerThread],
				bBuffer[Config::TasksPerThread],
				eBuffer[Config::TasksPerThread];

		Random random;

		mutex.lock();
		for(size_t n = 0; n < t; ++n) {
			aBuffer[n] = a[n];
			bBuffer[n] = b[n];
		}

		std::cout << "Starting thread-local simulation\n";
		mutex.unlock();

		auto *eB = eBuffer;
		const auto *aB = aBuffer, *bB = bBuffer;
		while(t--)
			*eB++ = energyForParameters(*aB++, *bB++, random);

		mutex.lock();
		std::cout << "Finishing thread-local simulation\n";
		for(const auto *eW = eBuffer; eW != eB;)
			*e++ = *eW++;

		progress -= Config::TasksPerThread;
		mutex.unlock();
	}

	void runSingle(const float *a, const float *b, float *e, size_t t)
	{
		Random random;

		while(t--)
			*e++ = energyForParameters(*a++, *b++, random);
	}

	void spawnThreads(std::vector<std::thread> &threads,
					  const float *a, const float *b, float *e,
					  std::mutex &mutex, size_t &progress)
	{
		std::cout
				<< "Spawning " << Config::ThreadCount << " threads, each with " << Config::TasksPerThread
				<< " out of " << Config::Trials << " tasks\n";

		for(size_t t = 0; t < Config::ThreadCount; ++t) {
			const auto tasksStart = t * Config::TasksPerThread;
			const auto tasksEnd = tasksStart + Config::TasksPerThread - 1;

			std::cout << "Thread " << t << " for tasks: " << tasksStart << "-" << tasksEnd << "\n";

			threads.emplace_back(threadHandler,
								 a + tasksStart,
								 b + tasksStart,
								 e + tasksStart,
								 Config::TasksPerThread,
								 std::ref(mutex),
								 std::ref(progress)
			);
		}
	}

	void waitForFinish(std::mutex &mutex, const size_t &progress)
	{
		size_t missing = Config::Trials, missingOld = 0;

		while(missing) {
			mutex.lock();
			missing = progress;
			mutex.unlock();

			if(missing != missingOld) {
				const auto done = Config::Trials - missing;
				std::cout
						<< "Progress: " << done << "/" << Config::Trials
						<< " (" << 100 * done / Config::Trials << " %)\n";
			}

			missingOld = missing;

			std::this_thread::sleep_for(std::chrono::milliseconds(5000));
		}
	}

	void join(std::vector<std::thread> &threads)
	{
		for(auto &thread : threads)
			thread.join();
	}

	size_t outputData(const float *A, const float *B, const float *E)
	{
		size_t min = 0;

		std::ofstream f("out.txt");
		for(size_t n = 0; n < Config::Trials; ++n) {
			f << A[n] << " " << B[n] << " " << E[n] << "\n";

			if(E[n] < E[min])
				min = n;
		}
		f.close();

		return min;
	}
}

int main()
{
	float A[Config::Trials], B[Config::Trials], E[Config::Trials];

	init();
	prepare(A, B);

	std::vector<std::thread> threads;
	std::mutex mutex;
	size_t progress = Config::Trials;

	//runSingle(A, B, E, Config::Trials);
	spawnThreads(threads, A, B, E, mutex, progress);

	//waitForFinish(mutex, progress);
	join(threads);

	size_t min = outputData(A, B, E);

	std::cout << "Minimum E = " << E[min] << " at A = " << A[min] << ", B = " << B[min] << "\n";

	return 0;
}