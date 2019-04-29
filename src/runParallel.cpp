//
// Created by xgallom on 4/29/19.
//

#include "runParallel.h"
#include "config.h"
#include "runSingle.h"

#include <mutex>
#include <thread>
#include <vector>
#include <iostream>

static void threadHandler(const float *a, const float *b, float *e, size_t t, std::mutex &mutex)
{
	float eBuffer[Config::TasksPerThread];

	runSingle(a, b, eBuffer, t);

	std::lock_guard<std::mutex> lock(mutex);

	std::cout << "Finishing thread-local simulation\n";

	const auto *eB = eBuffer;
	while(t--)
		*e++ = *eB++;
}

static void spawnThreads(std::vector<std::thread> &threads, const float *a, const float *b, float *e, std::mutex &mutex)
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
							 std::ref(mutex)
		);
	}
}

static void join(std::vector<std::thread> &threads)
{
	for(auto &thread : threads)
		thread.join();
}

void runParallel(const float *a, const float *b, float *e)
{
	std::vector<std::thread> threads;
	std::mutex mutex;

	spawnThreads(threads, a, b, e, mutex);
	join(threads);
}

