#include "config.h"
#include "runParallel.h"
#include "runSingle.h"

#include <fstream>
#include <iostream>

static void prepare(float *a, float *b)
{
	for(size_t n = 0; n < Config::Trials; ++n) {
		*b++ = Config::Beta1 + Config::dBeta * (n / Config::TrialsAlpha);
		*a++ = Config::Alpha1 + Config::dAlpha * (n % Config::TrialsAlpha);
	}
}

static void outputData(const float *A, const float *B, const float *E)
{
	size_t min = 0;

	std::ofstream f("out.txt");
	for(size_t n = 0; n < Config::Trials; ++n) {
		f << A[n] << " " << B[n] << " " << E[n] << "\n";

		if(E[n] < E[min])
			min = n;
	}
	f.close();

	std::cout << "Minimum E = " << std::scientific << E[min] << std::defaultfloat
			  << " at A = " << A[min] << ", B = " << B[min] << "\n";
}

int main(int argc, char *[])
{
	float A[Config::Trials], B[Config::Trials], E[Config::Trials];
	prepare(A, B);

	if(argc > 2)
		runSingle(&A[Config::Trials - 1], &B[Config::Trials - 1], E, 1);
	else if(argc > 1)
		runSingle(A, B, E, Config::Trials);
	else
		runParallel(A, B, E);

	outputData(A, B, E);

	return 0;
}