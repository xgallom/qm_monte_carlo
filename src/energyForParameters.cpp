//
// Created by xgallom on 4/24/19.
//

#include "energyForParameters.h"
#include "random.h"
#include "config.h"
#include <memory>

template<typename T>
using owner = std::unique_ptr<T>;

struct alignas(16) D {
	double value;

	static constexpr size_t count()
	{ return 1; }
};

struct D3 {
	D value[3];

	static constexpr size_t count()
	{ return 3; }
};

template<typename T, size_t N>
struct Array {
	T data[N];

	static constexpr size_t count()
	{ return N * T::count(); }
};

template<size_t N>
using ArrayD = Array<D, N>;

template<size_t N>
using ArrayD3 = Array<D3, N>;

using WArrayD = ArrayD<Config::WalkerCount>;
using WArrayD3 = ArrayD3<Config::WalkerCount>;

struct Context {
	D3
			*x1,
			*x2;
	D
			*r1,
			*r2,
			*wv1,
			*wv2,
			*prob,
			*energy;
};

double energyForParameters(double a, double b)
{
	owner<WArrayD3> x1, x2;
	owner<WArrayD> r1, r2, wv1, wv2, prob, energy;

	owner<WArrayD3> xn1, xn2;
	owner<WArrayD> rn1, rn2, wvn1, wvn2, probn;

	Context context = {
			x1->data,
			x2->data,
			r1->data,
			r2->data,
			wv1->data,
			wv2->data,
			prob->data,
			energy->data
	};

	return a * b;
}
