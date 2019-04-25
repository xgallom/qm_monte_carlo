//
// Created by xgallom on 4/25/19.
//

#ifndef QM_MONTE_CARLO_CONTEXT_H
#define QM_MONTE_CARLO_CONTEXT_H

#include <xmmintrin.h>

static constexpr size_t Align = 16;
static constexpr size_t Double = 64;
static constexpr size_t Batch = 128;
static constexpr size_t DPerBatch = Batch / Double;

struct alignas(Align) D {
	double value;

	static constexpr size_t Count()
	{ return 1; }

	static constexpr size_t Batches()
	{ return Count() / DPerBatch; }
};

struct D3 {
	double
			*x,
			*y,
			*z;
};

template<typename T>
inline D *d(T *t)
{ return reinterpret_cast<D *>(t); }

template<typename T>
inline const D *d(const T *t)
{ return reinterpret_cast<const D *>(t); }

using Mmd = __m128d;
struct Mmd3 {
	Mmd
			*x,
			*y,
			*z;
};

struct CMmd3 {
	const Mmd
			*x,
			*y,
			*z;
};

template<typename T>
inline Mmd *mmd(T *t)
{ return reinterpret_cast<Mmd *>(t); }

template<typename T>
inline const Mmd *mmd(const T *t)
{ return reinterpret_cast<const Mmd *>(t); }

template<typename T>
inline Mmd3 mmd3(T t)
{
	return {
			reinterpret_cast<Mmd *>(t.x),
			reinterpret_cast<Mmd *>(t.y),
			reinterpret_cast<Mmd *>(t.z)
	};
}

template<typename T>
inline CMmd3 cmmd3(T t)
{
	return {
			reinterpret_cast<const Mmd *>(t.x),
			reinterpret_cast<const Mmd *>(t.y),
			reinterpret_cast<const Mmd *>(t.z)
	};
}

template<typename T, size_t N>
struct alignas(Align) Array {
	T data[N];

	static constexpr size_t Count()
	{ return N; }

	static constexpr size_t Batches()
	{ return Count() / DPerBatch; }
};

template<size_t N>
using ArrayD = Array<double, N>;
using WArrayD = ArrayD<Config::WalkerCount>;
using WArrayD = ArrayD<Config::WalkerCount>;

template<typename T>
double *dbl(T *t)
{ return reinterpret_cast<double *>(t); }

template<typename T>
const double *dbl(const T *t)
{ return reinterpret_cast<const double *>(t); }

struct Context {
	D3
			x1,
			x2;
	double
			*r1,
			*r2,
			*wv1,
			*wv2,
			*prob;
};

#endif //QM_MONTE_CARLO_CONTEXT_H
