//
// Created by xgallom on 4/25/19.
//

#ifndef QM_MONTE_CARLO_CONTEXT_H
#define QM_MONTE_CARLO_CONTEXT_H

#include "config.h"
#include "avx_mathfun.h"

using Mm = __m256;

static Mm *mm(float *x)
{ return reinterpret_cast<Mm *>(x); }

static const Mm *mm(const float *x)
{ return reinterpret_cast<const Mm *>(x); }

static constexpr size_t
		Align = 32,
		BatchSize = 256,
		PerBatch = BatchSize / (8 * sizeof(float)),
		Batches = Config::WalkerCount / PerBatch;


struct CMm3 {
	const Mm
			*x,
			*y,
			*z;
};

struct Mm3 {
	Mm
			*x,
			*y,
			*z;

	operator CMm3() const
	{ return {x, y, z}; }
};

struct CD3 {
	const float
			*x,
			*y,
			*z;

	operator CMm3() const
	{ return {mm(x), mm(y), mm(z)}; }
};

struct D3 {
	float
			*x,
			*y,
			*z;

	operator CD3() const
	{ return {x, y, z}; }

	operator CMm3() const
	{ return {mm(x), mm(y), mm(z)}; }

	operator Mm3() const
	{ return {mm(x), mm(y), mm(z)}; }
};

template<typename T, size_t N>
struct alignas(Align) Array {
	T data[N];

	static constexpr size_t Count()
	{ return N; }
};

template<size_t N>
using ArrayD = Array<float, N>;
using WArrayD = ArrayD<Config::WalkerCount>;
using WArrayD = ArrayD<Config::WalkerCount>;

struct Context {
	D3
			x1,
			x2;
	float
			*r1,
			*r2,
			*wv1,
			*wv2,
			*prob;
};

#endif //QM_MONTE_CARLO_CONTEXT_H
