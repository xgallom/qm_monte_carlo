//
// Created by xgallom on 4/25/19.
//

#ifndef QM_MONTE_CARLO_CONTEXT_H
#define QM_MONTE_CARLO_CONTEXT_H

static constexpr size_t Align = 32;

struct CD3 {
	const float
			*x,
			*y,
			*z;
};

struct D3 {
	float
			*x,
			*y,
			*z;

	operator CD3() const { return {x, y, z}; }
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
