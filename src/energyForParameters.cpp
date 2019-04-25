//
// Created by xgallom on 4/24/19.
//

#include "energyForParameters.h"
#include "random.h"
#include "config.h"
#include "context.h"
#include <memory>
#include <iostream>
#include <cmath>
#include <numeric>

static size_t Out1 = 0, Out2 = 0;

template<typename T>
using owner = std::unique_ptr<T>;

void generate(double *to, size_t n)
{
	while(n--)
		*to++ = Random::dist();
}

void generate(Context &c)
{
	generate(c.x1.x, WArrayD::Count());
	generate(c.x1.y, WArrayD::Count());
	generate(c.x1.z, WArrayD::Count());
	generate(c.x2.x, WArrayD::Count());
	generate(c.x2.y, WArrayD::Count());
	generate(c.x2.z, WArrayD::Count());
}

void step(double *to, double *from, size_t n)
{
	while(n--)
		*to++ = *from++ + Random::get() * Config::dR;
}

void step(Context &cn, const Context &c)
{
	step(cn.x1.x, c.x1.x, WArrayD::Count());
	step(cn.x1.y, c.x1.y, WArrayD::Count());
	step(cn.x1.z, c.x1.z, WArrayD::Count());
	step(cn.x2.x, c.x2.x, WArrayD::Count());
	step(cn.x2.y, c.x2.y, WArrayD::Count());
	step(cn.x2.z, c.x2.z, WArrayD::Count());
}

void abs(Mmd *to, CMmd3 from, size_t n)
{
	while(n--) {
		const auto
				&x = *from.x++,
				&y = *from.y++,
				&z = *from.z++;

		*to++ =
				_mm_sqrt_pd(
						_mm_add_pd(
								_mm_add_pd(
										_mm_mul_pd(x, x),
										_mm_mul_pd(y, y)),
								_mm_mul_pd(z, z)
						)
				);
	}
}

void abs(Mmd *to, CMmd3 fromr, CMmd3 froml, size_t n)
{
	while(n--) {
		const auto
				x = *fromr.x++ - *froml.x++,
				y = *fromr.y++ - *froml.y++,
				z = *fromr.z++ - *froml.z++;

		*to++ =
				_mm_sqrt_pd(
						_mm_add_pd(
								_mm_add_pd(
										_mm_mul_pd(x, x),
										_mm_mul_pd(y, y)),
								_mm_mul_pd(z, z)
						)
				);
	}
}

void wave(double *to, const Mmd *r1, const Mmd *r2, size_t n, const double an, const double bn)
{
	const auto
			An = _mm_load_pd1(&an),
			Bn = _mm_load_pd1(&bn);

	while(n--) {
		const auto x =
				_mm_add_pd(
						_mm_mul_pd(An, *r1++),
						_mm_mul_pd(Bn, *r2++)
				);

		*to++ = exp(x[0]);
		*to++ = exp(x[1]);
	}
}

void update(Context &c, const double a, const double b)
{
	abs(mmd(c.r1), cmmd3(c.x1), WArrayD::Batches());
	abs(mmd(c.r2), cmmd3(c.x2), WArrayD::Batches());
	wave(c.wv1, mmd(c.r1), mmd(c.r2), WArrayD::Batches(), -a, -b);
	wave(c.wv2, mmd(c.r1), mmd(c.r2), WArrayD::Batches(), -b, -a);
}

void prob(Mmd *to, const Mmd *wv1, const Mmd *wv2, size_t n)
{
	while(n--) {
		const auto wv = _mm_add_pd(*wv1++, *wv2++);
		*to++ = _mm_mul_pd(wv, wv);
	}
}

int *cmp(int *updates, const Mmd *p, const Mmd *pn, size_t n)
{
	int i = 0;
	while(n--) {
		const auto result =
				_mm_cmpge_pd(
						_mm_div_pd(*pn++, *p++),
						Mmd{Random::norm(), Random::norm()}
				);
		if(result[0])
			*updates++ = i;
		++i;

		if(result[1])
			*updates++ = i;
		++i;
	}

	return updates;
}

void transfer(Context &to, const Context &from, int *updates, size_t n)
{
	while(n--) {
		const auto i = *updates++;
		to.x1.x[i] = from.x1.x[i];
		to.x1.y[i] = from.x1.y[i];
		to.x1.z[i] = from.x1.z[i];
		to.x2.x[i] = from.x2.x[i];
		to.x2.y[i] = from.x2.y[i];
		to.x2.z[i] = from.x2.z[i];
		to.r1[i] = from.r1[i];
		to.r2[i] = from.r2[i];
		to.wv1[i] = from.wv1[i];
		to.wv2[i] = from.wv2[i];
		to.prob[i] = from.prob[i];
	}
}

void run(Context &c, Context &cn, size_t n, const double a, const double b)
{
	while(n--) {
		step(cn, c);
		update(cn, a, b);

		prob(mmd(cn.prob), mmd(cn.wv1), mmd(cn.wv2), WArrayD::Batches());

		int updates[Config::WalkerCount];
		const auto size = static_cast<size_t>(cmp(updates, mmd(c.prob), mmd(cn.prob), WArrayD::Batches()) - updates);
		transfer(c, cn, updates, size);
	}
}

void addEnergy(Mmd *energy, const Context &context, const double a, const double b)
{
	const auto
			A = _mm_load_pd1(&a),
			B = _mm_load_pd1(&b),
			C1 = Mmd{2, 2},
			C2 = Mmd{-1, -1},
			C3 = Mmd{-0.5, -0.5};

	WArrayD rar;
	Mmd *ra = mmd(&rar);
	abs(ra, cmmd3(context.x2), cmmd3(context.x1), WArrayD::Batches());

	const auto
			*r1 = mmd(context.r1),
			*r2 = mmd(context.r2),
			*wv1 = mmd(context.wv1),
			*wv2 = mmd(context.wv2);

	for(size_t n = 0; n < Config::WalkerCount; ++n) {
		const auto V =
				_mm_add_pd(
						_mm_div_pd(C1, *r1),
						_mm_sub_pd(
								_mm_div_pd(C1, *r2),
								_mm_div_pd(C2, *ra)
						)
				),
				wv = _mm_add_pd(*wv1, *wv2),

				e1 =
				_mm_mul_pd(
						C3,
						_mm_add_pd(
								_mm_mul_pd(A, A),
								_mm_mul_pd(B, B)
						)
				),

				e2 =
				_mm_div_pd(
						_mm_add_pd(
								_mm_div_pd(
										_mm_add_pd(
												_mm_mul_pd(A, *wv1),
												_mm_mul_pd(B, *wv2)
										),
										*r1
								),
								_mm_div_pd(
										_mm_add_pd(
												_mm_mul_pd(B, *wv1),
												_mm_mul_pd(A, *wv2)
										),
										*r2
								)
						),
						wv
				);

		*energy =
				_mm_add_pd(
						*energy,
						_mm_add_pd(
								e1,
								_mm_sub_pd(e2, V)
						)
				);

		++energy;
		++r1;
		++r2;
		++ra;
		++wv1;
		++wv2;
	}
}

double accumulate(const Mmd *energy, size_t n)
{
	Mmd sum = {0, 0};
	while(n--)
		sum = _mm_add_pd(sum, *energy++);

	return sum[0] + sum[1];
}

double energyForParameters(const double a, const double b)
{
	Out1 = Out2 = 0;

	WArrayD
			x1x,
			x1y,
			x1z,
			x2x,
			x2y,
			x2z,
			r1,
			r2 ,
			wv1,
			wv2,
			probability,
			energy;

	for(auto &e : energy.data)
		e = 0;

	WArrayD
			xn1x,
			xn1y,
			xn1z,
			xn2x,
			xn2y,
			xn2z,
			rn1,
			rn2,
			wvn1,
			wvn2,
			probn;

	Context context = {
			{dbl(&x1x), dbl(&x1y), dbl(&x1z)},
			{dbl(&x2x), dbl(&x2y), dbl(&x2z)},
			dbl(&r1),
			dbl(&r2),
			dbl(&wv1),
			dbl(&wv2),
			dbl(&probability)
	};

	Context newContext = {
			{dbl(&xn1x), dbl(&xn1y), dbl(&xn1z)},
			{dbl(&xn2x), dbl(&xn2y), dbl(&xn2z)},
			dbl(&rn1),
			dbl(&rn2),
			dbl(&wvn1),
			dbl(&wvn2),
			dbl(&probn)
	};

	generate(context);
	update(context, a, b);
	prob(mmd(context.prob), mmd(context.wv1), mmd(context.wv2), WArrayD::Batches());

	run(context, newContext, Config::Therm, a, b);

	for(size_t n = 0; n < Config::PointsCount; ++n) {
		addEnergy(mmd(&energy), context, a, b);

		run(context, newContext, Config::SkipSteps, a, b);
	}

	std::cout << Config::PointsCount << " " << Config::WalkerCount << "\n";

	return accumulate(mmd(&energy), Config::WalkerCount) / Config::EnergiesCount;
}
