//
// Created by xgallom on 4/24/19.
//

#include "energyForParameters.h"
#include "random.h"
#include "config.h"
#include "context.h"
#include <cmath>
#include <iostream>

static void generate(double *to, size_t n)
{
	while(n--)
		*to++ = Random::dist();
}

static void generate(Context &c)
{
	generate(c.x1.x, WArrayD::Count());
	generate(c.x1.y, WArrayD::Count());
	generate(c.x1.z, WArrayD::Count());
	generate(c.x2.x, WArrayD::Count());
	generate(c.x2.y, WArrayD::Count());
	generate(c.x2.z, WArrayD::Count());
}

static void step(Mmd *to, const Mmd *from, size_t n)
{
	while(n--)
		*to++ =
				_mm_add_pd(
						*from++,
						Mmd{Random::get() * Config::dR, Random::get() * Config::dR}
				);
}

static void step(Context &cn, const Context &c)
{
	step(mmd(cn.x1.x), mmd(c.x1.x), WArrayD::Batches());
	step(mmd(cn.x1.y), mmd(c.x1.y), WArrayD::Batches());
	step(mmd(cn.x1.z), mmd(c.x1.z), WArrayD::Batches());
	step(mmd(cn.x2.x), mmd(c.x2.x), WArrayD::Batches());
	step(mmd(cn.x2.y), mmd(c.x2.y), WArrayD::Batches());
	step(mmd(cn.x2.z), mmd(c.x2.z), WArrayD::Batches());
}

static void abs(Mmd *to, CMmd3 from, size_t n)
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

static void abs(Mmd *to, CMmd3 fromr, CMmd3 froml, size_t n)
{
	while(n--) {
		const auto
				x = _mm_sub_pd(*fromr.x++, *froml.x++),
				y = _mm_sub_pd(*fromr.y++, *froml.y++),
				z = _mm_sub_pd(*fromr.z++, *froml.z++);

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

static void wave(double *to, const Mmd *r1, const Mmd *r2, size_t n, const double an, const double bn)
{
	const Mmd
			An = {an, an},
			Bn = {bn, bn};

	while(n--) {
		const auto
				l = _mm_mul_pd(An, *r1++),
				r = _mm_mul_pd(Bn, *r2++),

				x = _mm_add_pd(l, r);

		//std::cout << (r1 - 1)[0][0] << " " << (r2 - 1)[0][0] << " " << x[0] << "\n";
		//std::cout << (r1 - 1)[0][1] << " " << (r2 - 1)[0][1] << " " << x[1] << "\n";

		*to++ = std::exp(x[0]);
		*to++ = std::exp(x[1]);
	}

	//std::cout.flush();
}

static void update(Context &c, const double a, const double b)
{
	abs(mmd(c.r1), cmmd3(c.x1), WArrayD::Batches());
	abs(mmd(c.r2), cmmd3(c.x2), WArrayD::Batches());
	wave(c.wv1, mmd(c.r1), mmd(c.r2), WArrayD::Batches(), -a, -b);
	wave(c.wv2, mmd(c.r1), mmd(c.r2), WArrayD::Batches(), -b, -a);
}

static void prob(Mmd *to, const Mmd *wv1, const Mmd *wv2, size_t n)
{
	while(n--) {
		const auto wv = _mm_add_pd(*wv1++, *wv2++);
		*to++ = _mm_mul_pd(wv, wv);
	}
}

static int *cmp(int *updates, const Mmd *p, const Mmd *pn, size_t n)
{
	int i = 0;
	while(n--) {
		const auto result =
				_mm_cmpge_pd(
						_mm_div_pd(*pn++, *p++),
						Mmd{Random::norm(), Random::norm()}
				);

		if(result[0] != 0)
			*updates++ = i;
		++i;

		if(result[1] != 0)
			*updates++ = i;
		++i;
	}

	return updates;
}

static void transfer(Context &to, const Context &from, int *updates, size_t n)
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

static void run(Context &c, Context &cn, size_t n, const double a, const double b)
{
	while(n--) {
		step(cn, c);
		update(cn, a, b);

		prob(mmd(cn.prob), mmd(cn.wv1), mmd(cn.wv2), WArrayD::Batches());

		int updates[Config::WalkerCount];
		const auto size =
				static_cast<size_t>(
						cmp(updates, mmd(c.prob), mmd(cn.prob), WArrayD::Batches()) - updates
				);
		transfer(c, cn, updates, size);
	}
}

static void addEnergy(Mmd *energy, const Context &context, const double a, const double b)
{
	const Mmd
			A = {a, a},
			B = {b, b},
			C1 = {2, 2},
			C2 = {-1, -1},
			C3 = {-0.5, -0.5};

	WArrayD rar;
	Mmd *ra = mmd(&rar);
	abs(ra, cmmd3(context.x2), cmmd3(context.x1), WArrayD::Batches());

	const auto
			*r1 = mmd(context.r1),
			*r2 = mmd(context.r2),
			*wv1 = mmd(context.wv1),
			*wv2 = mmd(context.wv2);

	size_t m = 0;
	size_t n = WArrayD::Batches();
	while(n--) {
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

		std::cout
				<< context.x1.x[m] << " " << context.x1.y[m] << " " << context.x1.z[m] << " "
				<< context.x2.x[m] << " " << context.x2.y[m] << " " << context.x2.z[m] << " "

				<< r1[0][0] << " " << r2[0][0] << " " << rar.data[m] << " " << wv[0] << " "
				<< wv1[0][0] << " " << wv2[0][0] << " "
				<< V[0] << " " << e1[0] << " " << e2[0] << " " << energy[0][0] << "\n";
		m += 1;
		std::cout
				<< context.x1.x[m] << " " << context.x1.y[m] << " " << context.x1.z[m] << " "
				<< context.x2.x[m] << " " << context.x2.y[m] << " " << context.x2.z[m] << " "

				<< r1[0][1] << " " << r2[0][1] << " " << rar.data[m] << " " << wv[1] << " "
				<< wv1[0][1] << " " << wv2[0][1] << " "
				<< V[1] << " " << e1[1] << " " << e2[1] << " " << energy[0][1] << "\n";
		m += 1;
	}
}

static double accumulate(const Mmd *energy, size_t n)
{
	Mmd sum = {0, 0};
	while(n--)
		sum = _mm_add_pd(sum, *energy++);

	return sum[0] + sum[1];
}

double energyForParameters(const double a, const double b)
{
	WArrayD
			x1x,
			x1y,
			x1z,
			x2x,
			x2y,
			x2z,
			r1,
			r2,
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

	addEnergy(mmd(&energy), context, a, b);

	std::cout.flush();

	return 0.0;
	run(context, newContext, Config::Therm, a, b);

	for(size_t n = 0; n < Config::PointsCount; ++n) {
		addEnergy(mmd(&energy), context, a, b);

		run(context, newContext, Config::SkipSteps, a, b);
	}

	return accumulate(mmd(&energy), WArrayD::Batches()) / Config::EnergiesCount;
}
