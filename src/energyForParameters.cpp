//
// Created by xgallom on 4/24/19.
//

#include "energyForParameters.h"
#include "random.h"
#include "config.h"
#include "context.h"
#include <cmath>
#include <iostream>

static double sqr(const double x)
{ return x * x; }

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

static void step(double *to, const double *from, size_t n)
{
	while(n--)
		*to++ = *from++ + Random::get() * Config::dR;
}

static void step(Context &cn, const Context &c)
{
	step(cn.x1.x, c.x1.x, WArrayD::Count());
	step(cn.x1.y, c.x1.y, WArrayD::Count());
	step(cn.x1.z, c.x1.z, WArrayD::Count());
	step(cn.x2.x, c.x2.x, WArrayD::Count());
	step(cn.x2.y, c.x2.y, WArrayD::Count());
	step(cn.x2.z, c.x2.z, WArrayD::Count());
}

static double abs(const double x, const double y, const double z)
{ return sqrt(sqr(x) + sqr(y) + sqr(z)); }

static void abs(double *to, CD3 from, size_t n)
{
	while(n--) {
		const auto
				&x = *from.x++,
				&y = *from.y++,
				&z = *from.z++;

		*to++ = abs(x, y, z);
	}
}

static void abs(double *to, CD3 fromr, CD3 froml, size_t n)
{
	while(n--) {
		const auto
				x = *fromr.x++ - *froml.x++,
				y = *fromr.y++ - *froml.y++,
				z = *fromr.z++ - *froml.z++;

		*to++ = abs(x, y, z);
	}
}

static void wave(double *to, const double *r1, const double *r2, size_t n, const double a, const double b)
{
	while(n--)
		*to++ = std::exp(-a * *r1++ - b * *r2++);
}

static void update(Context &c, const double a, const double b)
{
	abs(c.r1, c.x1, WArrayD::Count());
	abs(c.r2, c.x2, WArrayD::Count());
	wave(c.wv1, c.r1, c.r2, WArrayD::Count(), a, b);
	wave(c.wv2, c.r1, c.r2, WArrayD::Count(), b, a);
}

static void prob(double *to, const double *wv1, const double *wv2, size_t n)
{
	while(n--)
		*to++ = sqr(*wv1++ + *wv2++);
}

static int *cmp(int *updates, const double *p, const double *pn, size_t n)
{
	int i = 0;
	while(n--) {
		const auto result = *pn++ / *p++;

		if(result >= Random::norm())
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

		prob(cn.prob, cn.wv1, cn.wv2, WArrayD::Count());

		int updates[Config::WalkerCount];
		const auto size = static_cast<size_t>(cmp(updates, c.prob, cn.prob, WArrayD::Count()) - updates);

		transfer(c, cn, updates, size);
	}
}

static void addEnergy(double *energy, const Context &context, const double a, const double b)
{
	WArrayD rar;
	double *ra = rar.data;
	abs(ra, context.x2, context.x1, WArrayD::Count());

	const auto
			*r1 = context.r1,
			*r2 = context.r2,
			*wv1 = context.wv1,
			*wv2 = context.wv2;

	const auto e1 = -0.5 * (sqr(a) + sqr(b));

	size_t n = WArrayD::Count();
	while(n--) {
		const auto
				V = 2. / *r1 + 2. / *r2 - 1. / *ra,
				wv = *wv1 + *wv2,
				e2 = (
							 (a * *wv1 + b * *wv2) / *r1 +
							 (b * *wv1 + a * *wv2) / *r2
					 ) / wv;

		*energy += e1 + e2 - V;

		++energy;
		++r1;
		++r2;
		++ra;
		++wv1;
		++wv2;
	}
}

static double accumulate(const double *energy, size_t n)
{
	double sum = 0.;
	while(n--)
		sum += *energy++;

	return sum;
}

double energyForParameters(const double a, const double b)
{
	WArrayD
			x1x, x1y, x1z,
			x2x, x2y, x2z,
			r1, r2,
			wv1, wv2,
			probability, energy,

			xn1x, xn1y, xn1z,
			xn2x, xn2y, xn2z,
			rn1, rn2,
			wvn1, wvn2,
			probn;

	Context
			context = {
			{x1x.data, x1y.data, x1z.data},
			{x2x.data, x2y.data, x2z.data},
			r1.data,
			r2.data,
			wv1.data,
			wv2.data,
			probability.data
	},
			newContext = {
			{xn1x.data, xn1y.data, xn1z.data},
			{xn2x.data, xn2y.data, xn2z.data},
			rn1.data,
			rn2.data,
			wvn1.data,
			wvn2.data,
			probn.data
	};

	// Prepare
	for(auto &e : energy.data)
		e = 0;

	generate(context);
	update(context, a, b);
	prob(context.prob, context.wv1, context.wv2, WArrayD::Count());

	// Thermalize
	run(context, newContext, Config::Therm, a, b);

	// Simulate
	for(size_t n = 0; n < Config::PointsCount; ++n) {
		addEnergy(energy.data, context, a, b);

		run(context, newContext, Config::SkipSteps, a, b);
	}

	// Average
	return accumulate(energy.data, WArrayD::Count()) / Config::EnergiesCount;
}
