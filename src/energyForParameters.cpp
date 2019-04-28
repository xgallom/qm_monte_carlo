//
// Created by xgallom on 4/24/19.
//

#include "energyForParameters.h"
#include "config.h"
#include "context.h"

#include <cmath>

static float sqr(const float x)
{ return x * x; }

static void generate(float *to, size_t n, Random &random)
{
	while(n--)
		*to++ = random.dist();
}

static void generate(Context &c, Random &random)
{
	generate(c.x1.x, WArrayD::Count(), random);
	generate(c.x1.y, WArrayD::Count(), random);
	generate(c.x1.z, WArrayD::Count(), random);
	generate(c.x2.x, WArrayD::Count(), random);
	generate(c.x2.y, WArrayD::Count(), random);
	generate(c.x2.z, WArrayD::Count(), random);
}

static void step(float *to, const float *from, size_t n, Random &random)
{
	while(n--)
		*to++ = *from++ + random.get() * Config::dR;
}

static void step(Context &cn, const Context &c, Random &random)
{
	step(cn.x1.x, c.x1.x, WArrayD::Count(), random);
	step(cn.x1.y, c.x1.y, WArrayD::Count(), random);
	step(cn.x1.z, c.x1.z, WArrayD::Count(), random);
	step(cn.x2.x, c.x2.x, WArrayD::Count(), random);
	step(cn.x2.y, c.x2.y, WArrayD::Count(), random);
	step(cn.x2.z, c.x2.z, WArrayD::Count(), random);
}

static float abs(const float x, const float y, const float z)
{ return sqrtf(sqr(x) + sqr(y) + sqr(z)); }

static void abs(float *to, CD3 from, size_t n)
{
	while(n--) {
		const auto
				&x = *from.x++,
				&y = *from.y++,
				&z = *from.z++;

		*to++ = abs(x, y, z);
	}
}

static void abs(float *to, CD3 fromr, CD3 froml, size_t n)
{
	while(n--) {
		const auto
				x = *fromr.x++ - *froml.x++,
				y = *fromr.y++ - *froml.y++,
				z = *fromr.z++ - *froml.z++;

		*to++ = abs(x, y, z);
	}
}

static void wave(float *to, const float *r1, const float *r2, size_t n, const float a, const float b)
{
	while(n--)
		*to++ = std::exp(-a * *r1++ - b * *r2++);
}

static void update(Context &c, const float a, const float b)
{
	abs(c.r1, c.x1, WArrayD::Count());
	abs(c.r2, c.x2, WArrayD::Count());
	wave(c.wv1, c.r1, c.r2, WArrayD::Count(), a, b);
	wave(c.wv2, c.r1, c.r2, WArrayD::Count(), b, a);
}

static void prob(float *to, const float *wv1, const float *wv2, size_t n)
{
	while(n--)
		*to++ = sqr(*wv1++ + *wv2++);
}

static int *cmp(int *updates, const float *p, const float *pn, size_t n, Random &random)
{
	int i = 0;
	while(n--) {
		const auto result = *pn++ / *p++;

		if(result >= random.norm())
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

static void run(Context &c, Context &cn, size_t n, const float a, const float b, Random &random)
{
	while(n--) {
		step(cn, c, random);
		update(cn, a, b);

		prob(cn.prob, cn.wv1, cn.wv2, WArrayD::Count());

		int updates[Config::WalkerCount];
		const auto size = static_cast<size_t>(cmp(updates, c.prob, cn.prob, WArrayD::Count(), random) - updates);

		transfer(c, cn, updates, size);
	}
}

static void addEnergy(float *energy, const Context &context, const float a, const float b)
{
	WArrayD rar;
	float *ra = rar.data;
	abs(ra, context.x2, context.x1, WArrayD::Count());

	const auto
			*r1 = context.r1,
			*r2 = context.r2,
			*wv1 = context.wv1,
			*wv2 = context.wv2;

	const auto e1 = -0.5f * (sqr(a) + sqr(b));

	size_t n = WArrayD::Count();
	while(n--) {
		const auto
				V = 2.f / *r1 + 2.f / *r2 - 1.f / *ra,
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

static float accumulate(const float *energy, size_t n)
{
	float sum = 0.;
	while(n--)
		sum += *energy++;

	return sum;
}

float energyForParameters(const float a, const float b, Random &random)
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

	generate(context, random);
	update(context, a, b);
	prob(context.prob, context.wv1, context.wv2, WArrayD::Count());

	// Thermalize
	run(context, newContext, Config::Therm, a, b, random);

	// Simulate
	for(size_t n = 0; n < Config::PointsCount; ++n) {
		addEnergy(energy.data, context, a, b);

		run(context, newContext, Config::SkipSteps, a, b, random);
	}

	// Average
	return accumulate(energy.data, WArrayD::Count()) / Config::EnergiesCount;
}
