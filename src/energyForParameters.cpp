//
// Created by xgallom on 4/24/19.
//

#include "energyForParameters.h"
#include "config.h"
#include "context.h"
#include "avx_mathfun.h"

static constexpr size_t Batches = Config::WalkerCount / PerBatch;

static float sqr(const float x)
{ return x * x; }

static void clear(Mm *to, size_t n)
{
	const auto zero = _mm256_setzero_ps();
	while(n--)
		*to++ = zero;
}

static void generate(Mm *to, size_t n, Random &random)
{
	while(n--)
		*to++ = _mm256_set_ps(
				random.dist(), random.dist(), random.dist(), random.dist(),
				random.dist(), random.dist(), random.dist(), random.dist()
		);
}

static void generate(Context &c, Random &random)
{
	generate(mm(c.x1.x), Batches, random);
	generate(mm(c.x1.y), Batches, random);
	generate(mm(c.x1.z), Batches, random);
	generate(mm(c.x2.x), Batches, random);
	generate(mm(c.x2.y), Batches, random);
	generate(mm(c.x2.z), Batches, random);
}

static void step(Mm *to, const Mm *from, size_t n, Random &random)
{
	while(n--)
		*to++ = _mm256_add_ps(*from++, _mm256_set_ps(
				random.get(), random.get(), random.get(), random.get(),
				random.get(), random.get(), random.get(), random.get()
		));
}

static void step(Context &cn, const Context &c, Random &random)
{
	step(mm(cn.x1.x), mm(c.x1.x), Batches, random);
	step(mm(cn.x1.y), mm(c.x1.y), Batches, random);
	step(mm(cn.x1.z), mm(c.x1.z), Batches, random);
	step(mm(cn.x2.x), mm(c.x2.x), Batches, random);
	step(mm(cn.x2.y), mm(c.x2.y), Batches, random);
	step(mm(cn.x2.z), mm(c.x2.z), Batches, random);
}

static Mm abs(const Mm x, const Mm y, const Mm z)
{
	return _mm256_sqrt_ps(
			_mm256_add_ps(
					_mm256_mul_ps(x, x),
					_mm256_add_ps(
							_mm256_mul_ps(y, y),
							_mm256_mul_ps(z, z)
					)
			)
	);
}

static void abs(Mm *to, CMm3 from, size_t n)
{
	while(n--) {
		const auto
				&x = *from.x++,
				&y = *from.y++,
				&z = *from.z++;

		*to++ = abs(x, y, z);
	}
}

static void abs(Mm *to, CMm3 fromr, CMm3 froml, size_t n)
{
	while(n--) {
		const auto
				x = _mm256_sub_ps(*fromr.x++, *froml.x++),
				y = _mm256_sub_ps(*fromr.y++, *froml.y++),
				z = _mm256_sub_ps(*fromr.z++, *froml.z++);

		*to++ = abs(x, y, z);
	}
}

static void wave(Mm *to, const Mm *r1, const Mm *r2, size_t n, const float a, const float b)
{
	const Mm
			an = _mm256_set1_ps(-a),
			bn = _mm256_set1_ps(-b);
	while(n--)
		*to++ = exp256_ps(_mm256_add_ps(_mm256_mul_ps(an, *r1++), _mm256_mul_ps(bn, *r2++)));
}

static void update(Context &c, const float a, const float b)
{
	abs(mm(c.r1), c.x1, Batches);
	abs(mm(c.r2), c.x2, Batches);
	wave(mm(c.wv1), mm(c.r1), mm(c.r2), Batches, a, b);
	wave(mm(c.wv2), mm(c.r1), mm(c.r2), Batches, b, a);
}

static void prob(Mm *to, const Mm *wv1, const Mm *wv2, size_t n)
{
	while(n--) {
		const auto wv = _mm256_add_ps(*wv1++, *wv2++);
		*to++ = _mm256_mul_ps(wv, wv);
	}
}

static int *cmp(int *updates, const Mm *p, const Mm *pn, size_t n, Random &random)
{
	int i = 0;
	while(n--) {
		const auto norm = _mm256_set_ps(random.norm(), random.norm(), random.norm(), random.norm(),
										random.norm(), random.norm(), random.norm(), random.norm()),

				compare = _mm256_cmp_ps(_mm256_div_ps(*pn++, *p++), norm, 13 /* GE */);

		for(size_t m = 0; m < PerBatch; ++m) {
			if(compare[m] != 0)
				*updates++ = i;

			++i;
		}
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

		prob(mm(cn.prob), mm(cn.wv1), mm(cn.wv2), Batches);

		int updates[Config::WalkerCount];
		const auto size = static_cast<size_t>(cmp(updates, mm(c.prob), mm(cn.prob), Batches, random) - updates);

		transfer(c, cn, updates, size);
	}
}

static void addEnergy(Mm *energy, const Context &context, const float a, const float b)
{
	WArrayD rar;
	abs(mm(rar.data), context.x2, context.x1, Batches);

	const auto
			*r1 = mm(context.r1),
			*r2 = mm(context.r2),
			*wv1 = mm(context.wv1),
			*wv2 = mm(context.wv2),
			*ra = mm(rar.data);

	const auto
			e1 = _mm256_set1_ps(-0.5f * (sqr(a) + sqr(b))),
			C1 = _mm256_set1_ps(2.f),
			C2 = _mm256_set1_ps(-1.f),
			Ca = _mm256_set1_ps(a),
			Cb = _mm256_set1_ps(b);

	size_t n = WArrayD::Count();
	while(n--) {
		const auto
				r1R = _mm256_rcp_ps(*r1),
				r2R = _mm256_rcp_ps(*r2),
				raR = _mm256_rcp_ps(*ra),

				V =
				_mm256_add_ps(
						_mm256_mul_ps(C1, r1R),
						_mm256_add_ps(
								_mm256_mul_ps(C1, r2R),
								_mm256_mul_ps(C2, raR)
						)
				),
				wvR = _mm256_rcp_ps(_mm256_add_ps(*wv1, *wv2)),
				e2 =
				_mm256_mul_ps(
						_mm256_add_ps(
								_mm256_mul_ps(
										_mm256_add_ps(
												_mm256_mul_ps(Ca, *wv1),
												_mm256_mul_ps(Cb, *wv2)
										),
										r1R
								),
								_mm256_mul_ps(
										_mm256_add_ps(
												_mm256_mul_ps(Cb, *wv1),
												_mm256_mul_ps(Ca, *wv2)
										),
										r2R
								)
						),
						wvR
				);

		*energy =
				_mm256_add_ps(
						*energy,
						_mm256_add_ps(
								e1,
								_mm256_sub_ps(e2, V)
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

static float accumulate(const Mm *energy, size_t n)
{
	Mm sum = _mm256_setzero_ps();

	while(n--)
		sum = _mm256_add_ps(sum, *energy++);

	return sum[0] + sum[1] + sum[2] + sum[3] +
		   sum[4] + sum[5] + sum[6] + sum[7];
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

	auto *e = mm(energy.data);

	// Prepare
	clear(e, Batches);

	generate(context, random);
	update(context, a, b);
	prob(mm(context.prob), mm(context.wv1), mm(context.wv2), Batches);

	// Thermalize
	run(context, newContext, Config::Therm, a, b, random);

	// Simulate
	for(size_t n = 0; n < Config::PointsCount; ++n) {
		addEnergy(e, context, a, b);

		run(context, newContext, Config::SkipSteps, a, b, random);
	}

	// Average
	return accumulate(e, Batches) / Config::EnergiesCount;
}
