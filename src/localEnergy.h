//
// Created by xgallom on 4/29/19.
//

#ifndef QM_MONTE_CARLO_LOCALENERGY_H
#define QM_MONTE_CARLO_LOCALENERGY_H

#include "arithmetic.h"

static float computeE1(const float a, const float b)
{
	return -0.5f * (sqr(a) + sqr(b));
}

static Mm computePotential(const Mm C1, const Mm C2, const Mm r1R, const Mm r2R, const Mm raR)
{
	return
			_mm256_add_ps(
					_mm256_mul_ps(C1, r1R),
					_mm256_add_ps(
							_mm256_mul_ps(C1, r2R),
							_mm256_mul_ps(C2, raR)
					)
			);
}

static Mm computeE2(const Mm &r1R, const Mm &r2R, const Mm &Ca, const Mm &Cb, const Mm &wv1, const Mm &wv2)
{
	return
			_mm256_mul_ps(
					_mm256_add_ps(
							_mm256_mul_ps(
									_mm256_add_ps(
											_mm256_mul_ps(Ca, wv1),
											_mm256_mul_ps(Cb, wv2)
									),
									r1R
							),
							_mm256_mul_ps(
									_mm256_add_ps(
											_mm256_mul_ps(Cb, wv1),
											_mm256_mul_ps(Ca, wv2)
									),
									r2R
							)
					),
					_mm256_rcp_ps(_mm256_add_ps(wv1, wv2))
			);
}

static void addEnergy(Mm &e, const Mm &e1, const Mm &e2, const Mm &V)
{
	e = _mm256_add_ps(e, _mm256_add_ps(e1, _mm256_sub_ps(e2, V)));
}

static void computeLocalEnergies(Mm *e, const Context &context, const float a, const float b)
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
			e1 = _mm256_set1_ps(computeE1(a, b)),
			C1 = _mm256_set1_ps(2.f),
			C2 = _mm256_set1_ps(-1.f),
			Ca = _mm256_set1_ps(a),
			Cb = _mm256_set1_ps(b);

	size_t n = WArrayD::Count();
	while(n--) {
		const auto
				r1R = _mm256_rcp_ps(*r1++),
				r2R = _mm256_rcp_ps(*r2++),
				raR = _mm256_rcp_ps(*ra++),

				V = computePotential(C1, C2, r1R, r2R, raR),
				e2 = computeE2(r1R, r2R, Ca, Cb, *wv1++, *wv2++);

		addEnergy(*e++, e1, e2, V);
	}
}

#endif //QM_MONTE_CARLO_LOCALENERGY_H
