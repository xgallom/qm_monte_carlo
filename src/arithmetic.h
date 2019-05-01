//
// Created by xgallom on 4/29/19.
//

#ifndef QM_MONTE_CARLO_ARITHMETIC_H
#define QM_MONTE_CARLO_ARITHMETIC_H

static float sqr(const float x)
{ return x * x; }

#include <iostream>

static Mm abs(const Mm &x, const Mm &y, const Mm &z)
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

#endif //QM_MONTE_CARLO_ARITHMETIC_H
