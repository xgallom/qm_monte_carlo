//
// Created by xgallom on 4/29/19.
//

#ifndef QM_MONTE_CARLO_TRANSLATIONS_H
#define QM_MONTE_CARLO_TRANSLATIONS_H

#include "context.h"

static void clear(Mm *to, size_t n)
{
	const auto zero = _mm256_setzero_ps();
	while(n--)
		*to++ = zero;
}

static void generateRandomPositions(Mm *to, size_t n, Random &random)
{
	while(n--)
		*to++ = _mm256_set_ps(
				random.dist(), random.dist(), random.dist(), random.dist(),
				random.dist(), random.dist(), random.dist(), random.dist()
		);
}

static void generateRandomPositions(Context &c, Random &random)
{
	generateRandomPositions(mm(c.x1.x), Batches, random);
	generateRandomPositions(mm(c.x1.y), Batches, random);
	generateRandomPositions(mm(c.x1.z), Batches, random);
	generateRandomPositions(mm(c.x2.x), Batches, random);
	generateRandomPositions(mm(c.x2.y), Batches, random);
	generateRandomPositions(mm(c.x2.z), Batches, random);
}

static void stepPositions(Mm *to, const Mm *from, size_t n, Random &random)
{
	while(n--)
		*to++ = _mm256_add_ps(*from++, _mm256_set_ps(
				random.get(), random.get(), random.get(), random.get(),
				random.get(), random.get(), random.get(), random.get()
		));
}

static void stepPositions(Context &cn, const Context &c, Random &random)
{
	stepPositions(mm(cn.x1.x), mm(c.x1.x), Batches, random);
	stepPositions(mm(cn.x1.y), mm(c.x1.y), Batches, random);
	stepPositions(mm(cn.x1.z), mm(c.x1.z), Batches, random);
	stepPositions(mm(cn.x2.x), mm(c.x2.x), Batches, random);
	stepPositions(mm(cn.x2.y), mm(c.x2.y), Batches, random);
	stepPositions(mm(cn.x2.z), mm(c.x2.z), Batches, random);
}

#endif //QM_MONTE_CARLO_TRANSLATIONS_H
