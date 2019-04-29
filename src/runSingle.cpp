//
// Created by xgallom on 4/29/19.
//

#include "runSingle.h"
#include "energyForParameters.h"

void runSingle(const float *a, const float *b, float *e, size_t t)
{
	Random random;

	while(t--)
		*e++ = energyForParameters(*a++, *b++, random);
}

