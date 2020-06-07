//
// Created by xgallom on 4/21/20.
//

#include "Parameters.h"
#include "constants.h"

#include <cstddef>
#include <iostream>
#include <cmath>
#include <limits>

static double generateParameterA(double protonDistance)
{
	double a = BohrRadius, oldA = std::numeric_limits<double>::max();

	while(std::fabs(a - oldA) >= std::numeric_limits<double>::epsilon()) {
		oldA = a;

		const double
				exponential = exp(-protonDistance / a),
				positiveExponential = exp(protonDistance / a),
				parentheses = 1. + exponential,
				positiveParentheses = 1. + positiveExponential,
				derivativeParentheses = a * positiveParentheses;

		a = a
			+ (BohrRadius / parentheses - a)
			  /
			  (BohrRadius * protonDistance * positiveExponential / (derivativeParentheses * derivativeParentheses)
			   + 1.
			  );
	}

	return a;
}

Parameters generateParameters(double protonDistance)
{
	return {{
					.a = generateParameterA(protonDistance),
					.alpha = 2. * BohrRadius,
					.beta = 0.25 / BohrRadius
			}};
}
