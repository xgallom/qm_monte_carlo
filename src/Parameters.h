//
// Created by xgallom on 4/21/20.
//

#ifndef QM_MONTE_CARLO_SRC_PARAMETERS_H
#define QM_MONTE_CARLO_SRC_PARAMETERS_H

#include "format.h"

namespace Parameter {
	enum Enum {
		A = 0,
		Alpha,
		Beta,

		Count
	};
}

union Parameters {
	struct {
		double a, alpha, beta;
	};
	double data[Parameter::Count];
};

Parameters generateParameters(double protonDistance);

namespace format {
	struct Parameters {
		using Object = ::Parameters;
		static constexpr auto Name = "Parameters";

		inline static void PrintMembers(std::ostream &out, const Object &object, const std::string &offset)
		{
			PrintMember(out, object.a, "a", offset);
			PrintMember(out, object.alpha, "alpha", offset);
			PrintMember(out, object.beta, "beta", offset);
		}
	};
}

inline std::ostream &operator<<(std::ostream &out, const Parameters &object)
{
	return format::PrintObject<format::Parameters>(out, object);
}

#endif //QM_MONTE_CARLO_SRC_PARAMETERS_H
