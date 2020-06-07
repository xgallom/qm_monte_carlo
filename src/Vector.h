//
// Created by xgallom on 4/21/20.
//

#ifndef QM_MONTE_CARLO_SRC_VECTOR_H
#define QM_MONTE_CARLO_SRC_VECTOR_H

#include "format.h"

namespace Dimension {
	enum Enum : size_t {
		X = 0,
		Y,
		Z,

		Count
	};
}

union Vector {
	struct {
		double x, y, z;
	};
	double data[Dimension::Count];
};

namespace format
{
	struct Vector {
		using Object = ::Vector;
		static constexpr auto Name = "Vector";

		inline static void PrintMembers(std::ostream &out, const Object &object, const std::string &offset)
		{
			PrintMember(out, object.x, "x", offset);
			PrintMember(out, object.y, "y", offset);
			PrintMember(out, object.z, "z", offset);
		}
	};
}

inline std::ostream &operator<<(std::ostream &out, const Vector &object)
{
	return format::PrintObject<format::Vector>(out, object);
}

#endif //QM_MONTE_CARLO_SRC_VECTOR_H
