//
// Created by xgallom on 4/21/20.
//

#ifndef QM_MONTE_CARLO_SRC_CONFIGURATION_H
#define QM_MONTE_CARLO_SRC_CONFIGURATION_H

#include "Vector.h"

#include <random>
#include <cstring>

namespace Electron {
	enum Enum : size_t {
		E1 = 0,
		E2,

		Count
	};
}

namespace ExtendedElectron {
	enum Enum : size_t {
		E1 = 0,
		E2,
		F,

		Count
	};
}

namespace Proton {
	enum Enum : size_t {
		P1 = 0,
		P2,

		Count
	};
}

union Configuration {
	struct {
		Vector e1, e2;
	};
	double data[Electron::Count][Dimension::Count];
};

union ProtonConfiguration {
	struct {
		Vector p1, p2;
	};
	double data[Proton::Count][Dimension::Count];
};

Configuration generateConfiguration(
		std::mt19937 &generator,
		std::uniform_real_distribution<double> &distribution,
		double protonOffset
);
Configuration generateConfiguration(
		const Configuration &oldConfiguration,
		std::mt19937 &generator,
		std::uniform_real_distribution<double> &distribution,
		double protonDistance
);
ProtonConfiguration generateConfiguration(double protonOffset);

namespace format
{
	struct Configuration {
		using Object = ::Configuration;
		static constexpr auto Name = "Configuration";

		inline static void PrintMembers(std::ostream &out, const Object &object, const std::string &offset)
		{
			PrintMember(out, object.e1, "e1", offset);
			PrintMember(out, object.e2, "e2", offset);
		}
	};

	struct ProtonConfiguration {
		using Object = ::ProtonConfiguration;
		static constexpr auto Name = "ProtonConfiguration";

		inline static void PrintMembers(std::ostream &out, const Object &object, const std::string &offset)
		{
			PrintMember(out, object.p1, "p1", offset);
			PrintMember(out, object.p2, "p2", offset);
		}
	};
}

inline std::ostream &operator<<(std::ostream &out, const Configuration &object)
{
	return format::PrintObject<format::Configuration>(out, object);
}

inline std::ostream &operator<<(std::ostream &out, const ProtonConfiguration &object)
{
	return format::PrintObject<format::ProtonConfiguration>(out, object);
}

#endif //QM_MONTE_CARLO_SRC_CONFIGURATION_H
