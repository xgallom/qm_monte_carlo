//
// Created by xgallom on 4/21/20.
//

#ifndef QM_MONTE_CARLO_SRC_SCALARCONFIGURATION_H
#define QM_MONTE_CARLO_SRC_SCALARCONFIGURATION_H

#include "Configuration.h"

union ScalarConfiguration {
	struct {
		double e1, e2;
	};
	double data[Electron::Count];
};

union ExtendedScalarConfiguration {
	struct {
		double e1, e2, f;
	};
	double data[ExtendedElectron::Count];
};

union ElectronProtonConfiguration {
	struct {
		ScalarConfiguration p1, p2;
	};
	double data[Proton::Count][Electron::Count];
};

double distance(const Vector &vector);
double distance(const Vector &vector, const Vector &origin);
ScalarConfiguration distance(const Configuration &configuration, const Vector &origin);
ElectronProtonConfiguration distances(const Configuration &configuration, const ProtonConfiguration &origin);

double waveFunction(const double &distance, const double &a);
ScalarConfiguration waveFunction(const ScalarConfiguration &distanceConfiguration, const double &a);
ElectronProtonConfiguration waveFunction(const ElectronProtonConfiguration &distanceConfiguration, const double &a);
ExtendedScalarConfiguration waveFunction(
		const ElectronProtonConfiguration &electronProtonConfiguration,
		const double &electronDistance, const double &alpha, const double &beta
);

namespace format {
	struct ScalarConfiguration {
		using Object = ::ScalarConfiguration;
		static constexpr auto Name = "ScalarConfiguration";

		inline static void PrintMembers(std::ostream &out, const Object &object, const std::string &offset)
		{
			PrintMember(out, object.e1, "e1", offset);
			PrintMember(out, object.e2, "e2", offset);
		}
	};

	struct ExtendedScalarConfiguration {
		using Object = ::ExtendedScalarConfiguration;
		static constexpr auto Name = "ExtendedScalarConfiguration";

		inline static void PrintMembers(std::ostream &out, const Object &object, const std::string &offset)
		{
			PrintMember(out, object.e1, "e1", offset);
			PrintMember(out, object.e2, "e2", offset);
			PrintMember(out, object.f, "f", offset);
		}
	};

	struct ElectronProtonConfiguration {
		using Object = ::ElectronProtonConfiguration;
		static constexpr auto Name = "ElectronProtonConfiguration";

		inline static void PrintMembers(std::ostream &out, const Object &object, const std::string &offset)
		{
			PrintMember(out, object.p1, "p1", offset);
			PrintMember(out, object.p2, "p2", offset);
		}
	};
}

inline std::ostream &operator<<(std::ostream &out, const ScalarConfiguration &object)
{
	return format::PrintObject<format::ScalarConfiguration>(out, object);
}

inline std::ostream &operator<<(std::ostream &out, const ExtendedScalarConfiguration &object)
{
	return format::PrintObject<format::ExtendedScalarConfiguration>(out, object);
}

inline std::ostream &operator<<(std::ostream &out, const ElectronProtonConfiguration &object)
{
	return format::PrintObject<format::ElectronProtonConfiguration>(out, object);
}

#endif //QM_MONTE_CARLO_SRC_SCALARCONFIGURATION_H
