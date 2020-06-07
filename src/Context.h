//
// Created by xgallom on 4/21/20.
//

#ifndef QM_MONTE_CARLO_SRC_CONTEXT_H
#define QM_MONTE_CARLO_SRC_CONTEXT_H

#include "Configuration.h"
#include "ScalarConfiguration.h"
#include "Parameters.h"

struct DistanceContext {
	ElectronProtonConfiguration distances;
	double electronDistance;
};

struct WaveContext {
	ElectronProtonConfiguration electronProtonWaves;
	ExtendedScalarConfiguration electronWaves;
	double wave, waveSquared;
};

struct Context {
	Configuration electronConfiguration;
	ProtonConfiguration protonConfiguration;

	DistanceContext distanceContext;
	WaveContext waveContext;
};

DistanceContext generateContext(const Configuration &configuration, const ProtonConfiguration &protonConfiguration);
WaveContext generateContext(
		const DistanceContext &distanceContext,
		const Parameters &parameters
);
Context generateContext(
		std::mt19937 &generator,
		std::uniform_real_distribution<double> &distribution,
		const Parameters &parameters,
		double S
);
Context generateContext(
		const Context &oldContext,
		std::mt19937 &generator,
		std::uniform_real_distribution<double> &distribution,
		const Parameters &parameters,
		double S
);

namespace format {
	struct DistanceContext {
		using Object = ::DistanceContext;
		static constexpr auto Name = "DistanceContext";

		inline static void PrintMembers(std::ostream &out, const Object &object, const std::string &offset)
		{
			PrintMember(distances);
			PrintMember(electronDistance);
		}
	};

	struct WaveContext {
		using Object = ::WaveContext;
		static constexpr auto Name = "WaveContext";

		inline static void PrintMembers(std::ostream &out, const Object &object, const std::string &offset)
		{
			PrintMember(electronWaves);
			PrintMember(electronProtonWaves);
			PrintMember(waveSquared);
			PrintMember(wave);
		}
	};

	struct Context {
		using Object = ::Context;
		static constexpr auto Name = "Context";

		inline static void PrintMembers(std::ostream &out, const Object &object, const std::string &offset)
		{
			PrintMember(electronConfiguration);
			PrintMember(protonConfiguration);
			PrintMember(distanceContext);
			PrintMember(waveContext);
		}
	};
}

inline std::ostream &operator<<(std::ostream &out, const DistanceContext &object)
{
	return format::PrintObject<format::DistanceContext>(out, object);
}

inline std::ostream &operator<<(std::ostream &out, const WaveContext &object)
{
	return format::PrintObject<format::WaveContext>(out, object);
}

inline std::ostream &operator<<(std::ostream &out, const Context &object)
{
	return format::PrintObject<format::Context>(out, object);
}

#endif //QM_MONTE_CARLO_SRC_CONTEXT_H
