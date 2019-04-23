//
// Created by xgallom on 4/16/19.
//

#include "random.h"
#include "config.h"
#include <random>
#include <memory>

namespace Random
{
	namespace
	{
		struct Context {
			Context(std::mt19937 &&a_engine, std::uniform_real_distribution<double> &&a_dist) :
					engine(a_engine),
					dist(a_dist)
			{}

			std::mt19937 engine;
			std::uniform_real_distribution<double> dist;
		};

		std::unique_ptr<Context> s_context;
	}

	void init()
	{
		std::random_device device;

		s_context = std::make_unique<Context>(std::mt19937(device()), std::uniform_real_distribution<double>(-1., 1.));
	}

	double get()
	{ return s_context->dist(s_context->engine); }

	double dist()
	{ return Config::Dist * get(); }

	Vector3D vector()
	{ return {dist(), dist(), dist()}; }

	Vector3D vectorOffset()
	{ return {get(), get(), get()}; }
}
