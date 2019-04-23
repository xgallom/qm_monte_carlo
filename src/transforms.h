//
// Created by xgallom on 4/16/19.
//

#ifndef QM_MONTE_CARLO_TRANSFORMS_H
#define QM_MONTE_CARLO_TRANSFORMS_H

#include <cmath>
#include <functional>

namespace tr
{
	using tr = std::function<double(const double &)>;
	using op = std::function<double(const double &, const double &)>;

	constexpr double sqr(const double &v)
	{ return v * v; }

	inline double sqrt(const double &x)
	{ return std::sqrt(x); }

	inline tr pow(const double &p)
	{ return [p](const double &v) -> double { return std::pow(v, p); }; }

	inline double add(const double &l, const double &r) { return l + r; }
	inline double sub(const double &l, const double &r) { return l - r; }

	inline op mul() { return [](const double &l, const double &r) { return l * r;}; }
	inline tr mul(const double &r) { return [r](const double &l) { return l * r; }; }
	inline op div() { return [](const double &l, const double &r) { return l / r; }; }
	inline tr div(const double &r) { return [r](const double &l) { return l / r; }; }
}

using tr::sqr;

#endif //QM_MONTE_CARLO_TRANSFORMS_H
