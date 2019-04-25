//
// Created by xgallom on 4/24/19.
//

#include "energyForParameters.h"
#include "random.h"
#include "config.h"

static size_t A = 0, B = 0;

static inline void abs(VectorD &rs,
					   const Vector3D *x)
{
	for(auto &r : rs)
		r = abs(*x++);
}

static inline void wave(VectorD &wvs,
						const double *r1, const double *r2,
						const double a, const double b)
{
	for(auto &wv : wvs)
		wv = exp(-a * *r1++ - b * *r2++);
}

static inline void update(VectorD &r1s, VectorD &r2s, VectorD &wv1s, VectorD &wv2s,
						  const Vector3D *x1, const Vector3D *x2,
						  const double a, const double b)
{
	abs(r1s, x1);
	abs(r2s, x2);
	wave(wv1s, r1s.data(), r2s.data(), a, b);
	wave(wv2s, r1s.data(), r2s.data(), b, a);
}

static inline void step(VectorV3D &xns,
						const Vector3D *x)
{
	for(auto &xn : xns)
		xn = *x++ + Random::vectorOffset() * Config::dR;
}

static inline void step(VectorV3D &xn1s, VectorV3D &xn2s,
						const Vector3D *x1s, const Vector3D *x2s)
{
	step(xn1s, x1s);
	step(xn2s, x2s);
}

static inline void prob(VectorD &ps,
						const double *wv1, const double *wv2)
{
	for(auto &p : ps)
		p = sqr(*wv1++ + *wv2++);
}

static inline double energy(const Vector3D *x1, const Vector3D *x2, const double *r1,
							const double *r2, const double *wv1, const double *wv2,
							const double a, const double b)
{
	++A;
	const auto ra = abs(*x2 - *x1);
	const auto V = 2. / *r1 + 2. / *r2 - 1. / ra;
	const auto wv = *wv1 + *wv2;

	return -0.5 * (sqr(a) + sqr(b))
		   + (
					 (a * *wv1 + b * *wv2) / *r1 +
					 (b * *wv1 + a * *wv2) / *r2
			 ) / wv
		   - V;
}

static inline void energy(VectorD &es,
						  const Vector3D *x1, const Vector3D *x2, const double *r1,
						  const double *r2, const double *wv1, const double *wv2,
						  const double a, const double b)
{
	for(auto &e : es) {
		e = energy(x1++, x2++, r1++, r2++, wv1++, wv2++, a, b);
	}
}

static inline void cmp(VectorIdx &updates,
					   const double *p, const double *pn)
{
	updates.clear();
	for(size_t n = 0; n < Config::WalkerCount; ++n) {
		if(*pn++ / *p++ >= Random::norm())
			updates.push_back(n);
	}
}

static inline void transfer(Vector3D *x1s, Vector3D *x2s, double *r1s,
							double *r2s, double *wv1s, double *wv2s, double *ps,
							const Vector3D *xn1s, const Vector3D *xn2s, const double *rn1s,
							const double *rn2s, const double *wvn1s, const double *wvn2s, const double *pns,
							const VectorIdx &updates)
{
	for(const auto idx : updates) {
		x1s[idx] = xn1s[idx];
		x2s[idx] = xn2s[idx];
		r1s[idx] = rn1s[idx];
		r2s[idx] = rn2s[idx];
		wv1s[idx] = wvn1s[idx];
		wv2s[idx] = wvn2s[idx];
		ps[idx] = pns[idx];
	}
}

static inline void add(VectorD &energies,
					   const double *e)
{
	for(auto &energy : energies) {
		++B;
		energy += *e++;
	}
}

static inline void run(VectorV3D &x1s, VectorV3D &x2s, VectorD &r1s, VectorD &r2s,
					   VectorD &wv1s, VectorD &wv2s, VectorD &ps,
					   VectorV3D &xn1s, VectorV3D &xn2s, VectorD &rn1s, VectorD &rn2s,
					   VectorD &wvn1s, VectorD &wvn2s, VectorD &pns,
					   VectorIdx &updates, const double a, const double b)
{
	step(xn1s, xn2s, x1s.data(), x2s.data());
	update(rn1s, rn2s, wvn1s, wvn2s, xn1s.data(), xn2s.data(), a, b);
	prob(pns, wvn1s.data(), wvn2s.data());

	cmp(updates, ps.data(), pns.data());
	transfer(x1s.data(), x2s.data(), r1s.data(), r2s.data(), wv1s.data(), wv2s.data(), ps.data(),
			 xn1s.data(), xn2s.data(), rn1s.data(), rn2s.data(), wvn1s.data(), wvn2s.data(), pns.data(),
			 updates);
}

double energyForParameters(const double &a, const double &b)
{
	A = B = 0;

	VectorV3D
			x1s(Config::WalkerCount, Random::vector()),
			x2s(Config::WalkerCount, Random::vector());
	VectorD
			r1s(Config::WalkerCount),
			r2s(Config::WalkerCount),
			wv1s(Config::WalkerCount),
			wv2s(Config::WalkerCount),
			ps(Config::WalkerCount),
			es(Config::WalkerCount);

	VectorIdx updates;
	updates.reserve(Config::WalkerCount);

	update(r1s, r2s, wv1s, wv2s, x1s.data(), x2s.data(), a, b);
	prob(ps, wv1s.data(), wv2s.data());

	VectorV3D
			xn1s(Config::WalkerCount),
			xn2s(Config::WalkerCount);
	VectorD
			rn1s(Config::WalkerCount),
			rn2s(Config::WalkerCount),
			wvn1s(Config::WalkerCount),
			wvn2s(Config::WalkerCount),
			pns(Config::WalkerCount);

	VectorD energies(Config::WalkerCount);

	for(size_t n = 0; n < Config::Therm; ++n)
		run(x1s, x2s, r1s, r2s, wv1s, wv2s, ps, xn1s, xn2s, rn1s, rn2s, wvn1s, wvn2s, pns, updates, a, b);

	for(size_t n = 0; n < Config::PointsCount; ++n) {
		energy(es, x1s.data(), x2s.data(), r1s.data(), r2s.data(), wv1s.data(), wv2s.data(), a, b);
		add(energies, es.data());

		for(size_t s = 0; s < Config::SkipSteps; ++s)
			run(x1s, x2s, r1s, r2s, wv1s, wv2s, ps, xn1s, xn2s, rn1s, rn2s, wvn1s, wvn2s, pns, updates, a, b);
	}

	std::cout << A << " " << B << "\n";
	return accumulate(energies) / Config::EnergiesCount;
}
