//
// Created by xgallom on 4/16/19.
//

#ifndef QM_MONTE_CARLO_VECTOR_H
#define QM_MONTE_CARLO_VECTOR_H

#include "transforms.h"

#include <array>
#include <vector>
#include <algorithm>
#include <numeric>

using std::transform;
using std::accumulate;
using std::swap;

template<typename T>
using Vector = std::vector<T>;

static constexpr size_t SpaceDimensions = 3;
using Vector3D = std::array<double, SpaceDimensions>;

using VectorD = Vector<double>;
using VectorV3D = Vector<Vector3D>;
using VectorIdx = Vector<size_t>;

#include "vector_algorithms.h"
#include "vector_output.h"
#include "vector_operations.h"

inline Vector3D sqr(const Vector3D &vector)
{ return transform(vector, tr::sqr); }

inline double absSqr(const Vector3D &vector)
{ return accumulate(sqr(vector)); }

inline double abs(const Vector3D &vector)
{ return sqrt(absSqr(vector)); }

#endif //QM_MONTE_CARLO_VECTOR_H
