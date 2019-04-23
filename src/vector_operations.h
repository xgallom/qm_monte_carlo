//
// Created by xgallom on 4/17/19.
//

#ifndef QM_MONTE_CARLO_VECTOR_OPERATIONS_H
#define QM_MONTE_CARLO_VECTOR_OPERATIONS_H

#include "vector.h"

inline Vector3D operator+(const Vector3D &l, const Vector3D &r) { return transform(l, r, tr::add); }
inline Vector3D operator-(const Vector3D &l, const Vector3D &r) { return transform(l, r, tr::sub); }
inline Vector3D operator*(const Vector3D &l, const Vector3D &r) { return transform(l, r, tr::mul()); }
inline Vector3D operator*(const Vector3D &l, const double &r) { return transform(l, tr::mul(r)); }
inline Vector3D operator/(const Vector3D &l, const Vector3D &r) { return transform(l, r, tr::div()); }
inline Vector3D operator/(const Vector3D &l, const double &r) { return transform(l, tr::div(r)); }

#endif //QM_MONTE_CARLO_VECTOR_OPERATIONS_H
