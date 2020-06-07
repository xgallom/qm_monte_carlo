//
// Created by xgallom on 4/21/20.
//

#ifndef QM_MONTE_CARLO_SRC_TRANSFORMATIONS_H
#define QM_MONTE_CARLO_SRC_TRANSFORMATIONS_H

#include "ScalarConfiguration.h"

double sum(const double *data, size_t count);
double sum(const ExtendedScalarConfiguration &object);

double sumOfReciprocals(const double *data, size_t count);
double sumOfReciprocals(const ElectronProtonConfiguration &object);

double square(const double &value);

double dot(const Vector &left, const Vector &right);
double projection(Vector direction, const Vector &vector, const double &directionLength);
void normalize(Vector &vector, const double &length);

#endif //QM_MONTE_CARLO_SRC_TRANSFORMATIONS_H
