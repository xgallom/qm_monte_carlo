//
// Created by xgallom on 4/21/20.
//

#include "transformations.h"

double sum(const double *data, size_t count)
{
	double result = 0.;

	while(count--)
		result += *data++;

	return result;
}

double sum(const ExtendedScalarConfiguration &object)
{
	return sum(object.data, ExtendedElectron::Count);
}

double sumOfReciprocals(const double *data, size_t count)
{
	double result = 0.;

	while(count--)
		result += 1. / *data++;

	return result;
}

double sumOfReciprocals(const ElectronProtonConfiguration &object)
{
	return sumOfReciprocals(reinterpret_cast<const double *>(object.data), Proton::Count * Electron::Count);
}

double square(const double &value)
{
	return value * value;
}

double dot(const Vector &left, const Vector &right)
{
	return left.x * right.x + left.y * right.y + left.z * right.z;
}

double projection(Vector direction, const Vector &vector, const double &directionLength)
{
	normalize(direction, directionLength);
	return dot(direction, vector);
}

void normalize(Vector &vector, const double &length)
{
	vector.x /= length;
	vector.y /= length;
	vector.z /= length;
}
