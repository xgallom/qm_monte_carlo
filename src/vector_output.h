//
// Created by xgallom on 4/17/19.
//

#ifndef QM_MONTE_CARLO_VECTOR_OUTPUT_H
#define QM_MONTE_CARLO_VECTOR_OUTPUT_H

#include <iostream>

static int s_offset = 0;

inline void outputOffset(std::ostream &os)
{
	auto offset = s_offset;
	while(offset--)
		os << "  ";
}

inline std::ostream &operator<<(std::ostream &os, const Vector3D &vector)
{
	os << "{ ";

	for(const auto &v : vector) {
		os << v;

		if(&v + 1 != vector.end())
			os << ", ";
	}

	os << " }";

	return os;
}

template<typename T>
inline std::ostream &operator<<(std::ostream &os, const Vector<T> &vector)
{
	os << "{\n";

	++s_offset;

	for(const auto &v : vector) {
		outputOffset(os);
		os << v;

		if(&v + 1 != &*vector.end())
			os << ',';
		os << '\n';
	}

	--s_offset;

	outputOffset(os);
	os << '}';

	return os;
}

#endif //QM_MONTE_CARLO_VECTOR_OUTPUT_H
