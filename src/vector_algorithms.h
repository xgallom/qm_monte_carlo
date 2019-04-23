//
// Created by xgallom on 4/17/19.
//

#ifndef QM_MONTE_CARLO_VECTOR_ALGORITHMS_H
#define QM_MONTE_CARLO_VECTOR_ALGORITHMS_H

template<typename T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
inline T internal_floating_type() { return {}; }
template<typename T, std::enable_if_t<!std::is_floating_point<T>::value, bool> = true>
inline auto internal_floating_type() { return internal_floating_type<typename T::value_type>(); }

template<typename F = Vector3D (const Vector3D &)>
inline Vector3D transform(const Vector3D &vector, F f)
{
	Vector3D transformed = {};

	transform(vector.begin(), vector.end(), transformed.begin(), f);

	return transformed;
}

template<typename F = Vector3D (const Vector3D &, const Vector3D &)>
inline Vector3D transform(const Vector3D &left, const Vector3D &right, F f)
{
	Vector3D transformed = {};

	auto l = left.begin(), r = right.begin();
	auto o = transformed.begin();
	const auto end = left.end();

	for(; l != end;)
		*o++ = f(*l++, *r++);

	return transformed;
}

template<typename T, typename F = T (const T &)>
inline Vector<T> transform(const Vector<T> &vector, F f)
{
	Vector<T> transformed = {};
	transformed.resize(vector.size());

	transform(vector.begin(), vector.end(), transformed.begin(), f);

	return transformed;
}

template<typename U, typename T, typename F = U (const T &)>
inline Vector<U> transformTo(const Vector<T> &vector, F f)
{
	Vector<U> transformed = {};
	transformed.resize(vector.size());

	transform<decltype(vector.begin()), decltype(transformed.begin()), F>
	        (vector.begin(), vector.end(), transformed.begin(), f);

	return transformed;
}

inline double accumulate(const Vector3D &vector, double init)
{ return accumulate(vector.begin(), vector.end(), init); }

inline double accumulate(const Vector3D &vector)
{ return accumulate(vector.begin(), vector.end(), 0.); }

template<typename T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
inline T accumulate(const Vector<T> &vector, T init = {})
{ return accumulate(vector.begin(), vector.end(), init); }

template<typename T, std::enable_if_t<!std::is_floating_point<T>::value, bool> = true>
inline auto accumulate(const Vector<T> &vector, decltype(internal_floating_type<T>()) init = {})
{
	for(const auto &v : vector)
		init = accumulate(v, init);

	return init;
}

#endif //QM_MONTE_CARLO_VECTOR_ALGORITHMS_H
