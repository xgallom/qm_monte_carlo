//
// Created by xgallom on 4/17/19.
//

#ifndef QM_MONTE_CARLO_MATRIX_H
#define QM_MONTE_CARLO_MATRIX_H

#include <array>

template<size_t Rows, size_t Cols>
struct Matrix : public std::array<double, Rows * Cols> {
	static constexpr auto
			rows = Rows,
			cols = Cols;

	struct Coords { size_t row, col; };

	Coords operator()(size_t n) { return {n / Cols, n % Cols}; }
	const double &operator()(size_t y, size_t x) const { return this->operator[](y * Cols + x); }
	double &operator()(size_t y, size_t x) { return this->operator[](y * Cols + x); }
};

#endif //QM_MONTE_CARLO_MATRIX_H
