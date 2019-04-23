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

template<size_t Rows, size_t Cols>
inline std::ostream &operator<<(std::ostream &os, const Matrix<Rows, Cols> &matrix)
{
	os << "{\n";

	++s_offset;

	for(size_t y = 0; y < Rows; ++y) {
		outputOffset(os);

		os << "{ ";
		for(size_t x = 0; x < Cols; ++x) {
			os << matrix(y, x);

			if(x + 1 != Cols)
				os << ", ";
		}
		os << " }";

		if(y + 1 != Rows)
			os << ',';
		os << '\n';
	}

	--s_offset;

	outputOffset(os);
	os << '}';

	return os;
}

#endif //QM_MONTE_CARLO_MATRIX_H
