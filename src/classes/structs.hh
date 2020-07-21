#ifndef EMC_structs_hh_included
#define EMC_structs_hh_included

#include <fstream>
#include <math.h>
#include <map>

namespace EMC {

enum matrixType {
	mtNumber, mtProbability, mtDensity
};

typedef long double real;

struct ValueError {
	inline ValueError() :
			value(0), err_sq(0) {
	}
	inline ValueError(real value) :
			value(value), err_sq(0) {
	}
	inline ValueError(real value, real err_sq) :
			value(value), err_sq(err_sq) {
	}

	real value;
	real err_sq;

	inline ValueError operator+(const ValueError& other) const {
		return ValueError(value + other.value, other.err_sq + err_sq);
	}
	inline ValueError operator-(const ValueError& other) const {
		return ValueError(value - other.value, other.err_sq + err_sq);
	}
	inline ValueError operator*(const ValueError& other) const {
		return ValueError(value * other.value,
				value * value * other.err_sq
						+ other.value * other.value * err_sq);
	}
	inline ValueError operator/(const ValueError& other) const {
		return ValueError(value / other.value,
				(other.err_sq
						+ (value * value / (other.value * other.value)) * err_sq)
						/ (other.value * other.value));
	}
	inline bool operator!=(const ValueError& other) const {
		return other.value != value;
	}

	inline ValueError operator+(real other) const {
		return ValueError(value + other, err_sq);
	}
	inline ValueError operator-(real other) const {
		return ValueError(value - other, err_sq);
	}
	inline ValueError operator*(real other) const {
		return ValueError(value * other, err_sq * other * other);
	}
	inline ValueError operator/(real other) const {
		return ValueError(value / other, err_sq / (other * other));
	}

	inline void operator+=(const ValueError& other) {
		err_sq += other.err_sq;
		value += other.value;
	}
	inline void operator-=(const ValueError& other) {
		err_sq += other.err_sq;
		value -= other.value;
	}
	inline void operator*=(const ValueError& other) {
		err_sq = value * value * other.err_sq
				+ other.value * other.value * err_sq;
		value *= other.value;
	}
	inline void operator/=(const ValueError& other) {
		err_sq = (other.err_sq
				+ (value * value / (other.value * other.value)) * err_sq)
				/ (other.value * other.value);
		value /= other.value;
	}

	inline void operator+=(real other) {
		value += other;
	}
	inline void operator-=(real other) {
		value += other;
	}
	inline void operator*=(real other) {
		err_sq *= other * other;
		value *= other;
	}
	inline void operator/=(real other) {
		err_sq /= other * other;
		value /= value;
	}
};

inline ValueError operator+(real a, const ValueError& b) {
	return ValueError(a + b.value, b.err_sq);
}
inline ValueError operator-(real a, const ValueError& b) {
	return ValueError(a - b.value, b.err_sq);
}
inline ValueError operator*(real a, const ValueError& b) {
	return ValueError(a * b.value, a * a * b.err_sq);
}
inline ValueError operator/(real a, const ValueError& b) {
	return ValueError(a / b.value, (b.err_sq / (b.value * b.value)) / (a * a));
}

} /* namespace EMC */

#endif /* EMC_structs_hh_included */
