#include "mm.hpp"

void matAdd(int, point, point,
		const mat, const mat, mat&);

void simd_matAdd(int, point, point,
		const mat, const mat, mat&);

void matSub(int, point, point,
		const mat, const mat, mat&);

void simd_matSub(int, point, point,
		const mat, const mat, mat&);

void matCopy(int, point, const mat, mat&);

void set_zero(int, mat&);
