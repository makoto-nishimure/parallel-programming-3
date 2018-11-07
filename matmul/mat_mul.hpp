#include "mm.hpp"
#include "mat_etc.hpp"
#include "mm_opt.hpp"

#ifdef TIME
extern double mtime;
extern double mmtime;
#endif

void matMul(int, const mat, const mat, mat&);
void simd_matMul(int, const mat, const mat, mat&);
void block_matMul(int, const mat, const mat, mat&);

void m1to4(int, int, mat, mat,
		mat, mat, mat, mat, mat, mat);

void m5to7(int, int, mat, mat,
		mat, mat, mat, mat, mat);

void str(int, int, const mat, const mat, mat&);
void _str(int, int, const mat, const mat, mat&);
void sub(int, int, mat, mat, mat, mat, mat);
