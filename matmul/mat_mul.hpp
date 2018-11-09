#include "mm.hpp"
#include "mat_etc.hpp"
#include "mm_opt.hpp"

#ifdef TIME
extern double mtime;
extern double mmtime;
#endif

void matMul(int, const mat, const mat, mat&);
void block_Mul(int, const mat, const mat, mat&);
void tb_Mul(int, int, int, const mat, const mat, mat&);

void thr_Mul(int, const mat, const mat, mat&);

void str(int, int, const mat, const mat, mat&);
