#include "mm.hpp"
#include "mat_etc.hpp"
#include "mm_opt.hpp"

#ifdef TIME
extern double mtime;
extern double mmtime;
#endif

typedef struct {
	mat m;
	arr a;
} adds;

adds* allocMat(int, int);
void freeMat(adds *p);

void naive_ma(int size,	const mat x, const arr y, mat& z);
void block_ma(int size,	const mat x, const arr y, mat& z);
void matarrMul(int, const mat, const arr, mat&);
void matarrMul_b(int, const mat, const arr, mat&);
void matarrMul_l(int, const mat, const arr, mat&);
void matarrMul_t(int, const mat, const arr, mat&);
void matarrMul_tl(int, int, int, const mat, const arr, mat&);
void matarrMul_tlb(int, int, int, const mat, const arr, mat&);

void matMul(int, const mat, const mat, mat&);
void block_Mul(int, const mat, const mat, mat&);
void tb_Mul(int, int, int, const mat, const mat, mat&);

void thr_Mul(int, const mat, const mat, mat&);

void str(int, int, const mat, const mat, mat&);
