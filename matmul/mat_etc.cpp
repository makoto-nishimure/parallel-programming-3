#include "mat_etc.hpp"

void matAdd(int size, point pa, point pb,
		const mat mat_a, const mat mat_b, mat& mat_c)
{
	assert(size % 4 == 0);
#ifdef DEBUG
	cout << "ADD" <<endl;
	print_pair(size, pa);
	print_pair(size, pb);
#endif
	int i, j;
	for (i = 0; i < size; ++i) {
		for (j = 0; j < size; ++j) {
			mat_c[i][j] =  mat_a[i][j] + mat_b[i][j];
		}
	}
}

void simd_matAdd(int size, point pa, point pb,
		const mat mat_a, const mat mat_b, mat& mat_c)
{
	assert(size % 4 == 0);
	const int end = size >> 2;
	const int pa4s = pa.second >> 2;
	const int pb4s = pb.second >> 2;
#ifdef DEBUG
	cout << "ADD" <<endl;
	print_pair(size, pa);
	print_pair(size, pb);
	cout << "PA4S:" << pa4s <<endl;
	cout << "PB4S:" << pb4s <<endl;
#endif
	int i, j;
	int tmp_pa, tmp_pb;
	for (i = 0; i < size; ++i) {
		__m256d* a = (__m256d *)mat_a[pa.first + i];
		__m256d* b = (__m256d *)mat_b[pb.first + i];
		__m256d* c = (__m256d *)mat_c[i];

		for (j = 0; j < end; ++j) {
			c[j] = _mm256_add_pd(a[pa4s + j], b[pb4s + j]);
		}
	}
}

void matSub(int size, point pa, point pb,
		const mat mat_a, const mat mat_b, mat& mat_c)
{
	assert(size % 4 == 0);
#ifdef DEBUG
	cout << "SUB" <<endl;
	print_pair(size, pa);
	print_pair(size, pb);
#endif
	int i, j;
	for (i = 0; i < size; ++i) {
		for (j = 0; j < size; ++j) {
			mat_c[i][j] = mat_a[i][j] + mat_b[i][j];
		}
	}
}

void simd_matSub(int size, point pa, point pb,
		const mat mat_a, const mat mat_b, mat& mat_c)
{
	assert(size % 4 == 0);
	const int end = size >> 2;
	const int pa4s = pa.second >> 2;
	const int pb4s = pb.second >> 2;
#ifdef DEBUG
	cout << "SUB" <<endl;
	print_pair(size, pa);
	print_pair(size, pb);
	cout << "PA4S:" << pa4s <<endl;
	cout << "PB4S:" << pb4s <<endl;
#endif
	int i, j;
	int tmp_pa, tmp_pb;
	for (i = 0; i < size; ++i) {
		__m256d* a = (__m256d *)mat_a[pa.first + i];
		__m256d* b = (__m256d *)mat_b[pb.first + i];
		__m256d* c = (__m256d *)mat_c[i];

		for (j = 0; j < end; ++j) {
			c[j] = _mm256_sub_pd(a[pa4s + j], b[pb4s + j]);
		}
	}
}

void matCopy(int size, point p, const mat from, mat& to)
{
	int tmp;
	for (int i = 0; i < size; ++i) {
		tmp = i + p.first;
		for (int j = 0; j < size; ++j) {
			to[i][j] = from[tmp][j+p.second];
		}
	}
}

void set_zero(int size, mat& matrix)
{
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			matrix[i][j] = 0;
}
