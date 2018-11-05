#include <iostream>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <utility>
#include <sys/time.h>
#include <immintrin.h>
#include <cassert>
#include <thread>

using namespace std;
typedef double* arr;
typedef double** mat;
typedef pair<int,int> point;

#define BASE 700

#ifdef TIME
double mtime = 0;
double mmtime = 0;
#endif 
int b_size = 50;

void str(int, int, const mat, const mat, mat&);
void _str(int, int, const mat, const mat, mat&);

double elapsed_time(struct timeval tp[2])
{
	return tp[1].tv_sec-tp[0].tv_sec+1e-6*(tp[1].tv_usec-tp[0].tv_usec);
}

void set_zero(int size, mat& matrix)
{
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			matrix[i][j] = 0;
}

void print_pair(int size, point pa)
{
	cout << pa.first << "," << pa.second << " to " <<
		pa.first + size - 1 << "," << pa.second + size - 1 << endl;
}

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

void matMul(int size,
		const mat mat_a, const mat mat_b, mat& mat_c)
{
#ifdef TIME
	struct timeval tp[2];
	gettimeofday(tp, 0);
#endif
	assert(size % 8 == 0);

	int i, j, k;
	int ii, jj, kk;
	double tmp;

	for (ii = 0; ii < size; ii+=b_size) {
		int i_limit = (ii + b_size < size) ? ii + b_size : size;
		for (kk = 0; kk < size; kk+=b_size) {
			int k_limit = (kk + b_size < size) ? kk + b_size : size;
			for (jj = 0; jj < size; jj+=b_size) {
				int j_limit = (jj + b_size < size) ? jj + b_size : size;

				for (i = ii; i < i_limit; ++i) {
					for (k = kk; k < k_limit; ++k) {
						tmp = mat_a[i][k];
						for (j = jj; j < j_limit; ++j) {
							mat_c[i][j] += tmp * mat_b[k][j];
						}
					}
				}
			}
		}
	}
#ifdef TIME
	gettimeofday(tp+1, 0);
	mmtime += (double)elapsed_time(tp);
#endif
}

void simd_matMul(int size,
		const mat mat_a, const mat mat_b, mat& mat_c)
{
#ifdef TIME
	struct timeval tp[2];
	gettimeofday(tp, 0);
#endif
	assert(size % 8 == 0);
	const int end = size >> 2;

	int i, j, k;
	double tmp;
	for (i = 0; i < size; ++i) {
		__m256d* c = (__m256d*)mat_c[i];

		for (k = 0; k < size; ++k) {
			__m256d* b = (__m256d*)mat_b[k];
			tmp = mat_a[i][k];	
			__m256d alpha = _mm256_set_pd(tmp, tmp, tmp, tmp);

			for (j = 0; j < end; ++j) {
				c[j] = _mm256_fmadd_pd(alpha, b[j], c[j]);
				++j;
				c[j] = _mm256_fmadd_pd(alpha, b[j], c[j]);
			}
		}
	}
#ifdef TIME
	gettimeofday(tp+1, 0);
	mmtime += (double)elapsed_time(tp);
#endif
}

void block_matMul(int size,
		const mat mat_a, const mat mat_b, mat& mat_c)
{
#ifdef TIME
	struct timeval tp[2];
	gettimeofday(tp, 0);
#endif
	assert(size % 8 == 0);
	const int end = size >> 2;

	int i, j, k;
	int ii, jj, kk;
	double tmp;

	for (ii = 0; ii < size; ii+=b_size) {
		int i_limit = (ii + b_size < size) ? ii + b_size : size;
		for (kk = 0; kk < size; kk+=b_size) {
			int k_limit = (kk + b_size < size) ? kk + b_size : size;

			for (i = ii; i < i_limit; ++i) {
				__m256d* c = (__m256d*)mat_c[i];

				for (k = kk; k < k_limit; ++k) {
					__m256d* b = (__m256d*)mat_b[k];
					tmp = mat_a[i][k];	
					__m256d alpha = _mm256_set_pd(tmp, tmp, tmp, tmp);

					for (j = 0; j < end; ++j) {
						c[j] = _mm256_fmadd_pd(alpha, b[j], c[j]);
						++j;
						c[j] = _mm256_fmadd_pd(alpha, b[j], c[j]);
					}
				}
			}
		}
	}

#ifdef TIME
	gettimeofday(tp+1, 0);
	mmtime += (double)elapsed_time(tp);
#endif
}

void m1to4(int d, int mid, mat mat_a, mat mat_b,
		mat x, mat y, mat m1, mat m2, mat m3, mat m4)
{
	point p0, p1, p2, p3;
	p0 = make_pair(0,0);
	p1 = make_pair(0,mid);
	p2 = make_pair(mid,0);
	p3 = make_pair(mid,mid);
	// M1
	matAdd(mid, p0, p3, mat_a, mat_a, x);
	matAdd(mid, p0, p3, mat_b, mat_b, y);
	set_zero(mid, m1);
	str(d+1,mid, x, y, m1);
#ifdef DEBUG
	for (int i = 0; i < d; ++i)
		cout << "\t";
	cout << "M1 FINISHED" <<endl;
#endif

	// M2
	matAdd(mid, p2, p3, mat_a, mat_a, x);
	matCopy(mid, p0, mat_b, y);
	set_zero(mid, m2);
	str(d+1,mid, x, y, m2);
#ifdef DEBUG
	for (int i = 0; i < d; ++i)
		cout << "\t";
	cout << "M2 FINISHED" <<endl;
#endif

	// M3
	matCopy(mid, p0, mat_a, x);
	matSub(mid, p1, p3, mat_b, mat_b, y);
	set_zero(mid, m3);
	str(d+1, mid, x, y, m3);
#ifdef DEBUG
	for (int i = 0; i < d; ++i)
		cout << "\t";
	cout << "M3 FINISHED" <<endl;
#endif

	// M4
	matCopy(mid, p3, mat_a, x);
	matSub(mid, p2, p0, mat_b, mat_b, y);
	set_zero(mid, m4);
	str(d+1,mid, x, y, m4);
#ifdef DEBUG
	for (int i = 0; i < d; ++i)
		cout << "\t";
	cout << "M4 FINISHED" <<endl;
#endif
}

void m5to7(int d, int mid, mat mat_a, mat mat_b,
		mat x, mat y, mat m5, mat m6, mat m7)
{
	point p0, p1, p2, p3;
	p0 = make_pair(0,0);
	p1 = make_pair(0,mid);
	p2 = make_pair(mid,0);
	p3 = make_pair(mid,mid);

	// M5
	matAdd(mid, p0, p1, mat_a, mat_a, x);
	matCopy(mid, p3, mat_b, y);
	set_zero(mid, m5);
	str(d+1,mid, x, y, m5);
#ifdef DEBUG
	for (int i = 0; i < d; ++i)
		cout << "\t";
	cout << "M5 FINISHED" <<endl;
#endif

	// M6
	matSub(mid, p2, p0, mat_a, mat_a, x);
	matAdd(mid, p0, p1, mat_b, mat_b, y);
	set_zero(mid, m6);
	str(d+1,mid, x, y, m6);
#ifdef DEBUG
	for (int i = 0; i < d; ++i)
		cout << "\t";
	cout << "M6 FINISHED" <<endl;
#endif

	// M7
	matSub(mid, p1, p3, mat_a, mat_a, x);
	matAdd(mid, p2, p3, mat_b, mat_b, y);
	set_zero(mid, m7);
	str(d+1,mid, x, y, m7);
#ifdef DEBUG
	for (int i = 0; i < d; ++i)
		cout << "\t";
	cout << "M7 FINISHED" <<endl;
#endif
}

void str(int d, int size,
		const mat mat_a, const mat mat_b, mat& mat_c)
{
	struct timeval tp[2];
#ifdef DEBUG
	for (int i = 0; i < d; ++i)
		cout << "\t";
	cout << "STR START. SIZE:" << size <<endl;
#endif
	if (size <= BASE) {
		//block_matMul(size, mat_a, mat_b, mat_c);
		matMul(size, mat_a, mat_b, mat_c);
#ifdef DEBUG
		for (int i = 0; i < d; ++i)
			cout << "\t";
		cout << "matMul DONE. SIZE:" << size <<endl;
#endif
		return;
	}
	int i, j, k, mid;
	mid = size >> 1;
	point p0, p1, p2, p3;
	p0 = make_pair(0,0);
	p1 = make_pair(0,mid);
	p2 = make_pair(mid,0);
	p3 = make_pair(mid,mid);

	mat m1, m2, m3, m4, m5, m6, m7, x, y;

#ifdef TIME
	gettimeofday(tp, 0);
#endif
	m1 = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	m2 = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	m3 = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	m4 = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	m5 = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	m6 = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	m7 = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	x = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	y = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	for (i = 0; i < mid; ++i) {
		m1[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		m2[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		m3[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		m4[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		m5[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		m6[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		m7[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		x[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		y[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
	}
#ifdef TIME
	gettimeofday(tp+1, 0);
	mtime += (double)elapsed_time(tp);
#endif

	m1to4(d, mid, mat_a, mat_b, x, y, m1, m2, m3, m4);
	m5to7(d, mid, mat_a, mat_b, x, y, m5, m6, m7);

	// C1 and C2
	// C1
	int tmp_n, tmp_m;
	matAdd(mid, p0, p0, m1, m4, x);
	matSub(mid, p0, p0, x, m5, y);
	matAdd(mid, p0, p0, y, m7, x);
	for (i = 0; i < mid; ++i)
		for (j = 0; j < mid; ++j)
			mat_c[i][j] = x[i][j];
	// C2
	matAdd(mid, p0, p0, m3, m5, x);
	for (i = 0; i < mid; ++i)
		for (j = 0, tmp_m = mid; j < mid; ++j, ++tmp_m)
			mat_c[i][tmp_m] = x[i][j];

#ifdef DEBUG
	for (int i = 0; i < d; ++i)
		cout << "\t";
	cout << "C1 AND C2  FINISHED" <<endl;
#endif
	// C3 and C4
	// C3
	matAdd(mid, p0, p0, m2, m4, x);
	for (i = 0, tmp_n = mid; i < mid; ++i, ++tmp_n) 
		for (j = 0; j < mid; ++j)
			mat_c[tmp_n][j] = x[i][j];

	// C4
	matSub(mid, p0, p0, m1, m2, x);
	matAdd(mid, p0, p0, x, m3, y);
	matAdd(mid, p0, p0, y, m6, x);
	for (i = 0, tmp_n = mid; i < mid; ++i, ++tmp_n)
		for (j = 0, tmp_m = mid; j < mid; ++j, ++tmp_m)
			mat_c[tmp_n][tmp_m] = x[i][j];

#ifdef DEBUG
	for (int i = 0; i < d; ++i)
		cout << "\t";
	cout << "C3 AND C4  FINISHED" <<endl;
	cout << "MID:" << mid <<endl;
#endif
#ifdef TIME
	gettimeofday(tp, 0);
#endif
	for (i = 0; i < mid; ++i) {
		_mm_free(m1[i]); 
		_mm_free(m2[i]);
		_mm_free(m3[i]);
		_mm_free(m4[i]);
		_mm_free(m5[i]);
		_mm_free(m6[i]);
		_mm_free(m7[i]);
		_mm_free(x[i]);
		_mm_free(y[i]);
	}
	_mm_free(m1); 
	_mm_free(m2);
	_mm_free(m3);
	_mm_free(m4);
	_mm_free(m5);
	_mm_free(m6);
	_mm_free(m7);
	_mm_free(x);
	_mm_free(y);
#ifdef TIME
	gettimeofday(tp+1, 0);
	mtime += (double)elapsed_time(tp);
#endif
}

void sub(int d, int mid, mat mat_a, mat mat_b,
		mat m5, mat m6, mat m7)
{
	int i;
	mat x, y;
	x = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	y = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	for (i = 0; i < mid; ++i) {
		x[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		y[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
	}
	m5to7(d, mid, mat_a, mat_b, x, y, m5, m6, m7);
	for (i = 0; i < mid; ++i) {
		_mm_free(x[i]);
		_mm_free(y[i]);
	}
	_mm_free(x);
	_mm_free(y);
}


void _str(int d, int size,
		const mat mat_a, const mat mat_b, mat& mat_c)
{
	struct timeval tp[2];
#ifdef DEBUG
	for (int i = 0; i < d; ++i)
		cout << "\t";
	cout << "STR START. SIZE:" << size <<endl;
#endif
	int i, j, k, mid;
	mid = size >> 1;
	point p0, p1, p2, p3;
	p0 = make_pair(0,0);
	p1 = make_pair(0,mid);
	p2 = make_pair(mid,0);
	p3 = make_pair(mid,mid);
	mat m1, m2, m3, m4, m5, m6, m7, x, y;

#ifdef TIME
	gettimeofday(tp, 0);
#endif
	m1 = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	m2 = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	m3 = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	m4 = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	m5 = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	m6 = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	m7 = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	x = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	y = (double**)_mm_malloc(sizeof(double*)*mid, 64);
	for (i = 0; i < mid; ++i) {
		m1[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		m2[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		m3[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		m4[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		m5[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		m6[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		m7[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		x[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
		y[i] = (double*)_mm_malloc(sizeof(double)*mid, 64);
	}
#ifdef TIME
	gettimeofday(tp+1, 0);
	mtime += (double)elapsed_time(tp);
#endif

	thread t(sub, d, mid, ref(mat_a), ref(mat_b), ref(m5), ref(m6), ref(m7)); 
	m1to4(d, mid, mat_a, mat_b, x, y, m1, m2, m3, m4);
	t.join();

	// C1 and C2
	// C1
	int tmp_n, tmp_m;
	matAdd(mid, p0, p0, m1, m4, x);
	matSub(mid, p0, p0, x, m5, y);
	matAdd(mid, p0, p0, y, m7, x);
	for (i = 0; i < mid; ++i)
		for (j = 0; j < mid; ++j)
			mat_c[i][j] = x[i][j];
	// C2
	matAdd(mid, p0, p0, m3, m5, x);
	for (i = 0; i < mid; ++i)
		for (j = 0, tmp_m = mid; j < mid; ++j, ++tmp_m)
			mat_c[i][tmp_m] = x[i][j];

#ifdef DEBUG
	for (int i = 0; i < d; ++i)
		cout << "\t";
	cout << "C1 AND C2  FINISHED" <<endl;
#endif
	// C3 and C4
	// C3
	matAdd(mid, p0, p0, m2, m4, x);
	for (i = 0, tmp_n = mid; i < mid; ++i, ++tmp_n) 
		for (j = 0; j < mid; ++j)
			mat_c[tmp_n][j] = x[i][j];

	// C4
	matSub(mid, p0, p0, m1, m2, x);
	matAdd(mid, p0, p0, x, m3, y);
	matAdd(mid, p0, p0, y, m6, x);
	for (i = 0, tmp_n = mid; i < mid; ++i, ++tmp_n)
		for (j = 0, tmp_m = mid; j < mid; ++j, ++tmp_m)
			mat_c[tmp_n][tmp_m] = x[i][j];

#ifdef DEBUG
	for (int i = 0; i < d; ++i)
		cout << "\t";
	cout << "C3 AND C4  FINISHED" <<endl;
	cout << "MID:" << mid <<endl;
#endif
#ifdef TIME
	gettimeofday(tp, 0);
#endif
	for (i = 0; i < mid; ++i) {
		_mm_free(m1[i]); 
		_mm_free(m2[i]);
		_mm_free(m3[i]);
		_mm_free(m4[i]);
		_mm_free(m5[i]);
		_mm_free(m6[i]);
		_mm_free(m7[i]);
		_mm_free(x[i]);
		_mm_free(y[i]);
	}
	_mm_free(m1); 
	_mm_free(m2);
	_mm_free(m3);
	_mm_free(m4);
	_mm_free(m5);
	_mm_free(m6);
	_mm_free(m7);
	_mm_free(x);
	_mm_free(y);
	gettimeofday(tp+1, 0);
#ifdef TIME
	mtime += (double)elapsed_time(tp);
#endif
}

void print_mat(int size, mat m)
{
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			cout << m[i][j] << " ";
		}
		cout << endl;
	}
}

int main (int argc, char **argv)
{
	mat m1, m2, m3;
	int n = 1504;
	int N = 1504;
	m1 = (double **)_mm_malloc(sizeof(double*)*N, 64);
	m2 = (double **)_mm_malloc(sizeof(double*)*N, 64);
	m3 = (double **)_mm_malloc(sizeof(double*)*N, 64);
	for (int i = 0; i < N; ++i) {
		m1[i] = (double *)_mm_malloc(sizeof(double)*N, 64);
		m2[i] = (double *)_mm_malloc(sizeof(double)*N, 64);
		m3[i] = (double *)_mm_malloc(sizeof(double)*N, 64);
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			m1[i][j] = ((i + j) * 100 / 7.0);
			m2[i][j] = ((i + j) * 300 / 7.0);
		}
	}
	struct timeval tp[2];
	int roop = 30;
	gettimeofday(tp, 0);
	for (int i = 0; i < roop; ++i) {
		set_zero(n, m3);
		_str(0, n, m1, m2, m3);
	}
	gettimeofday(tp+1, 0);
	//cout << "STR DONE" << endl;
	double one_exe_time = double(elapsed_time(tp)) / roop;
	cout << one_exe_time <<endl;;
	//print_mat(n, m3);
#ifdef TIME
	cout << "malloc and free TIME" <<endl;
	cout << mtime / roop <<endl;
	cout << (mtime / roop) / one_exe_time <<endl;;
	cout << "MM TIME" << endl;
	cout << mmtime / (2 * roop) <<endl;
	cout << (mmtime / (2 * roop)) / one_exe_time <<endl;
#endif

	for (int i = 0; i < N; ++i) {
		_mm_free(m1[i]);
		_mm_free(m2[i]);
		_mm_free(m3[i]);
	}
	_mm_free(m1);
	_mm_free(m2);
	_mm_free(m3);
	return 0;
}
