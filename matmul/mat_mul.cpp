#include "mat_mul.hpp"

#define BASE 700

#ifdef TIME
double mtime = 0;
double mmtime = 0;
#endif 
int block_size = 200;

typedef struct {
	mat m;
	arr a;
} adds;

adds* allocMat(int row, int column)
{
	adds *p;
	p = (adds*)malloc(sizeof(adds));
	p->m = (double**)_mm_malloc(sizeof(double*)*column, 64);
	p->a = (double*)_mm_malloc(sizeof(double)*row*column, 64);
	for (int i = 0; i < column; ++i)
		p->m[i] = p->a + i * row;
	return p;
}

void freeMat(adds *p)
{
	_mm_free(ad->m);
	_mm_free(ad->a);
	free(p);
}

void matMul(int size,
		const mat mat_a, const mat mat_b, mat& mat_c)
{
	assert(size % 4 == 0);
	const int end = size >> 2;

	int i, j, k;
	double tmp;
	for (i = 0; i < size; ++i) {
		__m256d* c = (__m256d*)mat_c[i];
		for (k = 0; k < size; ++k) {
			__m256d* b = (__m256d*)mat_b[k];
			tmp = mat_a[i][k];	
			const __m256d alpha = _mm256_set_pd(tmp, tmp, tmp, tmp);
			for (j = 0; j < end; ++j) {
				c[j] = _mm256_fmadd_pd(alpha, b[j], c[j]);
			}
		}
	}
}

void block_Mul(int size,
		const mat mat_a, const mat mat_b, mat& mat_c)
{
	assert(size % 4 == 0);
	const int end = size >> 2;

	int i, j, k;
	int ii, kk;
	double tmp;
	for (ii = 0; ii < size; ii+=block_size) {
		int i_limit = (ii + block_size < size) ? ii + block_size : size;
		for (kk = 0; kk < size; kk+=block_size) {
			int k_limit = (kk + block_size < size) ? kk + block_size : size;

			for (i = ii; i < i_limit; ++i) {
				__m256d* c = (__m256d*)mat_c[i];
				for (k = kk; k < k_limit; ++k) {
					__m256d* b = (__m256d*)mat_b[k];
					tmp = mat_a[i][k];	
					const __m256d alpha = _mm256_set_pd(tmp, tmp, tmp, tmp);
					for (j = 0; j < end; ++j) {
						c[j] = _mm256_fmadd_pd(alpha, b[j], c[j]);
					}
				}
			}
		}
	}
}

void tb_Mul(int size, int thr_s, int thr_e,
		const mat mat_a, const mat mat_b, mat& mat_c)
{
	assert(size % 4 == 0);
	const int end = size >> 2;

	int i, j, k;
	int ii, kk;
	double tmp;
	for (ii = thr_s; ii < thr_e; ii+=block_size) {
		int i_limit = (ii + block_size < thr_e) ? ii + block_size : thr_e;
		for (kk = 0; kk < size; kk+=block_size) {
			int k_limit = (kk + block_size < size) ? kk + block_size : size;

			for (i = ii; i < i_limit; ++i) {
				__m256d* c = (__m256d*)mat_c[i];
				for (k = kk; k < k_limit; ++k) {
					__m256d* b = (__m256d*)mat_b[k];
					tmp = mat_a[i][k];	
					const __m256d alpha = _mm256_set_pd(tmp, tmp, tmp, tmp);
					for (j = 0; j < end; ++j) {
						c[j] = _mm256_fmadd_pd(alpha, b[j], c[j]);
						++j;
						c[j] = _mm256_fmadd_pd(alpha, b[j], c[j]);
					}
				}
			}
		}
	}
}

void thr_Mul(int size,
		const mat a, const mat b, mat& c)
{
	assert(size % 8 == 0);
	thread t1(tb_Mul, size, 0, size >> 1, ref(a), ref(b), ref(c));
	tb_Mul(size, size >> 1, size, ref(a), ref(b), ref(c));
	t1.join();
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
		//		matMul(size, mat_a, mat_b, mat_c);
		//		block_Mul(size, mat_a, mat_b, mat_c);
		thr_Mul(size, mat_a, mat_b, mat_c);
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
	adds *ad1, *ad2, *ad3, *ad4, *ad5, *ad6, *ad7, *adx, *ady;

#ifdef TIME
	gettimeofday(tp, 0);
#endif
	ad1 = alloc_mat(N, N);
	ad2 = alloc_mat(N, N);
	ad3 = alloc_mat(N, N);
	ad4 = alloc_mat(N, N);
	ad5 = alloc_mat(N, N);
	ad6 = alloc_mat(N, N);
	ad7 = alloc_mat(N, N);
	adx = alloc_mat(N, N);
	ady = alloc_mat(N, N);
	m1 = ad1->m;
	m2 = ad2->m;
	m3 = ad3->m;
	m4 = ad4->m;
	m5 = ad5->m;
	m6 = ad6->m;
	m7 = ad7->m;
	x = adx->m;
	y = ady->m;

#ifdef TIME
	gettimeofday(tp+1, 0);
	mtime += (double)elapsed_time(tp);
#endif
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

	//	int tmp_n, tmp_m;
	//	for (i = 0; i < mid; ++i) {
	//		// C1
	//		for (j = 0; j < mid; ++j)
	//			mat_c[i][j] = m1[i][j] + m4[i][j] - m5[i][j] + m7[i][j];
	//		// C2
	//		for (j = 0, tmp_m = mid; j < mid; ++j, ++tmp_m) {
	//			mat_c[i][tmp_m] = m3[i][j] + m5[i][j];
	//		}
	//	}
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

	//	for (i = 0, tmp_n = mid; i < mid; ++i, ++tmp_n) {
	//		// C3
	//		for (j = 0; j < mid; ++j) {
	//			mat_c[tmp_n][j] = m2[i][j] + m4[i][j];
	//		}
	//		// C4
	//		for (j = 0, tmp_m = mid; j < mid; ++j, ++tmp_m) {
	//			mat_c[tmp_n][tmp_m] = m1[i][j] - m2[i][j] + m3[i][j] + m6[i][j];
	//		}
	//	}
#ifdef DEBUG
	for (int i = 0; i < d; ++i)
		cout << "\t";
	cout << "C3 AND C4  FINISHED" <<endl;
	cout << "MID:" << mid <<endl;
#endif
#ifdef TIME
	gettimeofday(tp, 0);
#endif
	freeMat(ad1);
	freeMat(ad2);
	freeMat(ad3);
	freeMat(ad4);
	freeMat(ad5);
	freeMat(ad6);
	freeMat(ad7);
	freeMat(adx);
	freeMat(ady);
	//for (i = 0; i < mid; ++i) {
	//	_mm_free(m1[i]); 
	//	_mm_free(m2[i]);
	//	_mm_free(m3[i]);
	//	_mm_free(m4[i]);
	//	_mm_free(m5[i]);
	//	_mm_free(m6[i]);
	//	_mm_free(m7[i]);
	//	_mm_free(x[i]);
	//	_mm_free(y[i]);
	//}
	//_mm_free(m1); 
	//_mm_free(m2);
	//_mm_free(m3);
	//_mm_free(m4);
	//_mm_free(m5);
	//_mm_free(m6);
	//_mm_free(m7);
	//_mm_free(x);
	//_mm_free(y);

	gettimeofday(tp+1, 0);
#ifdef TIME
	mtime += (double)elapsed_time(tp);
#endif
}