#include "mm.hpp"
//#include "mat_mul.hpp"

int block_size = 60;

typedef struct {
	mat m;
	arr a;
} adds;

adds* allocMat(int, int);
void freeMat(adds *p);

double elapsed_time(struct timeval tp[2])
{
	return tp[1].tv_sec-tp[0].tv_sec+1e-6*(tp[1].tv_usec-tp[0].tv_usec);
}

void print_mat(int size, mat m)
{
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			printf("%.5lf ", m[i][j]);
		}
		cout << endl;
	}
}

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
	_mm_free(p->m);
	_mm_free(p->a);
	free(p);
}

void _print_vec(int n, arr a)
{
  int size = 2 * n - 1;
  for (int i = size - n; 0 <= i; --i) {
    for (int j = i; j < i + n; ++j) {
      if (a[j] >= 0)
        printf(" %.5lf ", a[j]);
      else 
        printf("%.5lf ", a[j]);
    }
    cout << endl;
  }
}

/* 行列生成関数 */
void
genmat(int n, double** a, int k)
{
  int i,j;
  double d, aij;
  int nk = n/k;
  double nk1 = 1./nk;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      d = ((n - abs(i-j))/(double)n) / (abs(i - j) + 1);
      aij = d * ( (nk/2 - (abs(i - j) % nk))*nk1 + 1);
      a[i][j] = (0.1 * aij);
    }
    a[i][i] += 1;
  }
}

  void
genmat2vec(int n, double* a, int k)
{
  int i,j;
  double d, aij;
  int nk = n/k;
  double nk1 = 1./nk;
  for (i = n-1, j = 0; 0 <= i; --i, ++j) {
    d = ((n - i)/(double)n) / (i + 1);
    aij = d * ( (nk/2 - (i % nk))*nk1 + 1);
    a[j] = (0.1 * aij);
  }
  a[j-1] += 1;
  for (i = 1, j = n; i < n; ++i, ++j) {
    d = ((n - i)/(double)n) / (i + 1);
    aij = d * ( (nk/2 - (i % nk))*nk1 + 1);
    a[j] = (0.1 * aij);
  }
}

// length of B is 2 * n - 1
// length of B0 is n - 1
void genB01 (int n, arr B, arr B0, arr B1)
{
  double tmp;
  for (int i = 0, start = n >> 1; i < n - 1; ++i, ++start) {
    tmp = B[start];
    B0[i] = 2 * tmp;
    B1[i] = tmp;
  }
}

void genB24 (int n, arr B, arr B2, arr B4)
{
  for (int i = 0, B1s = n >> 1, Bs = n; i < n - 1; ++i, ++B1s, ++Bs) {
    B2[i] = B[Bs] - B[B1s];
    B4[i] = B[Bs] + B[B1s];
  }
}

void genB35 (int n, arr B, arr B3, arr B5)
{
  for (int i = 0, B1s = n >> 1; i < n - 1; ++i, ++B1s) {
    B3[i] = B[i] - B[B1s];
    B5[i] = B[i] + B[B1s];
  }
}

void genB (int n, arr B, adds *vec)
{
  mat m = vec->m;
  arr v0 = m[0];
  arr v1 = m[1];
  arr v2 = m[2];
  arr v3 = m[3];
  arr v4 = m[4];
  arr v5 = m[5];

  genB01(n, B, v0, v1);
  genB24(n, B, v2, v4);
  genB35(n, B, v3, v5);
}

void matAdd(int size, point pa, point pb,
		const mat a, const mat b, mat& c)
{
	int i, j;
	int tmp_pa, tmp_pb;
	for (i = 0; i < size; ++i) {
		tmp_pa = pa.first + i;
		tmp_pb = pb.first + i;
		for (j = 0; j < size; ++j) {
			c[i][j] = a[tmp_pa][pa.second + j] + b[tmp_pb][pb.second + j];
		}
	}
}

void matSub(int size, point pa, point pb,
		const mat a, const mat b, mat& c)
{
	int i, j;
	int tmp_pa, tmp_pb;
	for (i = 0; i < size; ++i) {
		tmp_pa = pa.first + i;
		tmp_pb = pb.first + i;
		for (j = 0; j < size; ++j) {
			c[i][j] = a[tmp_pa][pa.second + j] - b[tmp_pb][pb.second + j];
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

void matarrMul_tlb(int size, int thr_start, int thr_end,
    const mat x, const arr y, mat& z)
{
  assert(size % 4 == 0);
  const int end = size >> 2;
  //int arr_size = 2 * size - 1;
  int arr_size = (size << 1) - 1;

  int i, j, k;
  int ii, kk;
  int ar_p;
  double tmp;

  for (ii = thr_start; ii < thr_end; ii+=block_size) {
    int i_limit = (ii + block_size < thr_end) ? ii + block_size : thr_end;
    for (kk = 0; kk < size; kk+=block_size) {
      int k_limit = (kk + block_size < size) ? kk + block_size : size;
      
      for (i = ii; i < i_limit; ++i) {
        __m256d* c = (__m256d*)z[i];
        __m256d* d = (__m256d*)z[i+1];

        for (k = kk, ar_p = arr_size - size - kk; k < k_limit; ++k, ar_p--) {
          __m256d* b = (__m256d*)(y + ar_p);
          tmp = x[i][k];
          const __m256d alpha = _mm256_set_pd(tmp, tmp, tmp, tmp);

          tmp = x[i+1][k];
          const __m256d beta = _mm256_set_pd(tmp, tmp, tmp, tmp);

          for (j = 0; j < end; ++j) {
            c[j] = _mm256_fmadd_pd(alpha, b[j], c[j]);
            d[j] = _mm256_fmadd_pd(beta, b[j], d[j]);
            ++j;
            c[j] = _mm256_fmadd_pd(alpha, b[j], c[j]);
            d[j] = _mm256_fmadd_pd(beta, b[j], d[j]);
          }
        }
        ++i;
      }
    }
  }
}

void matarrMul_t(int size,
    const mat a, const arr b, mat& c)
{
  assert(size % 8 == 0);
  
  thread t1(matarrMul_tlb, size, 0, size >> 1, ref(a), ref(b), ref(c));
  matarrMul_tlb(size, size >> 1, size, ref(a), ref(b), ref(c));
  t1.join();
}


void ma_str (int size, const mat ma, const mat vecs, mat& mc,
    mat* m, mat w_area)
{
  //matarrMul_t(int size, const mat x, const arr y, mat& z);
  int i, j, k, mid;
  mid = size >> 1;

  point p0, p1, p2, p3;
  p0 = make_pair(0,0);
  p1 = make_pair(0,mid);
  p2 = make_pair(mid,0);
  p3 = make_pair(mid,mid);

  // M0
  matAdd(mid, p0, p3, ma, ma, w_area);
  matarrMul_t(mid, w_area, vecs[0], m[0]);

  // M1
  matAdd(mid, p2, p3, ma, ma, w_area);
  matarrMul_t(mid, w_area, vecs[1], m[1]);

  // M2
  matCopy(mid, p0, ma, w_area);
  matarrMul_t(mid, w_area, vecs[2], m[2]);

  // M3
  matCopy(mid, p3, ma, w_area);
  matarrMul_t(mid, w_area, vecs[3], m[3]);

  // M4
  matAdd(mid, p0, p1, ma, ma, w_area);
  matarrMul_t(mid, w_area, vecs[1], m[4]);

  // M5
  matSub(mid, p2, p0, ma, ma, w_area);
  matarrMul_t(mid, w_area, vecs[4], m[5]);

  // M6
  matSub(mid, p1, p3, ma, ma, w_area);
  matarrMul_t(mid, w_area, vecs[5], m[6]);

  // C1 and C2
  // C1
  int tmp_n, tmp_m;
  matAdd(mid, p0, p0, m[0], m[3], w_area);
  matSub(mid, p0, p0, w_area, m[4], w_area);
  matAdd(mid, p0, p0, w_area, m[6], w_area);
  for (i = 0; i < mid; ++i)
    for (j = 0; j < mid; ++j)
      mc[i][j] = w_area[i][j];
  // C2
  matAdd(mid, p0, p0, m[2], m[4], w_area);
  for (i = 0; i < mid; ++i)
    for (j = 0, tmp_m = mid; j < mid; ++j, ++tmp_m)
      mc[i][tmp_m] = w_area[i][j];

  // C3 and C4
  // C3
  matAdd(mid, p0, p0, m[1], m[3], w_area);
  for (i = 0, tmp_n = mid; i < mid; ++i, ++tmp_n) 
    for (j = 0; j < mid; ++j)
      mc[tmp_n][j] = w_area[i][j];

  // C4
  matSub(mid, p0, p0, m[0], m[1], w_area);
  matAdd(mid, p0, p0, w_area, m[2], w_area);
  matAdd(mid, p0, p0, w_area, m[5], w_area);
  for (i = 0, tmp_n = mid; i < mid; ++i, ++tmp_n)
    for (j = 0, tmp_m = mid; j < mid; ++j, ++tmp_m)
      mc[tmp_n][tmp_m] = w_area[i][j];
}

void set_zero(int size, mat& matrix)
{
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			matrix[i][j] = 0;
}

int main (void)
{
  int N = 1504;
  int n = 1504;
  int mid = n >> 1;

  // for matrix
  adds *a, *c;
  mat ma, mc;
  a = allocMat(n, n);
  c = allocMat(n, n);
  ma = a->m;
  mc = c->m;
  //

  // for B
  adds *ar1, *vec;
  ar1  = allocMat(2 * n - 1, 1);
  arr a1 = ar1->a;
  vec  = allocMat(n - 1, 6);
  //

  // for str
  adds **matrix, *work_area;
  matrix = (adds**)malloc(sizeof(adds*) * 7);
  mat *m, w_area;
  m = (mat*)malloc(sizeof(mat) * 7);
  for (int i = 0; i < 7; ++i) {
    matrix[i] = allocMat(mid, mid);
    m[i] = matrix[i]->m;
  }
  work_area = allocMat(mid, mid);
  w_area = work_area->m;
  //


  set_zero(n, mc);
  for (int i = 0; i < 7; ++i)
    set_zero(mid, m[i]);

  genmat(n, ma, 5);
  genmat2vec(n, a1, 6);
  //_print_vec(n, a1);
  //cout << endl;
  genB(n, a1, vec);
  struct timeval tp[2];

  gettimeofday(tp, 0);
  ma_str (n, ma, vec->m, mc, m, w_area);
  gettimeofday(tp+1, 0);
  double one_exe_time = double(elapsed_time(tp));
  printf("TIME : %lf\n", one_exe_time);

  print_mat(n, mc);

//void ma_str (int size, const mat ma, const mat vecs, mat& mc,
//    mat* m, mat w_area)

  free(ar1);
  free(vec);
}
