#include "mat_mul.hpp"

#define BASE 1024

int block_size = 60;

adds *allocMat(int row, int column) {
  adds *p;
  p = (adds *)malloc(sizeof(adds));
  p->m = (double **)_mm_malloc(sizeof(double *) * column, 64);
  p->a = (double *)_mm_malloc(sizeof(double) * row * column, 64);
  for (int i = 0; i < column; ++i)
    p->m[i] = p->a + i * row;
  return p;
}

void freeMat(adds *p) {
  _mm_free(p->m);
  _mm_free(p->a);
  free(p);
}

void matMul(int size, const mat mat_a, const mat mat_b, mat &mat_c) {
  assert(size % 4 == 0);
  const int end = size >> 2;

  int i, j, k;
  double tmp;
  for (i = 0; i < size; ++i) {
    __m256d *c = (__m256d *)mat_c[i];
    for (k = 0; k < size; ++k) {
      __m256d *b = (__m256d *)mat_b[k];
      tmp = mat_a[i][k];
      const __m256d alpha = _mm256_set_pd(tmp, tmp, tmp, tmp);
      for (j = 0; j < end; ++j) {
        c[j] = _mm256_fmadd_pd(alpha, b[j], c[j]);
      }
    }
  }
}

void naive_ma(int size, const mat x, const arr y, mat &z) {
  assert(size % 2 == 0);
  int h_size = size >> 1;
  int arr_size = 2 * size - 1;

  int i, j, k;
  int ar_p;
  double tmp;
  for (i = 0; i < h_size; ++i) {
    for (k = 0, ar_p = arr_size - size; k < size; ++k, ar_p--) {
      double *b = y + ar_p;
      tmp = x[i][k];
      for (j = 0; j < size; ++j) {
        z[i][j] += tmp * b[j];
      }
    }
  }

  for (i = 0; i < h_size; ++i) {
    int ii = size - i - 1;
    for (j = 0; j < size; ++j) {
      z[ii][size - j - 1] = z[i][j];
    }
  }
}

void block_ma(int size, const mat x, const arr y, mat &z) {
  int h_size = size >> 1;
  // int arr_size = 2 * size - 1;
  int arr_size = (size << 1) - 1;

  int i, j, k;
  int ii, kk;
  int ar_p;
  double tmp;

  for (ii = 0; ii < size; ii += block_size) {
    int i_limit = (ii + block_size < size) ? ii + block_size : size;

    for (i = ii; i < i_limit; ++i) {
      for (k = 0, ar_p = arr_size - size; k < size; ++k, ar_p--) {
        double *b = y + ar_p;
        tmp = x[i][k];
        for (j = 0; j < size; ++j) {
          z[i][j] += tmp * b[j];
        }
      }
    }
  }
}

void matarrMul(int size, const mat x, const arr y, mat &z) {
  assert(size % 4 == 0);
  const int end = size >> 2;
  // int arr_size = 2 * size - 1;
  int arr_size = (size << 1) - 1;

  int i, j, k;
  int ar_p;
  double tmp;
  for (i = 0; i < size; ++i) {
    __m256d *c = (__m256d *)z[i];
    for (k = 0, ar_p = arr_size - size; k < size; ++k, ar_p--) {
      __m256d *b = (__m256d *)(y + ar_p);
      tmp = x[i][k];
      const __m256d alpha = _mm256_set_pd(tmp, tmp, tmp, tmp);
      for (j = 0; j < end; ++j) {
        c[j] = _mm256_fmadd_pd(alpha, b[j], c[j]);
      }
    }
  }
}

void matarrMul_t(int size, const mat a, const arr b, mat &c) {
  assert(size % 8 == 0);

  thread t1(matarrMul_tlb, size, 0, size >> 1, ref(a), ref(b), ref(c));
  matarrMul_tlb(size, size >> 1, size, ref(a), ref(b), ref(c));
  t1.join();
}

void matarrMul_tlb(int size, int thr_start, int thr_end, const mat x,
                   const arr y, mat &z) {
  assert(size % 4 == 0);
  const int end = size >> 2;
  // int arr_size = 2 * size - 1;
  int arr_size = (size << 1) - 1;

  int i, j, k;
  int ii, kk;
  int ar_p;
  double tmp;

  for (ii = thr_start; ii < thr_end; ii += block_size) {
    int i_limit = (ii + block_size < thr_end) ? ii + block_size : thr_end;
    for (kk = 0; kk < size; kk += block_size) {
      int k_limit = (kk + block_size < size) ? kk + block_size : size;

      for (i = ii; i < i_limit; ++i) {
        __m256d *c = (__m256d *)z[i];
        __m256d *d = (__m256d *)z[i + 1];

        for (k = kk, ar_p = arr_size - size - kk; k < k_limit; ++k, ar_p--) {
          __m256d *b = (__m256d *)(y + ar_p);
          tmp = x[i][k];
          const __m256d alpha = _mm256_set_pd(tmp, tmp, tmp, tmp);

          tmp = x[i + 1][k];
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

void matarrMul_tl(int size, int thr_start, int thr_end, const mat x,
                  const arr y, mat &z) {
  assert(size % 4 == 0);
  const int end = size >> 2;
  // int arr_size = 2 * size - 1;
  int arr_size = (size << 1) - 1;

  int i, j, k;
  int ar_p;
  double tmp;
  for (i = thr_start; i < thr_end; ++i) {
    __m256d *c = (__m256d *)z[i];
    __m256d *d = (__m256d *)z[i + 1];

    for (k = 0, ar_p = arr_size - size; k < size; ++k, ar_p--) {
      __m256d *b = (__m256d *)(y + ar_p);
      tmp = x[i][k];
      const __m256d alpha = _mm256_set_pd(tmp, tmp, tmp, tmp);

      tmp = x[i + 1][k];
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

void matarrMul_l(int size, const mat x, const arr y, mat &z) {
  assert(size % 4 == 0);
  const int end = size >> 2;
  int arr_size = (size << 1) - 1;

  int i, j, k;
  int ar_p;
  double tmp;
  for (i = 0; i < size; ++i) {
    __m256d *c = (__m256d *)z[i];
    __m256d *d = (__m256d *)z[i + 1];

    for (k = 0, ar_p = arr_size - size; k < size; ++k, ar_p--) {
      __m256d *b = (__m256d *)(y + ar_p);
      tmp = x[i][k];
      const __m256d alpha = _mm256_set_pd(tmp, tmp, tmp, tmp);

      tmp = x[i + 1][k];
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

void matarrMul_b(int size, const mat x, const arr y, mat &z) {
  assert(size % 4 == 0);
  const int end = size >> 2;
  // int arr_size = 2 * size - 1;
  int arr_size = (size << 1) - 1;

  int i, j, k;
  int ii, kk;
  int ar_p;
  double tmp;

  for (ii = 0; ii < size; ii += block_size) {
    int i_limit = (ii + block_size < size) ? ii + block_size : size;
    for (i = ii; i < i_limit; ++i) {
      __m256d *c = (__m256d *)z[i];
      for (k = 0, ar_p = arr_size - size; k < size; ++k, ar_p--) {
        __m256d *b = (__m256d *)(y + ar_p);
        tmp = x[i][k];
        const __m256d alpha = _mm256_set_pd(tmp, tmp, tmp, tmp);
        for (j = 0; j < end; ++j) {
          c[j] = _mm256_fmadd_pd(alpha, b[j], c[j]);
        }
      }
    }
  }
}

void block_Mul(int size, const mat mat_a, const mat mat_b, mat &mat_c) {
  assert(size % 4 == 0);
  const int end = size >> 2;

  int i, j, k;
  int ii, kk;
  double tmp;
  for (ii = 0; ii < size; ii += block_size) {
    int i_limit = (ii + block_size < size) ? ii + block_size : size;
    for (kk = 0; kk < size; kk += block_size) {
      int k_limit = (kk + block_size < size) ? kk + block_size : size;

      for (i = ii; i < i_limit; ++i) {
        __m256d *c = (__m256d *)mat_c[i];
        for (k = kk; k < k_limit; ++k) {
          __m256d *b = (__m256d *)mat_b[k];
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

void tb_Mul(int size, int thr_s, int thr_e, const mat mat_a, const mat mat_b,
            mat &mat_c) {
  assert(size % 4 == 0);
  const int end = size >> 2;

  int i, j, k;
  int ii, kk;
  double tmp;
  for (ii = thr_s; ii < thr_e; ii += block_size) {
    int i_limit = (ii + block_size < thr_e) ? ii + block_size : thr_e;
    for (kk = 0; kk < size; kk += block_size) {
      int k_limit = (kk + block_size < size) ? kk + block_size : size;

      for (i = ii; i < i_limit; ++i) {
        __m256d *c = (__m256d *)mat_c[i];
        __m256d *d = (__m256d *)mat_c[i + 1];

        for (k = kk; k < k_limit; ++k) {
          __m256d *b = (__m256d *)mat_b[k];

          tmp = mat_a[i][k];
          const __m256d alpha = _mm256_set_pd(tmp, tmp, tmp, tmp);

          tmp = mat_a[i + 1][k];
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

void thr_Mul(int size, const mat a, const mat b, mat &c) {
  assert(size % 8 == 0);
  thread t1(tb_Mul, size, 0, size >> 1, ref(a), ref(b), ref(c));
  tb_Mul(size, size >> 1, size, ref(a), ref(b), ref(c));
  t1.join();
}

void str(int d, int size, const mat mat_a, const mat mat_b, mat &mat_c) {
  struct timeval tp[2];
#ifdef DEBUG
  for (int i = 0; i < d; ++i)
    cout << "\t";
  cout << "STR START. SIZE:" << size << endl;
#endif
  if (size <= BASE) {
    // matMul(size, mat_a, mat_b, mat_c);
    // block_Mul(size, mat_a, mat_b, mat_c);
    thr_Mul(size, mat_a, mat_b, mat_c);
#ifdef DEBUG
    for (int i = 0; i < d; ++i)
      cout << "\t";
    cout << "matMul DONE. SIZE:" << size << endl;
#endif
    return;
  }
  int i, j, k, mid;
  mid = size >> 1;

  point p0, p1, p2, p3;
  p0 = make_pair(0, 0);
  p1 = make_pair(0, mid);
  p2 = make_pair(mid, 0);
  p3 = make_pair(mid, mid);
  mat m1, m2, m3, m4, m5, m6, m7, x, y;
  adds *ad1, *ad2, *ad3, *ad4, *ad5, *ad6, *ad7, *adx, *ady;

  ad1 = allocMat(mid, mid);
  ad2 = allocMat(mid, mid);
  ad3 = allocMat(mid, mid);
  ad4 = allocMat(mid, mid);
  ad5 = allocMat(mid, mid);
  ad6 = allocMat(mid, mid);
  ad7 = allocMat(mid, mid);
  adx = allocMat(mid, mid);
  ady = allocMat(mid, mid);
  m1 = ad1->m;
  m2 = ad2->m;
  m3 = ad3->m;
  m4 = ad4->m;
  m5 = ad5->m;
  m6 = ad6->m;
  m7 = ad7->m;
  x = adx->m;
  y = ady->m;

  // M1
  matAdd(mid, p0, p3, mat_a, mat_a, x);
  matAdd(mid, p0, p3, mat_b, mat_b, y);
  set_zero(mid, m1);
  str(d + 1, mid, x, y, m1);
#ifdef DEBUG
  for (int i = 0; i < d; ++i)
    cout << "\t";
  cout << "M1 FINISHED" << endl;
#endif

  // M2
  matAdd(mid, p2, p3, mat_a, mat_a, x);
  matCopy(mid, p0, mat_b, y);
  set_zero(mid, m2);
  str(d + 1, mid, x, y, m2);
#ifdef DEBUG
  for (int i = 0; i < d; ++i)
    cout << "\t";
  cout << "M2 FINISHED" << endl;
#endif

  // M3
  matCopy(mid, p0, mat_a, x);
  matSub(mid, p1, p3, mat_b, mat_b, y);
  set_zero(mid, m3);
  str(d + 1, mid, x, y, m3);
#ifdef DEBUG
  for (int i = 0; i < d; ++i)
    cout << "\t";
  cout << "M3 FINISHED" << endl;
#endif

  // M4
  matCopy(mid, p3, mat_a, x);
  matSub(mid, p2, p0, mat_b, mat_b, y);
  set_zero(mid, m4);
  str(d + 1, mid, x, y, m4);
#ifdef DEBUG
  for (int i = 0; i < d; ++i)
    cout << "\t";
  cout << "M4 FINISHED" << endl;
#endif

  // M5
  matAdd(mid, p0, p1, mat_a, mat_a, x);
  matCopy(mid, p3, mat_b, y);
  set_zero(mid, m5);
  str(d + 1, mid, x, y, m5);
#ifdef DEBUG
  for (int i = 0; i < d; ++i)
    cout << "\t";
  cout << "M5 FINISHED" << endl;
#endif

  // M6
  matSub(mid, p2, p0, mat_a, mat_a, x);
  matAdd(mid, p0, p1, mat_b, mat_b, y);
  set_zero(mid, m6);
  str(d + 1, mid, x, y, m6);
#ifdef DEBUG
  for (int i = 0; i < d; ++i)
    cout << "\t";
  cout << "M6 FINISHED" << endl;
#endif

  // M7
  matSub(mid, p1, p3, mat_a, mat_a, x);
  matAdd(mid, p2, p3, mat_b, mat_b, y);
  set_zero(mid, m7);
  str(d + 1, mid, x, y, m7);
#ifdef DEBUG
  for (int i = 0; i < d; ++i)
    cout << "\t";
  cout << "M7 FINISHED" << endl;
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

#ifdef DEBUG
  for (int i = 0; i < d; ++i)
    cout << "\t";
  cout << "C1 AND C2  FINISHED" << endl;
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
  cout << "C3 AND C4  FINISHED" << endl;
  cout << "MID:" << mid << endl;
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
}
