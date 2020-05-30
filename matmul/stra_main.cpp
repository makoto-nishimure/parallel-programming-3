//#include "mm.hpp"
#include "mat_mul.hpp"
#include "mm_opt.hpp"

void genmat(int n, double **a, int k);

/* 行列生成関数 */
void genmat(int n, double **a, int k) {
  int i, j;
  double d, aij;
  int nk = n / k;
  double nk1 = 1. / nk;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      d = ((n - abs(i - j)) / (double)n) / (abs(i - j) + 1);
      aij = d * ((nk / 2 - (abs(i - j) % nk)) * nk1 + 1);
      a[i][j] = (0.1 * aij);
    }
    a[i][i] += 1;
  }
}

void genmat2vec(int n, double *a, int k) {
  int i, j;
  double d, aij;
  int nk = n / k;
  double nk1 = 1. / nk;
  for (i = n - 1, j = 0; 0 <= i; --i, ++j) {
    d = ((n - i) / (double)n) / (i + 1);
    aij = d * ((nk / 2 - (i % nk)) * nk1 + 1);
    a[j] = (0.1 * aij);
  }
  a[j - 1] += 1;
  for (i = 1, j = n; i < n; ++i, ++j) {
    d = ((n - i) / (double)n) / (i + 1);
    aij = d * ((nk / 2 - (i % nk)) * nk1 + 1);
    a[j] = (0.1 * aij);
  }
}

void _genmat(int n, mat a, mat b) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      a[i][j] = ((i + j) * 100 / 7.0);
      b[i][j] = ((i + j) * 300 / 7.0);
    }
  }
}

void _print_vec(int n, arr a) {
  int size = 2 * n - 1;
  for (int i = size - n; 0 <= i; --i) {
    for (int j = i; j < i + n; ++j) {
      printf("%.5lf ", a[j]);
    }
    cout << endl;
  }
}

void diff(int size, mat a, mat b) {
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      if (a[i][j] != b[i][j])
        cout << i << "," << j << endl;
      cout << "WORNG" << endl;
      break;
    }
  }
}

void diffmatvec(int n, mat a, arr b) {
  int size = 2 * n - 1;
  for (int i = size - n, s = 0; 0 <= i; --i, ++s) {
    for (int j = i, t = 0; j < i + n; ++j, ++t) {
      if (a[s][t] != b[j]) {
        cout << "DIFF MAT AND ARR" << endl;
        exit(1);
      }
    }
  }
}

int main(int argc, char **argv) {
  int n = 1024;
  int N = 1024;

  adds *mat1, *mat2, *mat3, *mat4, *ar1;
  mat1 = allocMat(N, N);
  ar1 = allocMat(1, 2 * N - 1);

  mat3 = allocMat(N, N);

  mat2 = allocMat(N, N);
  // mat4 = allocMat(N, N);

  mat m1, m2, m3, m4;
  m1 = mat1->m;

  m3 = mat3->m;

  m2 = mat2->m;
  // m4 = mat4->m;

  arr a1;
  a1 = ar1->a;

  set_zero(n, m3);
  // set_zero(n, m4);

  struct timeval tp[2];
  int roop = 10;
  genmat(n, m1, 5);
  genmat(n, m2, 6);
  genmat2vec(n, a1, 6);
  gettimeofday(tp, 0);
  for (int i = 0; i < roop; ++i) {
    // matMul(n, m1, m2, m3);
    // naive_ma(n, m1, a1, m3);
    // block_ma(n, m1, a1, m3);
    str(0, n, m1, m2, m3);
    // matarrMul(n, m1, a1, m3);
    // matarrMul_b(n, m1, a1, m3);
    // matarrMul_l(n, m1, a1, m3);
    // matarrMul_t(n, m1, a1, m3);
    // thr_Mul(n, m1, m2, m3);
  }
  gettimeofday(tp + 1, 0);
  // cout << "STR DONE" << endl;
  double one_exe_time = double(elapsed_time(tp)) / roop;
  printf("TIME : %lf\n", one_exe_time);
  // diff (n, m3, m4);
  // print_mat(n, m3);
#ifdef TIME
  cout << "malloc and free TIME" << endl;
  cout << mtime / roop << endl;
  cout << (mtime / roop) / one_exe_time << endl;
#endif
  freeMat(mat1);
  freeMat(mat3);

  freeMat(mat2);
  // freeMat(mat4);
  return 0;
}