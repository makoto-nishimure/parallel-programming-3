#include "mat_mul.hpp"
#include "mm.hpp"
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

void _genmat(int n, mat a, mat b) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      a[i][j] = ((i + j) * 100 / 7.0);
      b[i][j] = ((i + j) * 300 / 7.0);
    }
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

int main(int argc, char **argv) {
  int n = 1504;
  int N = 1504;

  adds *mat1, *mat2, *mat3;
  mat1 = allocMat(N, N);
  mat2 = allocMat(N, N);
  mat3 = allocMat(N, N);

  mat m1, m2, m3;
  m1 = mat1->m;
  m2 = mat2->m;
  m3 = mat3->m;

  set_zero(n, m3);

  struct timeval tp[2];
  int roop = 10;
  genmat(n, m1, 5);
  genmat(n, m2, 6);

  gettimeofday(tp, 0);
  for (int i = 0; i < roop; ++i) {
    // matMul(n, m1, m2, m3);
    str(0, n, m1, m2, m3);
    // thr_Mul(n, m1, m2, m3);
  }
  gettimeofday(tp + 1, 0);

  double one_exe_time = double(elapsed_time(tp)) / roop;
  printf("TIME : %lf\n", one_exe_time);

  freeMat(mat1);
  freeMat(mat2);
  freeMat(mat3);

  return 0;
}
