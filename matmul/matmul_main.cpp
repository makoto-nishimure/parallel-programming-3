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

bool isCorrect(int size, mat a, mat b) {
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      if (a[i][j] != b[i][j]) {
        cout << "Matrix calc is wrong." << endl;
        cout << i << "," << j << endl;
        return false;
      }
    }
  }
  return true;
}

int main(int argc, char **argv) {

  int n = 1504; // Default matrix size
  int N = 1504; // Default matrix size

  int roop = 3; // Default roop num

  if (argc != 1) {
    if (1 < argc) {
      n = N = atoi(argv[1]);
      if (N % 16 != 0) {
        for (int i = N + 1;; i++) {
          if (i % 16 == 0) {
            N = i;
            cout << "N: " << N << endl;
            break;
          }
        }
      }
    }
    if (2 < argc) {
      roop = atoi(argv[2]);
    }
  }

  adds *mat1, *mat2, *mat3;
  adds *mat4;
  mat1 = allocMat(N, N);
  mat2 = allocMat(N, N);
  mat3 = allocMat(N, N);

  mat4 = allocMat(N, N);

  mat m1, m2, m3;
  mat m4;
  m1 = mat1->m;
  m2 = mat2->m;
  m3 = mat3->m;

  m4 = mat4->m;

  set_zero(N, m1);
  set_zero(N, m2);
  set_zero(N, m3);
  set_zero(N, m4);

  genmat(n, m1, 5);
  genmat(n, m2, 6);

  struct timeval tp[2];
  double one_exe_time;

  // Naive MatMul
  gettimeofday(tp, 0);
  for (int i = 0; i < roop; ++i) {
    matMul(N, m1, m2, m3);
  }
  gettimeofday(tp + 1, 0);

  one_exe_time = double(elapsed_time(tp)) / roop;
  cout << "Naive MatMul." << endl;
  printf("TIME : %lf\n\n", one_exe_time);

  // Block MatMul
  gettimeofday(tp, 0);
  for (int i = 0; i < roop; ++i) {
    block_Mul(N, m1, m2, m4);
  }
  gettimeofday(tp + 1, 0);

  if (!isCorrect(n, m3, m4)) {
    cout << "MatMul using block is not correct." << endl;
  }
  one_exe_time = double(elapsed_time(tp)) / roop;
  cout << "MatMul using block." << endl;
  printf("TIME : %lf\n\n", one_exe_time);

  // Thread Block MatMul
  set_zero(N, m4);
  gettimeofday(tp, 0);
  for (int i = 0; i < roop; ++i) {
    thr_Mul(N, m1, m2, m4);
  }
  gettimeofday(tp + 1, 0);

  if (!isCorrect(n, m3, m4)) {
    cout << "MatMul using thread and block is not correct." << endl;
    cout << m3[0][0] << endl;
    cout << m4[0][0] << endl;
  }

  one_exe_time = double(elapsed_time(tp)) / roop;
  cout << "MatMul using thread and block." << endl;
  printf("TIME : %lf\n\n", one_exe_time);

  // Strassen
  set_zero(N, m4);
  gettimeofday(tp, 0);
  for (int i = 0; i < roop; ++i) {
    str(0, N, m1, m2, m4);
  }
  gettimeofday(tp + 1, 0);

  if (!isCorrect(n, m3, m4)) {
    cout << "Strassen is not correct." << endl;
  }

  one_exe_time = double(elapsed_time(tp)) / roop;
  cout << "Strassen" << endl;
  printf("TIME : %lf\n\n", one_exe_time);

  freeMat(mat1);
  freeMat(mat2);
  freeMat(mat3);
  freeMat(mat4);

  return 0;
}
