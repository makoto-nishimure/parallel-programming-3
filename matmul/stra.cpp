#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <sys/time.h>
#include <utility>
#include <vector>

using namespace std;
typedef vector<double> arr;
typedef vector<arr> mat;
typedef pair<int, int> point;

#define BASE 512

double elapsed_time(struct timeval tp[2]) {
  return tp[1].tv_sec - tp[0].tv_sec + 1e-6 * (tp[1].tv_usec - tp[0].tv_usec);
}

void print_pair(int size, point pa) {
  cout << pa.first << "," << pa.second << " to " << pa.first + size - 1 << ","
       << pa.second + size - 1 << endl;
}

void matAdd(int size, point pa, point pb, mat &mat_a, mat &mat_b, mat &mat_c) {
#ifdef DEBUG
  print_pair(size, pa);
  print_pair(size, pb);
#endif
  int i, j;
  int tmp_pa, tmp_pb;
  for (i = 0; i < size; i++) {
    tmp_pa = pa.first + i;
    tmp_pb = pb.first + i;
    for (j = 0; j < size; j++)
      mat_c[i][j] = mat_a[tmp_pa][j + pa.second] + mat_b[tmp_pb][j + pb.second];
  }
}

void matSub(int size, point pa, point pb, mat &mat_a, mat &mat_b, mat &mat_c) {
#ifdef DEBUG
  print_pair(size, pa);
  print_pair(size, pb);
#endif
  int i, j;
  int tmp_pa, tmp_pb;
  for (i = 0; i < size; i++) {
    tmp_pa = pa.first + i;
    tmp_pb = pb.first + i;
    for (j = 0; j < size; j++)
      mat_c[i][j] = mat_a[tmp_pa][j + pa.second] - mat_b[tmp_pb][j + pb.second];
  }
}

void matCopy(int size, point p, mat &from, mat &to) {
  int tmp;
  for (int i = 0; i < size; i++) {
    tmp = i + p.first;
    for (int j = 0; j < size; j++)
      to[i][j] = from[tmp][j + p.second];
  }
}

void matMul(int size, mat &mat_a, mat &mat_b, mat &mat_c) {
  int i, j, k;
  double tmp;
  for (i = 0; i < size; i++) {
    for (k = 0; k < size; k++) {
      tmp = mat_a[i][k];
      for (j = 0; j < size; j++) {
        mat_c[i][j] += tmp * mat_b[k][j];
      }
    }
  }
}

void str(int d, int size, mat &mat_a, mat &mat_b, mat &mat_c) {
#ifdef DEBUG
  for (int i = 0; i < d; i++)
    cout << "\t";
  cout << "STR START. SIZE:" << size << endl;
#endif
  if (size <= BASE) {
    matMul(size, mat_a, mat_b, mat_c);
#ifdef DEBUG
    for (int i = 0; i < d; i++)
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
  for (i = 0; i < mid; i++) {
    m1.push_back(arr(mid));
    m2.push_back(arr(mid));
    m3.push_back(arr(mid));
    m4.push_back(arr(mid));
    m5.push_back(arr(mid));
    m6.push_back(arr(mid));
    m7.push_back(arr(mid));
    x.push_back(arr(mid));
    y.push_back(arr(mid));
  }

  // M1
  matAdd(mid, p0, p3, mat_a, mat_a, x);
  matAdd(mid, p0, p3, mat_b, mat_b, y);
  str(d + 1, mid, x, y, m1);
#ifdef DEBUG
  for (int i = 0; i < d; i++)
    cout << "\t";
  cout << "M1 FINISHED" << endl;
#endif

  // M2
  matAdd(mid, p2, p3, mat_a, mat_a, x);
  matCopy(mid, p0, mat_b, y);
  str(d + 1, mid, x, y, m2);
#ifdef DEBUG
  for (int i = 0; i < d; i++)
    cout << "\t";
  cout << "M2 FINISHED" << endl;
#endif

  // M3
  matCopy(mid, p0, mat_a, x);
  matSub(mid, p1, p3, mat_b, mat_b, y);
  str(d + 1, mid, x, y, m3);
#ifdef DEBUG
  for (int i = 0; i < d; i++)
    cout << "\t";
  cout << "M3 FINISHED" << endl;
#endif

  // M4
  matCopy(mid, p3, mat_a, x);
  matSub(mid, p2, p0, mat_b, mat_b, y);
  str(d + 1, mid, x, y, m4);
#ifdef DEBUG
  for (int i = 0; i < d; i++)
    cout << "\t";
  cout << "M4 FINISHED" << endl;
#endif

  // M5
  matAdd(mid, p0, p1, mat_a, mat_a, x);
  matCopy(mid, p3, mat_b, y);
  str(d + 1, mid, x, y, m5);
#ifdef DEBUG
  for (int i = 0; i < d; i++)
    cout << "\t";
  cout << "M5 FINISHED" << endl;
#endif

  // M6
  matSub(mid, p2, p0, mat_a, mat_a, x);
  matAdd(mid, p0, p1, mat_b, mat_b, y);
  str(d + 1, mid, x, y, m6);
#ifdef DEBUG
  for (int i = 0; i < d; i++)
    cout << "\t";
  cout << "M6 FINISHED" << endl;
#endif

  // M7
  matSub(mid, p1, p3, mat_a, mat_a, x);
  matAdd(mid, p2, p3, mat_b, mat_b, y);
  str(d + 1, mid, x, y, m7);
#ifdef DEBUG
  for (int i = 0; i < d; i++)
    cout << "\t";
  cout << "M7 FINISHED" << endl;
#endif

  // C1 and C2
  int tmp_n, tmp_m;
  for (i = 0; i < mid; i++) {
    for (j = 0; j < mid; j++)
      mat_c[i][j] = m1[i][j] + m4[i][j] - m5[i][j] + m7[i][j];
    for (j = 0, tmp_m = mid; j < mid; j++, tmp_m++)
      mat_c[i][tmp_m] = m3[i][j] + m5[i][j];
  }
#ifdef DEBUG
  for (int i = 0; i < d; i++)
    cout << "\t";
  cout << "C1 AND C2  FINISHED" << endl;
#endif
  // C3 and C4
  for (i = 0, tmp_n = mid; i < mid; i++, tmp_n++) {
    for (j = 0; j < mid; j++)
      mat_c[tmp_n][j] = m2[i][j] + m4[i][j];
    for (j = 0, tmp_m = mid; j < mid; j++, tmp_m++)
      mat_c[tmp_n][tmp_m] = m1[i][j] - m2[i][j] + m3[i][j] + m6[i][j];
  }
#ifdef DEBUG
  for (int i = 0; i < d; i++)
    cout << "\t";
  cout << "C3 AND C4  FINISHED" << endl;
#endif
}

void print_mat(mat m) {
  for (int i = 0; i < m.size(); i++) {
    for (int j = 0; j < m[i].size(); j++) {
      cout << m[i][j] << " ";
    }
    cout << endl;
  }
}

int main(int argc, char **argv) {
  mat m1, m2, m3;
  mat m4;
  int n = 1600;
  int N = 1600;
  for (int i = 0; i < N; i++) {
    m1.push_back(arr(N));
    m2.push_back(arr(N));
    m3.push_back(arr(N));
    // m4.push_back(arr(N));
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      m1[i][j] = ((i + j) * 100 / 7.0);
      m2[i][j] = ((i + j) * 300 / 7.0);
    }
  }

  struct timeval tp[2];
  gettimeofday(tp, 0);
  str(0, n, m1, m2, m3);
  gettimeofday(tp + 1, 0);
  // cout << double(elapsed_time(tp)) <<endl;;
  print_mat(m3);

  // gettimeofday(tp, 0);
  // matMul(n, m1, m2, m4);
  // gettimeofday(tp+1, 0);
  // cout << double(elapsed_time(tp)) <<endl;;
  return 0;
}
