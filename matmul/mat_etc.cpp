#include "mat_etc.hpp"

void matAdd(int size, point pa, point pb, const mat a, const mat b, mat &c) {
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

void matSub(int size, point pa, point pb, const mat a, const mat b, mat &c) {
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

void matCopy(int size, point p, const mat from, mat &to) {
  int tmp;
  for (int i = 0; i < size; ++i) {
    tmp = i + p.first;
    for (int j = 0; j < size; ++j) {
      to[i][j] = from[tmp][j + p.second];
    }
  }
}

void set_zero(int size, mat &matrix) {
  for (int i = 0; i < size; ++i)
    for (int j = 0; j < size; ++j)
      matrix[i][j] = 0;
}
