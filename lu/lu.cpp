#include "lu.hpp"
int lu_block_size = 60;

void _print_mat(int n, mat a) {
  for (int i = 0; i < n; ++i) {
    if (i % 4 == 0 && i != 0) {
      for (int k = 0; k < 150; ++k)
        cout << "-";
      cout << endl;
    }
    for (int j = 0; j < n; ++j) {
      if (j % 4 == 0 && j != 0)
        cout << " | ";
      printf(" %6.2f ", a[i][j]);
    }
    cout << endl;
  }
}

void naive(int n, mat a) {
  int i, j, k;
  for (i = 0; i < n; i++) {
    for (k = 0; k < i; k++)
      for (j = k + 1; j < n; j++)
        a[i][j] -= a[i][k] * a[k][j];
    // reciprocal
    double rep = a[i][i];
    for (j = i + 1; j < n; j++)
      a[i][j] /= a[i][i];
  }
}

void init_vec(int n, double *vec) {
  for (int i = 0; i < n; ++i)
    vec[i] = 0.0;
}

void block(int n, point p, mat a) {
#ifdef DEBUG
  cout << "BLOCK" << endl;
  cout << "N:" << n << endl;
  cout << "LU_BLOCK_SIZE:" << lu_block_size << endl;
  cout << "POINT : P(" << p.first << "," << p.second << ")" << endl;
  cout << endl;
  assert(p.first + lu_block_size <= n);
  assert(p.second + lu_block_size <= n);
#endif
  int cLimit = (p.first + lu_block_size < n) ? p.first + lu_block_size : n;
  int rLimit = (p.second + lu_block_size < n) ? p.second + lu_block_size : n;
  for (int i = p.first; i < cLimit; ++i) {
    for (int k = p.second; k < i; ++k) {
      for (int j = k + 1; j < rLimit; ++j) {
        a[i][j] -= a[i][k] * a[k][j];
      }
    }
    for (int j = i + 1; j < rLimit; ++j) {
      a[i][j] /= a[i][i];
    }
    //	_print_mat(n, a);
    //	cout << endl;
  }
}

void calc_L_0(int n, point start, mat matrix, mat lu)
// void calc_L_0 (int n, point start, double matrix[3][3], double lu[3][3])
{
#ifdef DEBUG
  cout << "calc_L" << endl;
  cout << "START(" << start.first << "," << start.second << ")" << endl;
#endif
  // loop with column
  for (int i = 0; i < lu_block_size; ++i) {
    int c_i = start.first + i;

    double vec[lu_block_size];
    init_vec(lu_block_size, vec);
    vec[0] = matrix[c_i][start.second];
    int vecLen = 1;
    // loop with row
    for (int j = 1; j < lu_block_size; ++j) {
      int r_j = start.second + j;

      double elm = matrix[c_i][r_j];
      for (int k = 0; k < vecLen; ++k) {
        //#ifdef DEBUG
        //				cout << "ELM -= vec * lu" << endl;
        //				cout << elm << " -= " << vec[k] << " * "
        //<< lu[k][j] << endl; #endif
        elm -= vec[k] * lu[k][j];
      }
      vec[vecLen] = elm;
      ++vecLen;
    }
    for (int j = 0; j < lu_block_size; ++j) {
      matrix[c_i][start.second + j] = vec[j];
    }
  }
}

void calc_L(int n, int start, mat matrix, mat lu) {
  for (int i = start + lu_block_size; i < n; i += lu_block_size) {
    calc_L_0(n, make_pair(i, start), matrix, lu);
  }
#ifdef DEBUG
  _print_mat(n, matrix);
  cout << endl;
#endif
}

void calc_U_0(int n, point start, mat matrix, mat lu) {
#ifdef DEBUG
  cout << "calc_U" << endl;
  cout << "START(" << start.first << "," << start.second << ")" << endl;
#endif
  double tmp = lu[0][0];
  int i_lim = start.second + lu_block_size;
  for (int i = start.second; i < i_lim; ++i)
    matrix[start.first][i] /= tmp;

  for (int i = 1; i < lu_block_size; ++i) {
    int c_i = start.first + i;

    double vec[lu_block_size];
    init_vec(lu_block_size, vec);

    for (int j = 0; j < i; ++j) {
      int r_j = start.second + j;
      double scala = lu[i][j];

      for (int k = 0; k < lu_block_size; ++k) {
        // vec[k] += scala * matrix[c_i-1][start.second + k];
        vec[k] += scala * matrix[start.first + j][start.second + k];
      }
    }

    double scala = lu[i][i];
    for (int j = 0; j < lu_block_size; ++j) {
      int r_j = start.second + j;
      matrix[c_i][r_j] = (matrix[c_i][r_j] - vec[j]) / scala;
    }
  }
}

void calc_U(int n, int start, mat matrix, mat lu) {
  for (int i = start + lu_block_size; i < n; i += lu_block_size) {
    calc_U_0(n, make_pair(start, i), matrix, lu);
  }
#ifdef DEBUG
  _print_mat(n, matrix);
  cout << endl;
#endif
}

void matSub_L(int n, int start, int numBlock, mat matrix) {
#ifdef DEBUG
  cout << "matSub_L" << endl;
  cout << "START(" << start << "," << start << ")" << endl;
#endif
  adds *adc;
  adc = allocMat(lu_block_size, lu_block_size);

  // loop with column
  for (int c = start; c < n; c += lu_block_size) {
    int x = start;
    int y = start - lu_block_size;
    matCopy(lu_block_size, make_pair(x, y), matrix, adc->m);
    mat tmp = adc->m;
    // loop with row
    for (int r = start; r < n; r += lu_block_size) {
      // loop with block_i
      for (int i = 0; i < numBlock; ++i) {
        int cur_i = c + i;
        // loop with block_j
        for (int j = i; j <= i; ++j) {
          matrix[cur_i][r + j] -= tmp[i][j];
        }
      }
    }
  }
  freeMat(adc);
}

void _matMul(int size, mat l, mat u, mat lu) {
  double tmp;
  for (int i = 0; i < size; i++) {
    for (int k = 0; k < size; k++) {
      tmp = l[i][k];
      for (int j = 0; j < size; j++) {
        lu[i][j] += tmp * u[k][j];
      }
    }
  }
}

void _matSub(int n, point p, mat matrix, mat lu) {
  for (int i = 0; i < n; ++i) {
    int c_i = p.first + i;
    for (int j = 0; j < n; ++j) {
      matrix[c_i][p.second + j] -= lu[i][j];
    }
  }
}

void MatSub_0(int n, point start, int numBlock, mat matrix, adds **us) {
#ifdef DEBUG
  cout << "MatSub_0 : (" << start.first << "," << start.second << ")" << endl;
#endif
  adds *l, *lu;
  l = allocMat(lu_block_size, lu_block_size);
  lu = allocMat(lu_block_size, lu_block_size);
  int x = start.first;
  int y = start.second - lu_block_size;
  matCopy(lu_block_size, make_pair(x, y), matrix, l->m);

  for (int i = 0; i < numBlock; ++i) {
    set_zero(lu_block_size, lu->m);
    //_matMul(lu_block_size, l->m, us[i]->m, lu->m);
    matMul(lu_block_size, l->m, us[i]->m, lu->m);

    int r_i = start.second + lu_block_size * i;
    _matSub(lu_block_size, make_pair(start.first, r_i), matrix, lu->m);
  }
  freeMat(l);
  freeMat(lu);
}

void matMS(int size, point start, mat matrix, mat l, mat u) {
  assert(size % 4 == 0);
  const int end = size >> 2;

  int i, j, k;
  double tmp;
  for (i = 0; i < size; ++i) {
    __m256d *c = (__m256d *)matrix[start.first + i];
    for (k = 0; k < size; ++k) {
      __m256d *b = (__m256d *)u[k];
      tmp = l[i][k];
      const __m256d alpha = _mm256_set_pd(tmp, tmp, tmp, tmp);
      for (j = 0; j < end; ++j) {
        c[j] = _mm256_fmsub_pd(alpha, b[j], c[j]);
      }
    }
  }
}

void MatSub_1(int n, point start, int numBlock, mat matrix, adds **us) {
#ifdef DEBUG
  cout << "MatSub_0 : (" << start.first << "," << start.second << ")" << endl;
#endif
  adds *l;
  l = allocMat(lu_block_size, lu_block_size);
  int x = start.first;
  int y = start.second - lu_block_size;
  matCopy(lu_block_size, make_pair(x, y), matrix, l->m);

  for (int i = 0; i < numBlock; ++i) {
    int r_i = start.second + lu_block_size * i;
    matMS(lu_block_size, make_pair(start.first, r_i), matrix, l->m, us[i]->m);
  }
  freeMat(l);
}

void MatSub(int n, int start, mat matrix) {
#ifdef DEBUG
  cout << "MatSub" << endl;
#endif
  int numBlock = (n - start) / lu_block_size;
  adds **us;
  us = (adds **)malloc(sizeof(adds *) * numBlock);
  for (int i = 0; i < numBlock; ++i) {
    us[i] = allocMat(lu_block_size, lu_block_size);
  }

  // Copy matblocks
  for (int i = 0; i < numBlock; ++i) {
    int x = start - lu_block_size;
    int y = start + (i * lu_block_size);
    matCopy(lu_block_size, make_pair(x, y), matrix, us[i]->m);
  }

  for (int i = 0; i < numBlock; ++i) {
    int x = start + i * lu_block_size;
    int y = start;
    MatSub_0(n, make_pair(x, y), numBlock, matrix, us);
    // MatSub_1(n, make_pair(x, y), numBlock, matrix, us);
  }

  for (int i = 0; i < numBlock; ++i) {
    freeMat(us[i]);
  }
  free(us);
}

// int start; �Ϲ���κ���Υݥ����
void block_gauss_0(int n, int start, mat matrix) {
#ifdef DEBUG
  cout << "gauss_0" << endl;
  cout << "START(" << start << "," << start << ")" << endl;
#endif
  point p = make_pair(start, start);
  MatSub(n, start, matrix);
#ifdef DEBUG
  _print_mat(n, matrix);
  cout << endl;
#endif
  // LUʬ��
  block(n, p, matrix);
#ifdef DEBUG
  _print_mat(n, matrix);
  cout << endl;
#endif

  if (start + lu_block_size <= n) {
    adds *lu;
    lu = allocMat(lu_block_size, lu_block_size);
    matCopy(lu_block_size, p, matrix, lu->m);

    // thread
    calc_U(n, start, matrix, lu->m);
    // thread
    calc_L(n, start, matrix, lu->m);

    freeMat(lu);
  }
}

// �ǽ�Υ֥��å��Τ�
void block_gauss(int n, mat matrix) {
  point p = make_pair(0, 0);
  block(n, p, matrix);

  adds *lu;
  lu = allocMat(lu_block_size, lu_block_size);
  matCopy(lu_block_size, p, matrix, lu->m);
  // thread
  calc_U(n, 0, matrix, lu->m);
  // thread
  calc_L(n, 0, matrix, lu->m);
#ifdef DEBUG
  _print_mat(n, matrix);
  cout << endl;
#endif
  for (int i = lu_block_size; i < n; i += lu_block_size) {
    block_gauss_0(n, i, matrix);
  }
  freeMat(lu);
}
