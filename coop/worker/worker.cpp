#include "worker.hpp"
#include "../../matmul/mat_etc.hpp"
#include "../../matmul/mat_mul.hpp"
#include "../../matmul/mm.hpp"

int echo_client(char *serv, int port);
int tcp_connect(char *serv, int port);

/* 行列生成関数 */
void genmat2(int n, double *a, int k) {
  int i, j;
  double d, aij;
  int nk = n / k;
  double nk1 = 1. / nk;
  for (j = 0; j < n * 2; j++) {
    d = ((n - abs(n - j)) / (double)n) / (abs(n - j) + 1);
    aij = d * ((nk / 2 - (abs(n - j) % nk)) * nk1 + 1);
    a[j] = (0.1 * aij);
  }
  a[n] += 1;
}

void matrix_init2(double *a, int n) {
  int i;
  for (i = 0; i < n * 2; i++) {
    a[i] = 0;
  }
}

// n:size z:=x*y  x:掛けられる行列 y:掛ける行列
void mat_mul2(const size_t n, const size_t nex, const double end,
              const double num, double **z, double **x, double *y) {
  // static const size_t double_size = 4;

  //行列xのx[i][j]とx[i+1][j]要素をそれぞれ格納する  [i][j]と[i+1][j]など...
  __m256d a, b;

  //行列x[i]とx[i+1]へのポインタ
  double *ax, *bx;

  //行列z[i]とz[i+1]へのポインタ
  __m256d *ans1, *ans2;
  //行列y[i]へのポインタ
  __m256d *s, *t;

  for (size_t i = 0; i < nex; i++) {

    //行列x[i]とx[i+1]へのポインタ指定
    //行列z[i]とz[i+1]へのポインタ指定
    ax = x[i];
    ans1 = (__m256d *)z[i];
    bx = x[++i];
    ans2 = (__m256d *)z[i];

    for (size_t j = 0; j < n; j++) {

      size_t tmp = n - j;

      //行列y[i]へのポインタ   kループが上手く回らないのを防ぐためにsとtの二つ
      s = (__m256d *)&y[tmp];
      t = (__m256d *)&y[tmp];

      //行列x[i][j]の要素をaに4つコピー
      a = _mm256_broadcast_sd(ax + j);
      //行列x[i+1][j]の要素をaに4つコピー
      b = _mm256_broadcast_sd(bx + j);

      for (size_t k = 0; k <= end; k++) {
        //行列計算
        ans1[k] += _mm256_mul_pd(a, s[k]);
        ans2[k] += _mm256_mul_pd(b, t[k]);
      }
    }
  }
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j <= num; j++)
      z[n - i - 1][n - j - 1] = z[i][j];
}

void naive(int n, mat a) {
  int i, j, k;
  for (i = 0; i < n; i++) {
    for (k = 0; k < i; k++)
      for (j = k + 1; j < n; j++)
        a[i][j] -= a[i][k] * a[k][j];
    for (j = i + 1; j < n; j++)
      a[i][j] /= a[i][i];
    //	print_mat(n, a);
    //	cout << endl;
  }
}

void solveEqL(mat L, double *b, double *y, int size) {
  int i, j;
  double tmp;

  for (i = 0; i < size; i++) {
    tmp = b[i];
    for (j = 0; j < i; j++) {
      tmp -= L[i][j] * y[j];
    }
    y[i] = tmp / L[i][i];
  }
}

void solveEqU(mat U, double *y, double *x, int size) {
  int i, j;
  double tmp;

  for (i = size - 1; i >= 0; i--) {
    tmp = y[i];
    for (j = size - 1; j > i; j--) {
      tmp -= U[i][j] * x[j];
    }
    x[i] = tmp;
  }
}

void solveEq(mat m, double *b, double *x, int size) {
  double *y = (double *)malloc(sizeof(double) * size);
  solveEqL(m, b, y, size);
  solveEqU(m, y, x, size);
}

int main(int argc, char **argv) {
  char *serv;
  int port;
  if (argc != 3) {
    fprintf(stdout, "Usage: %s host port\n", argv[0]);
    exit(-1);
  }
  serv = argv[1];
  port = strtol(argv[2], 0, 10);
  echo_client(serv, port);

  return 0;
}

int echo_client(char *serv, int port) {
  int sock;
  FILE *in, *out;

  int i, j;
  int pi, n, k, ps[32], l, vecSeed[32];
  int nex;

  int max = 1600;
  int size = 1504;

  double end, num;

  sock = tcp_connect(serv, port);
  if (sock < 0)
    exit(-1);
  if (fdopen_sock(sock, &in, &out) < 0) {
    fprintf(stderr, "fdopen()\n");
    exit(1);
  }

  adds *ad, *adx, *ady;
  mat matrix, m1, m2;

  ad = allocMat(max, max);
  adx = allocMat(max, max);
  ady = allocMat(max, max);

  matrix = ad->m;
  m1 = adx->m;
  m2 = ady->m;

  double *b = (double *)malloc(sizeof(double) * max); //ベクトルb
  double *x = (double *)malloc(sizeof(double) * max); //ベクトルx

  for (;;) {
    set_zero(max, matrix);
    set_zero(max, m1);
    set_zero(max, m2);

    //接続先から問題受信
    fscanf(in, "%d", &pi);
    fscanf(in, "%d", &n);
    fscanf(in, "%d", &k);
    for (i = 0; i < k; i++) {
      fscanf(in, "%d", &ps[i]);
    }

    // 行列乗算 開始
    if (1 == k) {
      genmat(n, matrix, ps[0]);
    } else {
      genmat(n, m1, ps[0]);
      for (i = 1; i < k; ++i) {
        genmat(n, m2, ps[i]);
        str(0, size, m1, m2, matrix);
        matCopy(size, make_pair(0, 0), matrix, m1);
      }
    }

    // LU分解 開始
    naive(size, matrix);
    // LU分解 終了

    // LU完了通知
    i = 1;
    fprintf(out, "%d ", i);
    fflush(out);

    //　親プログラム　から l と　ベクトルの種受信
    fscanf(in, "%d", &l);
    for (i = 0; i < l; i++) {
      fscanf(in, "%d", &vecSeed[i]);
    }

    int tmp = n - 2;

    for (i = 0; i < l; i++) {
      // ベクトル生成し、解を求める
      genvec(n, b, vecSeed[i]);
      solveEq(matrix, b, x, n);

      for (j = 0; j <= tmp; j += 2) {
        fprintf(out, "%.15e\n%.15e\n", x[j], x[j + 1]);
        fflush(out);
      }
      for (; j < n; j++) {
        fprintf(out, "%.15e\n", x[j]);
        fflush(out);
      }
    }
  }
  fclose(in);
  fclose(out);
}

int tcp_connect(char *serv, int port) {
  struct addrinfo hints, *ai;
  struct sockaddr_in addr;
  int s;
  int err;

  if ((s = socket(PF_INET, SOCK_STREAM, 0)) < 0) {
    perror("socket");
    return (-1);
  }

  memset(&hints, 0, sizeof(hints));
  hints.ai_family = AF_INET;
  hints.ai_socktype = SOCK_STREAM;
  if ((err = getaddrinfo(serv, NULL, &hints, &ai))) {
    fprintf(stderr, "unknown host %s (%s)\n", serv, gai_strerror(err));
    close(s);
    return (-1);
  }
  if (ai->ai_addrlen > sizeof(addr)) {
    fprintf(stderr, "sockaddr too lage (%d) > (%d)\n", ai->ai_addrlen,
            sizeof(addr));
    freeaddrinfo(ai);
    close(s);
    return (-1);
  }
  memcpy(&addr, ai->ai_addr, ai->ai_addrlen);
  addr.sin_port = htons(port);
  freeaddrinfo(ai);

  if (connect(s, (struct sockaddr *)&addr, sizeof(addr)) < 0) {
    perror(serv);
    close(s);
    return (-1);
  }
  return s;
}
