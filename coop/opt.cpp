#include "opt.hpp"
#include <math.h>

#define BUFLEN 1048576
char msg_buf[BUFLEN];

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

void genvec(int n, double *x, int k) {
  int i;
  double dd = 2 * 3.141592653589793238462643383 / n * k;
  for (i = 0; i < n; i++) {
    x[i] = ((i % 19) + 1) * (sin(i * dd) + 1.1);
  }
}

int tcp_acc_port(int port) {
  struct sockaddr_in addr;
  int addr_len;
  int s;

  if ((s = socket(PF_INET, SOCK_STREAM, 0)) < 0) {
    perror("socket");
    return (-1);
  }

  memset(&addr, 0, sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_addr.s_addr = INADDR_ANY;
  addr.sin_port = htons(port);

  if (bind(s, (struct sockaddr *)&addr, sizeof(addr)) < 0) {
    perror("bine");
    fprintf(stderr, "port number %d is already used. wait a moment.\n", port);
    return (-1);
  }
  if (listen(s, 5) < 0) {
    perror("listen");
    close(s);
    return (-1);
  }
  return (s);
}

int fdopen_sock(int sock, FILE **in, FILE **out) {
  int sock2;
  if ((sock2 = dup(sock)) < 0) {
    return (-1);
  }
  if ((*in = fdopen(sock2, "r")) == NULL) {
    close(sock2);
    return (-1);
  }
  if ((*out = fdopen(sock, "w")) == NULL) {
    fclose(*in);
    *in = 0;
    return (-1);
  }
  setvbuf(*out, (char *)NULL, _IOLBF, 0);
  return (0);
}

/*
 期待される reply が返って来たかをチェックする関数。
 期待していない reply が帰って来たら、即座に終了。
 error 時にも、即座に終了。

  FILE *is : input stream 、クライアントからの応答がここに届く
  FILE *os : output stream、ここへ書き出すとクライアントに届く
  char *str : 期待される reply の文字列
 */
void expected_replay(FILE *is, FILE *os, const char *str) {
  fscanf(is, "%s", msg_buf); /* reply を読み込む */

  /* エラーだったら、エラーだと出力して終わる。 */
  if (!strcmp(msg_buf, "error")) {
    fprintf(stderr, "got an error. abort\n");
    fflush(stderr);
    exit(0);
  }

  if (!strcmp(msg_buf, "timeout")) {
    /* report コマンドでレポートをもらう */
    fprintf(os, "report\n");
    fflush(os);

    /* report reply が来る */
    expected_replay(is, os, "report");

    /* レポート本体を読み込む */
    fscanf(is, "%s", msg_buf); // the report body

    /* レポートをエスケープ解除しつつ表示 */
    fprintf(stderr, "reported: \n");
    print_escaped(msg_buf, stderr);
    fflush(stderr);

    /* クライアントを停止させる terminate コマンドを最後に送る */
    fprintf(os, "terminate\n");
    fflush(os);

  }

  /* 期待した返事かどうかをチェック */
  else if (strcmp(msg_buf, str)) {
    fprintf(stderr, "unexpected reply from the client: %s ('%s' expected)\n",
            msg_buf, str);
    fprintf(stderr, "abort.\n");
    fflush(stderr);

    exit(0);
  }
}

/* エスケープを解除しつつ表示する */
void print_escaped(const char *p, FILE *fp) {
  int c = 0;
  while (*p != 0) {
    if (*p == '%') {
      p++;
      if (*p == 0)
        return;
      else if (*p == '%')
        c = '%';
      else if (*p == 's')
        c = ' ';
      else if (*p == 'n')
        c = '\n';
      else {
        fprintf(stderr, "invalid escape: %s", p - 1);
        exit(0);
      }
    } else {
      c = *p;
    }
    fputc(c, fp);
    p++;
  }
}
