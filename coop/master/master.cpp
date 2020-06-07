//司令塔となるプログラム
#include <netdb.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

#include "master.hpp"

#define BUFFERSIZE 1024
#define BUFLEN 1048576

using namespace std;

char host[BUFLEN];

FILE *is, *os; // from/to server to/from client

mutex client;

bool START = false;

void create_thr(int port, int worker_num);
void each_thr(int com);

int main(int argc, char **argv) {
  int port, sub_port, worker_num;

  if (argc != 5) {
    fprintf(
        stdout,
        "Usage: %s [client_add] [client_port] [port_for_calc] [worker_num]\n",
        argv[0]);
    exit(0);
  } else {
    strncpy(host, argv[1], BUFLEN - 1); // hostname
    sscanf(argv[2], "%d", &port);       // port
    sub_port = strtol(argv[3], 0, 10);
    sscanf(argv[4], "%d", &worker_num); // worker_num
  }

  { /* 通信準備 = os と is を準備する */
    struct sockaddr_in addr;
    struct hostent *serv;
    int sd = socket(AF_INET, SOCK_STREAM, 0); // create a socket

    serv = gethostbyname(host); // resolve the hostname
    if (serv == NULL) {
      fprintf(stderr, "unknown host %s", host);
      exit(0);
    }

    addr.sin_family = AF_INET;
    addr.sin_port = htons(port);
    bcopy(serv->h_addr, (char *)&addr.sin_addr, serv->h_length);
    fprintf(stderr, "connecting to %s::%d\n", host, port);
    // connect to the client at host::port
    if (connect(sd, (struct sockaddr *)&addr, sizeof(addr)) == -1) {
      fprintf(stderr, "failed in connecting to %s::%d\n", host, port);
      fprintf(stderr, "abort.");
      exit(0);
    }
    fprintf(stderr, "connected.\n");
    os = fdopen(sd, "w"); // output stream を作る
    is = fdopen(sd, "r"); // input stream を作る
  }
  create_thr(sub_port, worker_num);
  return 0;
}

void create_thr(int port, int worker_num) {
  int acc, com;
  acc = tcp_acc_port(port);
  if (acc < 0) {
    perror("tcp_acc_port error");
    exit(-1);
  }

  for (int i = 0; i < worker_num; ++i) {
    com = accept(acc, 0, 0);
    if (com < 0) {
      perror("accept error");
      exit(-1);
    }
    //
    if (i != worker_num - 1) {
      thread t(each_thr, com);
      t.detach();
    }
    // last worker
    else {
      fprintf(os, "ready\n");
      fflush(os);
      expected_replay(is, os, "start");
      START = true;
      each_thr(com);
    }
  }
}

void getProblem(int &pi, int &n, int &k, int *ps) {
  lock_guard<mutex> lock(client);
  fprintf(os, "next\n");
  fflush(os);
  expected_replay(is, os, "problem");

  // read problem specification from stdio
  fscanf(is, "%d", &pi); // problem id
  fscanf(is, "%d", &n);  // size of matrices
  fscanf(is, "%d", &k);  // # of matrices
  for (int i = 0; i < k; ++i)
    fscanf(is, "%d", &ps[i]); // parameters for matrices
}

void sendProblem(int &pi, int &n, int &k, int *ps, FILE *out) {
  fprintf(out, "%d ", pi);
  fprintf(out, "%d ", n);
  fprintf(out, "%d ", k);
  for (int i = 0; i < k; ++i)
    fprintf(out, "%d ", ps[i]);
  fflush(out);
}

void getVecSeed(int &pi, int &l, int *vecSeed) {
  lock_guard<mutex> lock(client);
  //ベクトルの種受信
  fprintf(os, "lu_done %d\n", pi);
  fflush(os);
  expected_replay(is, os, "instances");

  fscanf(is, "%d", &pi);
  fscanf(is, "%d", &l);
  for (int i = 0; i < l; ++i) {
    fscanf(is, "%d", &vecSeed[i]);
  }
  //ベクトルの種受信完了
}

void sendVecSeed(int &l, int *vecSeed, FILE *out) {
  //他計算機にベクトルの種送信
  fprintf(out, "%d ", l);
  for (int i = 0; i < l; ++i) {
    fprintf(out, "%d ", vecSeed[i]);
  }
  fflush(out);
  //他計算機にベクトルの種送信完了
}

void getAns(int &n, int &l, double **x, FILE *in) {
  int i, j;
  for (i = 0; i < l; ++i) {
    //他計算機から解受信
    for (j = 0; j <= n - 2; j += 2) {
      fscanf(in, "%lf", &x[i][j]);
      fscanf(in, "%lf", &x[i][j + 1]);
    }
    for (; j < n; ++j)
      fscanf(in, "%lf", &x[i][j]);
    //他計算機から解受信完了
  }
}

void sendAns(int &pi, int &n, int &l, double **x) {
  int i, j;
  int s;
  lock_guard<mutex> lock(client);
  fprintf(os, "solved\n");
  fprintf(os, "%d\n%d\n", pi, l);
  fflush(os);

  for (i = 0; i < l; ++i) {
    //解送信
    fprintf(os, "%d\n", n);
    fflush(os);

    for (j = 0; j <= n - 2; j += 2) {
      fprintf(os, "%.15e\n%.15e\n", x[i][j], x[i][j + 1]);
      fflush(os);
    }
    for (; j < n; ++j) {
      fprintf(os, "%.15e\n", x[i][j]);
      fflush(os);
    }
  }

  expected_replay(is, os, "score");
  fscanf(is, "%d", &s);
}

void each_thr(int com) {
  FILE *in, *out; // from/to server to/from workers

  int pi, n, k, ps[32], vecSeed[32];
  int npi, nn;
  int l, s;
  int xSIZE = 8;

  int compLU = 0;

  if (fdopen_sock(com, &in, &out) < 0) {
    fprintf(stderr, "fdopen()\n");
    exit(1);
  }

  double **x = (double **)malloc(sizeof(double) * xSIZE); //ベクトルx
  for (int i = 0; i < xSIZE; ++i) {
    x[i] = (double *)malloc(sizeof(double) * 1600);
  }

  while (true) {
    if (START)
      break;
  }

  while (true) {
    getProblem(pi, n, k, ps);
    sendProblem(pi, n, k, ps, out);
    // LU分解終了フラグ
    fscanf(in, "%d", &compLU);

    if (compLU == 1) {
      getVecSeed(pi, l, vecSeed);
      sendVecSeed(l, vecSeed, out);

      getAns(n, l, x, in);
      sendAns(pi, n, l, x);
    }
    compLU = 0;
  }

  for (int i = 0; i < xSIZE; i++) {
    free(x[i]);
  }
  free(x);

  fclose(in);
  fclose(out);
}
