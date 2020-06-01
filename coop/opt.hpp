#include <iostream>
#include <netdb.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>

int tcp_acc_port(int port);
int fdopen_sock(int sock, FILE **in, FILE **out);

void expected_replay(FILE *is, FILE *os, const char *str);
void print_escaped(const char *p, FILE *fp);
