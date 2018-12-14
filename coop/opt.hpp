#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <string.h>

int tcp_acc_port(int port);
int fdopen_sock(int sock, FILE **in, FILE **out);

void expected_replay(FILE *is, FILE *os, const char *str);
void print_escaped(const char *p, FILE*fp);
