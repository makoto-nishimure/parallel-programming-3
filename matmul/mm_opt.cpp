#include "mm_opt.hpp"

double elapsed_time(struct timeval tp[2])
{
	return tp[1].tv_sec-tp[0].tv_sec+1e-6*(tp[1].tv_usec-tp[0].tv_usec);
}

void print_pair(int size, point pa)
{
	cout << pa.first << "," << pa.second << " to " <<
		pa.first + size - 1 << "," << pa.second + size - 1 << endl;
}

void print_mat(int size, mat m)
{
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			printf("%.5lf ", m[i][j]);
		}
		cout << endl;
	}
}
