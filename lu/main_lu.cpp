#include "lu.hpp"
int main (int argc, char **argv)
{
	int N = 1600;
	int n = 1500;

	adds *adx, *ady;
	adx = allocMat(N, N);
	ady = allocMat(N, N);

	mat matrix, lu;
	matrix = adx->m;
	lu = ady->m;

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i > j) {
				matrix[i][j] = ((i + j + 1) * 100 / 7);
			}
			else {
				matrix[i][j] = ((i + j + 1) * 100 / 3);
			}
		}
	}
	struct timeval tp[2];
	gettimeofday(tp, 0);
	//naive(n, matrix);
	block_gauss(n, matrix);
	gettimeofday(tp+1, 0);
	cout << double(elapsed_time(tp)) << endl;;
//#ifdef DEBUG
	//print_mat(n, matrix);
//#endif
	freeMat(adx);
	freeMat(ady);
	return 0;
}
