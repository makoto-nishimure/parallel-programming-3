//#include "mm.hpp"
#include "mat_mul.hpp"
#include "mm_opt.hpp"


int main (int argc, char **argv)
{
	mat m1, m2, m3;
	int n = 1504;
	int N = 1504;
	m1 = (double **)_mm_malloc(sizeof(double*)*N, 64);
	m2 = (double **)_mm_malloc(sizeof(double*)*N, 64);
	m3 = (double **)_mm_malloc(sizeof(double*)*N, 64);
	for (int i = 0; i < N; ++i) {
		m1[i] = (double *)_mm_malloc(sizeof(double)*N, 64);
		m2[i] = (double *)_mm_malloc(sizeof(double)*N, 64);
		m3[i] = (double *)_mm_malloc(sizeof(double)*N, 64);
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			m1[i][j] = ((i + j) * 100 / 7.0);
			m2[i][j] = ((i + j) * 300 / 7.0);
		}
	}
	struct timeval tp[2];
	int roop = 10;
	gettimeofday(tp, 0);
	for (int i = 0; i < roop; ++i) {
		set_zero(n, m3);
		str(0, n, m1, m2, m3);
	}
	gettimeofday(tp+1, 0);
	//cout << "STR DONE" << endl;
	double one_exe_time = double(elapsed_time(tp)) / roop;
	cout << one_exe_time <<endl;;
	//print_mat(n, m3);
#ifdef TIME
	cout << "malloc and free TIME" <<endl;
	cout << mtime / roop <<endl;
	cout << (mtime / roop) / one_exe_time <<endl;;
	cout << "set zero TIME" <<endl;
	cout << set_z / roop <<endl;
	cout << (set_z / roop) / one_exe_time <<endl;
#endif

	for (int i = 0; i < N; ++i) {
		_mm_free(m1[i]);
		_mm_free(m2[i]);
		_mm_free(m3[i]);
	}
	_mm_free(m1);
	_mm_free(m2);
	_mm_free(m3);
	return 0;
}
