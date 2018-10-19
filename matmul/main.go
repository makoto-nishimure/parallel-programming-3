package main

import (
	"fmt"
	"time"
	"sync"
)

type Matrix [][]float64

func printMat(n int, mat Matrix) {
	for i := 0; i < n; i++ {
			fmt.Println(mat[i])
	}
}
func matMul_ijk_loop(n int, m1 Matrix, m2 Matrix) Matrix {
	mat := make(Matrix, n)
	for i := 0; i < n; i++ {
		mat[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < n; k++ {
				mat[i][j] += m1[i][k] * m2[j][k]
			}
		}
	}
	return mat
}

// ikj(?) loop
func matMul_ikj_loop(n int, m1 Matrix, m2 Matrix) Matrix {
	mat := make(Matrix, n)
	for i := 0; i < n; i++ {
		mat[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			tmp_v := m1[i][j]
			for k := 0; k < n; k++ {
				mat[i][k] += tmp_v * m2[j][k]
			}
		}
	}
	return mat
}

// ikj_loop and used thread
func matMul_thr(n int, m1 Matrix, m2 Matrix) Matrix {

	mat := make(Matrix, n)
	for i := 0; i < n; i++ {
		mat[i] = make([]float64, n)
	}
	thrNum := 2
	wg := new(sync.WaitGroup)
	//isFin := make(chan bool, len(thrNum))

	for thrID := 0; thrID < thrNum; thrID++ {
		wg.Add(1)
		go func (tID int) {
			defer wg.Done()
			start := tID * n / thrNum
			end := (tID + 1) * n / thrNum
//			fmt.Println("thiID :", tID, "s :", start, "e :", end)
			for i := start; i < end; i++ {
				for j := 0; j < n; j++ {
					tmp_v := m1[i][j]
					for k := 0; k < n; k++ {
						mat[i][k] += tmp_v * m2[j][k]
					}
				}
			}
		}(thrID)
	}
	wg.Wait()
//	close(isFin)
	return mat
}

// ikj_loop and used thread
func matMul_thr_unw(n int, m1 Matrix, m2 Matrix) Matrix {

	mat := make(Matrix, n)
	for i := 0; i < n; i++ {
		mat[i] = make([]float64, n)
	}
	thrNum := 2
	wg := new(sync.WaitGroup)
	//isFin := make(chan bool, len(thrNum))

	for thrID := 0; thrID < thrNum; thrID++ {
		wg.Add(1)
		go func (tID int) {
			defer wg.Done()
			start := tID * n / thrNum
			end := (tID + 1) * n / thrNum
//			fmt.Println("thiID :", tID, "s :", start, "e :", end)
			for i := start; i < end; i+=2 {
				in := i+1
				for j := 0; j < n; j++ {
					t1 := m1[i][j]
					t2 := m1[in][j]
					for k := 0; k < n; k++ {
						mat[i][k] += t1 * m2[j][k]
						mat[in][k] += t2 * m2[j][k]
					}
				}
			}
		}(thrID)
	}
	wg.Wait()
//	close(isFin)
	return mat
}
// ikj_loop, block
func matMul_only_block (n int, m1 Matrix, m2 Matrix) Matrix {
	BLOCK := 50
	mat := make(Matrix, n)
	for i := 0; i < n; i++ {
		mat[i] = make([]float64, n)
	}
	for i := 0; i < n; i+=BLOCK {
		for j := 0; j < n; j+=BLOCK {
			for k := 0; k < n; k+=BLOCK {
				ib := i + BLOCK
				jb := j + BLOCK
				kb := k + BLOCK
				for ii := i; ii < ib; ii++ {
					for jj := j; jj < jb; jj++ {
						tmp_v := m1[ii][jj]
						for kk := k; kk < kb; kk++ {
							mat[ii][kk] += tmp_v * m2[jj][kk]
						}
					}
				}
			}
		}
	}
	return mat
}

// ikj_loop, block and used thread
func matMul_block (n int, m1 Matrix, m2 Matrix) Matrix {
	BLOCK := 10
	mat := make(Matrix, n)
	for i := 0; i < n; i++ {
		mat[i] = make([]float64, n)
	}
	thrNum := 2
	wg := new(sync.WaitGroup)

	for thrID := 0; thrID < thrNum; thrID++ {
		wg.Add(1)
		go func (tID int) {
			defer wg.Done()
			start := tID * n / thrNum
			end := (tID + 1) * n / thrNum
			//			fmt.Println("thiID :", tID, "s :", start, "e :", end)
			for i := start; i < end; i+=BLOCK {
				for j := 0; j < n; j+=BLOCK {
					for k := 0; k < n; k+=BLOCK {
						ib := i + BLOCK
						jb := j + BLOCK
						kb := k + BLOCK
						for ii := i; ii < ib; ii++ {
							for jj := j; jj < jb; jj++ {
								tmp_v := m1[ii][jj]
								for kk := k; kk < kb; kk++ {
									mat[ii][kk] += tmp_v * m2[jj][kk]
								}
							}
						}
					}
				}
			}
		}(thrID)
	}
	wg.Wait()
	return mat
}

// ikj_loop, block, loop unwinding and used thread
func matMul_unw (n int, m1 Matrix, m2 Matrix) Matrix {
	BLOCK := 100
	mat := make(Matrix, n)
	for i := 0; i < n; i++ {
		mat[i] = make([]float64, n)
	}
	thrNum := 2
	wg := new(sync.WaitGroup)

	for thrID := 0; thrID < thrNum; thrID++ {
		wg.Add(1)
		go func (tID int) {
			defer wg.Done()
			start := tID * n / thrNum
			end := (tID + 1) * n / thrNum
			//			fmt.Println("thiID :", tID, "s :", start, "e :", end)
			for i := start; i < end; i+=BLOCK {
				for j := 0; j < n; j+=BLOCK {
					for k := 0; k < n; k+=BLOCK {
						ib := i + BLOCK
						jb := j + BLOCK
						kb := k + BLOCK
						for ii := i; ii < ib; ii+=2 {
							in := ii + 1
							for jj := j; jj < jb; jj++ {
								t1 := m1[ii][jj]
								t2 := m1[in][jj]
								for kk := k; kk < kb; kk++ {
									mat[ii][kk] += t1 * m2[jj][kk]
									mat[in][kk] += t2 * m2[jj][kk]
								}
							}
						}
					}
				}
			}
		}(thrID)
	}
	wg.Wait()
	return mat
}

func main() {
	var n int
	fmt.Scan(&n)

	m1 := make(Matrix, n)
	m2 := make(Matrix, n)
	for i := 0; i < n; i++ {
		m1[i] = make([]float64, n)
		m2[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			//			m1[i][j] = (float64)(i + j)
			//			m2[i][j] = (float64)(i + j)
			m1[i][j] = float64((i + j) * 100 / 7.0)
			m2[i][j] = float64((i + j) * 300 / 7.0)
		}
	}
	//printMat(n, m1)
	//printMat(n, m2)

	st := time.Now()
	var _ = matMul_ikj_loop(n, m1, m2)
	et := time.Now()
	fmt.Println("thr_ikj_loop")
	fmt.Println(et.Sub(st))

	st = time.Now()
	var _ = matMul_only_block(n, m1, m2)
	et = time.Now()
	fmt.Println("thr_only_block")
	fmt.Println(et.Sub(st))

	st = time.Now()
	var _ = matMul_thr(n, m1, m2)
	et = time.Now()
	fmt.Println("thr")
	fmt.Println(et.Sub(st))

	st = time.Now()
	var _ = matMul_thr_unw(n, m1, m2)
	et = time.Now()
	fmt.Println("thr_unw")
	fmt.Println(et.Sub(st))

	st = time.Now()
	var _ = matMul_unw(n, m1, m2)
	et = time.Now()
	fmt.Println("thr_blo_unw")
	fmt.Println(et.Sub(st))
	//printMat(n, m3)
}
