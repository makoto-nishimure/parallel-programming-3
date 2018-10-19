// Strassen algorithm
package main

import (
	"fmt"
	"time"
	"sync"
)

type Matrix [][]float64

func matAdd(n int, i1 int, j1 int, i2 int, j2 int, m1 Matrix) Matrix {
	mat := make(Matrix, n)
	for i := 0; i < n; i++ {
		mat[i] = make([]float64, n)
	}
	for i := 0; i < n; i++, i1++, i2++ {
		for j := 0; j < n; j++, j1++, j2++ {
			mat[i][j] = m1[i1][j1] + m2[i2][j2]
		}
	}
	return mat
}

func matSub(n int, i1 int, j1 int, i2 int, j2 int, m1 Matrix) Matrix {
	mat := make(Matrix, n)
	for i := 0; i < n; i++ {
		mat[i] = make([]float64, n)
	}
	for i := 0; i < n; i++, i1++, i2++ {
		for j := 0; j < n; j++, j1++, j2++ {
			mat[i][j] = m1[i1][j1] - m2[i2][j2]
		}
	}
	return mat
}

// ikj(?) loop
func matMul_ikj_loop(n int, i1 int, j1 int, i2 int, j2 int, m1 Matrix, m2 Matrix) Matrix
{
	mat := make(Matrix, n)
	for i := 0; i < n; i++ {
		mat[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			tmp_v := m1[i1 + i][j1 + j]
			for k := 0; k < n; k++ {
				mat[i][k] += tmp_v * m2[i2 + j][j2 + k]
			}
		}
	}
	return mat
}

func sa(n int, mat_a Matrix, mat_b Matrix) Matrix {
	half_n := n/2
	// a11 + a22
	tmp1 := matAdd(half_n, 0, 0, half_n, half_n, mat_a)
	// b11 + b22
	tmp2 := matAdd(half_n, 0, 0, half_n, half_n, mat_a)
	m1 := matMul(half_n, 0, 0, 0, 0, mat_a, mat_b)
	m5 := matMul(half_n, 0, 0, half_n, half_n, tmp1, mat_b)
	m2 := matMul(half_n, 0, 0, 0, 0, matAdd(half_n, half_n, 0, half_n, half_n, mat_a, mat_a), mat_b)
	m3 := matMul(half_n, 0, 0, 0, 0, mat_a, mat)

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
	//var _ = matMul_thr(n, m1, m2)
	et := time.Now()
	fmt.Println("thr")
	fmt.Println(et.Sub(st))
}
