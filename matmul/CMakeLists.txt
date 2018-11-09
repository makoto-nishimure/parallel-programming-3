# CMakeのバージョンを設定
cmake_minimum_required(VERSION 3.12.4)
# プロジェクト名と使用する言語を設定
project(matmul CXX)


set(CMAKE_CXX_FLAGS "-std=c++11 -mavx -mfma")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O3 -D=TIME")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# NEED : -std=c++11
add_executable(stra
	stra_main.cpp # NEED : None, OPTIME : -O3, DEBUG : DEBUG, TIME
	mat_mul.cpp   # NEED : -mavx -mfma, OPTIME : -O3 : DEBUG, TIME
	mat_etc.cpp   # NEED : -mavx, OPTIME : -O3, DEBUG : DEBUG
	mm_opt.cpp		# NEED : None, OPTIME : -O3 DEBUG : None
	)
