# CMakeのバージョンを設定
cmake_minimum_required(VERSION 3.12.4)
# プロジェクト名と使用する言語を設定
project(matmul CXX)


set(CMAKE_CXX_FLAGS "-std=c++11 -mavx -mfma -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "-g -D=DEBUG")
#cmake -DCMAKE_BUILD_TYPE=Debug
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#cmake -DCMAKE_BUILD_TYPE=Release

# NEED : -std=c++11
add_executable(lu
	main_lu.cpp # NEED : None, OPTIME : -O3, DEBUG : DEBUG, TIME
	lu.cpp # NEED : None, OPTIME : -O3, DEBUG : DEBUG, TIME
	../matmul/mat_mul.cpp   # NEED : -mavx -mfma, OPTIME : -O3 : DEBUG, TIME
	../matmul/mat_etc.cpp   # NEED : -mavx, OPTIME : -O3, DEBUG : DEBUG
	../matmul/mm_opt.cpp		# NEED : None, OPTIME : -O3 DEBUG : None
	)
