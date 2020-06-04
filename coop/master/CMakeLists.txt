# CMakeのバージョンを設定
cmake_minimum_required(VERSION 3.12.4)
# プロジェクト名と使用する言語を設定
project(matmul CXX)

set(CMAKE_CXX_FLAGS "-std=c++11 -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
#cmake -DCMAKE_BUILD_TYPE=Debug
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#cmake -DCMAKE_BUILD_TYPE=Release

# NEED : -std=c++11
add_executable(master
    master.cpp
    ../opt.cpp
	)
