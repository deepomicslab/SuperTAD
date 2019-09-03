cmake_minimum_required(VERSION 3.14)
project(program)

set(CMAKE_CXX_STANDARD 11)

add_executable(program main.cpp data.cpp data.h params.cpp params.h detectorBinary.cpp detectorBinary.h binaryTree.cpp binaryTree.h inputAndOutput.cpp inputAndOutput.h utils.cpp utils.h detectorMulti.cpp detectorMulti.h)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(lib)