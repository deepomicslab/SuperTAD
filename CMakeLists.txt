cmake_minimum_required(VERSION 3.5.1)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER gcc)
set(CMAKE_C_FLAGS -O3)

project(program)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable(SuperTAD main.cpp data.cpp data.h params.cpp params.h detectorBinary.cpp detectorBinary.h
        binaryTree.cpp binaryTree.h inputAndOutput.cpp inputAndOutput.h utils.cpp utils.h
        detectorMulti.cpp detectorMulti.h multiTree.cpp multiTree.h detectorH.cpp detectorH.h)
