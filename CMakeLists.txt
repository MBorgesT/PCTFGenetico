cmake_minimum_required(VERSION 3.8)
project(PCTFGenetico)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc -fopenmp")

add_executable(PCTFGenetico pctfga.cpp)