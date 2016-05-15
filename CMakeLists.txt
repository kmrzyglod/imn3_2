cmake_minimum_required(VERSION 3.5)
project(imn3_2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp FlagMatrix.h imnmath.h Point.h imnmath.h Point.h Relaxation.h RelaxationNeumann.h IntegralGraphPoint.h RelaxationDirichlet.h)
add_executable(imn3_2 ${SOURCE_FILES})