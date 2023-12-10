cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

project(NumericalMethods CXX)

add_executable(NumericalMethods
    #Library executable
    CPP/include/gauss/gauss.cpp
    #CPP/include/neutone/neuton.cpp

    #main
    CPP/src/main.cpp
)

target_include_directories(NumericalMethods PUBLIC CPP/include/)