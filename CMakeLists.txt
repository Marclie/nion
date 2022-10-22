cmake_minimum_required(VERSION 3.22)
project(nion)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_executable(nion main.cpp nion.hpp)
