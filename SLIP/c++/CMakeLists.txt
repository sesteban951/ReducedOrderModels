cmake_minimum_required(VERSION 3.10)

# Project name
project(examples)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)  # Need c++ 17 for inline variables
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O3 -pthread")

add_executable(ode src/example_ode.cpp)
