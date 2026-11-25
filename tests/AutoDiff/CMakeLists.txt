cmake_minimum_required(VERSION 3.9)
#set(CMAKE_C_COMPILER "gcc")
#set(CMAKE_CXX_COMPILER "g++")
project(TestAutoDiff)
set(CMAKE_BUILD_TYPE Release)

# The automatic differentiation headers of TDSolve
include_directories(../../include)

# Add source files
file(GLOB SOURCES "TestAutoDiff.cc" "../../src/AutoDiff.cc")

# Executable
add_executable(TestAutoDiff ${SOURCES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
