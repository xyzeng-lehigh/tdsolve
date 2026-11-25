cmake_minimum_required(VERSION 3.9...3.31)
#set(CMAKE_C_COMPILER "gcc")
#set(CMAKE_CXX_COMPILER "g++")
project(TestAutoDiff)
set(CMAKE_BUILD_TYPE Release)

# C++17 without GNU extensions (-std=c++17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# macOS with Apple clang: make sure that the C++ standard library is found
include(../../cmake/AppleClangStdlib.cmake)

# The automatic differentiation headers of TDSolve
include_directories(../../include)

# Add source files
file(GLOB SOURCES "TestAutoDiff.cc" "../../src/AutoDiff.cc")

# Executable
add_executable(TestAutoDiff ${SOURCES})
