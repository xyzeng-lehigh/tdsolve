cmake_minimum_required(VERSION 3.9...3.31)
#set(CMAKE_C_COMPILER "gcc")
#set(CMAKE_CXX_COMPILER "g++")
project(TestBandSolve)
set(CMAKE_BUILD_TYPE Release)

# C++17 without GNU extensions (-std=c++17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# macOS with Apple clang: make sure that the C++ standard library is found
include(../../cmake/AppleClangStdlib.cmake)

# The headers of TDSolve
include_directories(../../include)

# Add source files: the solver to be tested and the linear algebra of
# TDSolve that it is compared with
file(GLOB SOURCES "TestBandSolve.cc"
     "../../src/Alg_BandSolve.cc" "../../src/Alg_Solve.cc" "../../src/Mat_Entry.cc"
     "../../src/HighDim_Alg.cc" "../../src/LowDim_Alg.cc" "../../src/LeastSquares_Solve.cc"
     "../../src/Quad_Coef.cc" "../../src/Comput_Choices.cc")

# Executable
add_executable(TestBandSolve ${SOURCES})

# LAPACK and BLAS
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(TestBandSolve ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# The Fortran runtime is needed by the reference BLAS/LAPACK on Linux, but
# not by the Accelerate framework of macOS (and clang has no libgfortran).
if(NOT APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lgfortran")
endif()
