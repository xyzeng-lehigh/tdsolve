cmake_minimum_required(VERSION 3.9)
#set(CMAKE_C_COMPILER "gcc")
#set(CMAKE_CXX_COMPILER "g++")
project(TDSolve)
set(CMAKE_BUILD_TYPE Release)

# The major version number: Number of problems supported
set(TDS_MAJOR_VERSION 2)
# The minor version number: Order of accuracy
set(TDS_MINOR_VERSION 2)
# The growth version number:
set(TDS_GRWTH_VERSION 1)
# The full version
set(TDS_VERSION
    ${TDS_MAJOR_VERSION}.${TDS_MINOR_VERSION}.${TDS_GRWTH_VERSION}
    )

# Bring the headers into the project
include_directories(include)

# Add source files
file(GLOB SOURCES "main/TumorDyn.cc" "src/*.cc")
file(GLOB PARSOURCES "main/ParAnal.cc" "src/*.cc")

# Executable
add_executable(TumorDynRun ${SOURCES})
add_executable(ParAnalRun ${PARSOURCES})

# LAPACK and BLAS
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
include_directories(${BLAS_INCLUDE_DIR} ${LAPACK_INCLUDE_DIR})
target_link_libraries(TumorDynRun ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Automatic differentiation
# By default the built-in implementation is used (include/AutoDiff*.h and
# src/AutoDiff.cc) and Trilinos is not needed. To use the Sacado package of
# Trilinos instead, e.g. for cross-checking, configure with
#   -DUSE_TRILINOS=ON -DTrilinos_DIR:PATH=<trilinos install directory>
option(USE_TRILINOS "Use Sacado of Trilinos for automatic differentiation" OFF)
if(USE_TRILINOS)
  SET(CMAKE_PREFIX_PATH ${Trilinos_DIR} ${CMAKE_PREFIX_PATH})
  SET(CMAKE_PREFIX_PATH ${Trilinos_PATH} ${CMAKE_PREFIX_PATH})
  find_package(Trilinos REQUIRED)
  # Trilinos 12.7 has many duplicates and causes problem in ld on some Linux systems
  list(REMOVE_DUPLICATES Trilinos_LIBRARIES)
  list(REMOVE_DUPLICATES Trilinos_TPL_INCLUDE_DIRS)
  list(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
  include_directories(${Trilinos_INCLUDE_DIRS})
  include_directories(${Trilinos_TPL_INCLUDE_DIRS})
  MESSAGE("\nFound Trilinos!  Here are the details: ")
  MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
  MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
  MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
  MESSAGE("End of Trilinos details\n")
  set(CMAKE_CXX_FLAGS "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
  string (REPLACE ";" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  add_definitions(-DTDS_USE_TRILINOS)

  target_link_libraries(TumorDynRun ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -lgfortran")

# OpenMP support (avaiable after CMake 3.9)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(TumorDynRun OpenMP::OpenMP_CXX)
  target_link_libraries(ParAnalRun OpenMP::OpenMP_CXX)
  #target_link_libraries(ParAnalRun PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(ParAnalRun ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
if(USE_TRILINOS)
  target_link_libraries(ParAnalRun ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
endif()

