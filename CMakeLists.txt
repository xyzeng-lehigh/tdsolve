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
#file(GLOB PARSOURCES "main/ParAnal.cc" "src/*.cc")

# Executable
add_executable(TumorDynRun ${SOURCES})
#add_executable(ParAnalRun ${PARSOURCES})

# LAPACK and BLAS
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
include_directories(${BLAS_INCLUDE_DIR} ${LAPACK_INCLUDE_DIR})
target_link_libraries(TumorDynRun ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Trilinos
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

target_link_libraries(TumorDynRun ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -lgfortran")

# Configure file that pass CMake settings to the source
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
#include(FindOpenMP)
#if(OPENMP_FOUND)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#else(OPENMP_FOUND)
#  message("ERROR: OpenMP could not be found.")
#endif(OPENMP_FOUND)

## OpenMP support (avaiable after CMake 3.9)
#find_package(OpenMP)
#if(OpenMP_CXX_FOUND)
#  target_link_libraries(ParAnalRun OpenMP::OpenMP_CXX)
#  #target_link_libraries(ParAnalRun PUBLIC OpenMP::OpenMP_CXX)
#endif()
#
#target_link_libraries(ParAnalRun ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
#target_link_libraries(ParAnalRun ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

