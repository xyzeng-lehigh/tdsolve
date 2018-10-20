cmake_minimum_required(VERSION 2.8)
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

# Configure file that pass CMake settings to the source
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# Executable
add_executable(TumorDynRun ${SOURCES})
add_executable(ParAnalRun ${PARSOURCES})
