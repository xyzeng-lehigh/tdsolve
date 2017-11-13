cmake_minimum_required(VERSION 3.5)
project(TDSolve)
set(CMAKE_BUILD_TYPE Release)

# The major version number: Number of problems supported
set(TDS_MAJOR_VERSION 2)
# The minor version number: Order of accuracy
set(TDS_MINOR_VERSION 2)
# The growth version number:
set(TDS_GRWTH_VERSION 0)
# The full version
set(TDS_VERSION
    ${TDS_MAJOR_VERSION}.${TDS_MINOR_VERSION}.${TDS_GRWTH_VERSION}
    )

# Bring the headers into the project
include_directories(include)

# Add source files
file(GLOB SOURCES "main/*.cc" "src/*.cc")

# Configure file that pass CMake settings to the source

# Executable
add_executable(TumorDynRun ${SOURCES})
