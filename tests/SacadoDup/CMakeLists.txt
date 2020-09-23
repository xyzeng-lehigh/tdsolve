cmake_minimum_required(VERSION 3.9)
#set(CMAKE_C_COMPILER "gcc")
#set(CMAKE_CXX_COMPILER "g++")
project(TestSacadoDup)
set(CMAKE_BUILD_TYPE Release)

# Add source files
file(GLOB SOURCES "TestSacadoDup.cc")

# Executable
add_executable(TestSacadoDup ${SOURCES})

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

target_link_libraries(TestSacadoDup ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

