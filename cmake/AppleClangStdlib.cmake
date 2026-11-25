# macOS with Apple clang: an outdated /Library/Developer/CommandLineTools/
# usr/include/c++/v1 left behind by an old installation of the Command Line
# Tools hides the C++ standard library of the SDK ("'iostream' file not
# found"). If the compiler cannot find <iostream>, point it to the SDK.
if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  include(CheckCXXSourceCompiles)
  check_cxx_source_compiles("#include <iostream>\nint main(){return 0;}"
                            TDS_HAVE_CXX_STDLIB)
  if(NOT TDS_HAVE_CXX_STDLIB)
    execute_process(COMMAND xcrun --show-sdk-path
                    OUTPUT_VARIABLE TDS_OSX_SDK
                    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(EXISTS "${TDS_OSX_SDK}/usr/include/c++/v1/iostream")
      message(STATUS "Using the C++ standard library headers of ${TDS_OSX_SDK}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -nostdinc++ -isystem ${TDS_OSX_SDK}/usr/include/c++/v1")
    endif()
  endif()
endif()
