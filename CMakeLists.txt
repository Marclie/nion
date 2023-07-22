cmake_minimum_required(VERSION 3.22)

project(nion VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

# test if concepts are supported
include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
        "
        #include <concepts>
        int main() {}
        "
        HAS_CONCEPTS
)
if ( NOT HAS_CONCEPTS )
    message(FATAL_ERROR "Your compiler does not support concepts. Concepts are required for this library")
endif()


add_library(nion STATIC nion_lib.cpp)

# Optional: If you want to install the library
install(TARGETS nion
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

install(FILES niolloc.hpp nion.hpp nion.cpp DESTINATION include)
