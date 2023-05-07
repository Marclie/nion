cmake_minimum_required(VERSION 3.22)
project(nion)

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

if (TEST_LARGE_DEGREE)
    add_definitions(-DTEST_LARGE_DEGREE)
endif()

if (FULL_RECURSION)
    add_definitions(-DFULL_RECURSION)
endif()

if ( TEST_LARGE_DEGREE )
    # fetch content from Douglas Wilhelm Harder's NumberSystems library for comparison of our algorithms
    include(FetchContent)
    FetchContent_Declare(
            numbersystems
            URL https://ece.uwaterloo.ca/~dwharder/C%2B%2B/CQOST/src/NumberSystems.zip
            DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/_deps
            SOURCE_DIR ${CMAKE_BINARY_DIR}/_deps/numbersystems
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            LOG_DOWNLOAD ON
    )

    FetchContent_GetProperties(numbersystems)
    if(NOT numbersystems_POPULATED)
        FetchContent_Populate(numbersystems)
    endif()
endif()

add_executable(nion ../nion.hpp main.cpp)

if ( TEST_LARGE_DEGREE )
    target_include_directories(nion PRIVATE ${CMAKE_BINARY_DIR}/_deps/numbersystems)
    file(GLOB_RECURSE number_systems_files ${CMAKE_BINARY_DIR}/_deps/numbersystems/*.cpp)
    target_sources(nion PRIVATE ${number_systems_files})
endif()
