cmake_minimum_required(VERSION 3.20)
project(nion_tests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to Release build with -O3 and native SIMD if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing -Wno-psabi")
if(NOT MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -mtune=native")
endif()

# Test if C++20 concepts are supported
include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
    "
    #include <concepts>
    int main() {}
    "
    HAS_CONCEPTS
)
if(NOT HAS_CONCEPTS)
    message(FATAL_ERROR "Your compiler does not support concepts. Concepts are required for this library")
endif()

option(TEST_LARGE_DEGREE "Test Sedenion (N=16) and Trigintaduonion (N=32) against Douglas Wilhelm Harder's NumberSystems library" OFF)
if (TEST_LARGE_DEGREE AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(TEST_LARGE_DEGREE OFF)
endif()
if(TEST_LARGE_DEGREE)
    add_compile_definitions(TEST_LARGE_DEGREE)
endif()

option(FULL_RECURSION "Force recursive Cayley-Dickson multiplication" OFF)
if(FULL_RECURSION)
    add_compile_definitions(FULL_RECURSION)
endif()


if(TEST_LARGE_DEGREE)
    # Fetch content from Douglas Wilhelm Harder's NumberSystems library for comparison of our algorithms
    include(FetchContent)
    if(POLICY CMP0135)
        cmake_policy(SET CMP0135 NEW)
    endif()
    if(POLICY CMP0169)
        cmake_policy(SET CMP0169 OLD)
    endif()
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

# Determine target name (avoid collision if included via add_subdirectory from root CMakeLists.txt)
if(TARGET nion)
    set(NION_TEST_TARGET nion_benchmarks)
else()
    set(NION_TEST_TARGET nion)
endif()

add_executable(${NION_TEST_TARGET}
    ../niolloc.hpp
    ../nion_simd.hpp
    ../nion.hpp
    ../nion_trig.hpp
    main.cpp
)

set_target_properties(${NION_TEST_TARGET} PROPERTIES OUTPUT_NAME "nion")
target_include_directories(${NION_TEST_TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

# Locate Boost headers (including Conda environments if present)
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
endif()
find_package(Boost QUIET)
if(Boost_FOUND)
    target_include_directories(${NION_TEST_TARGET} PRIVATE ${Boost_INCLUDE_DIRS})
else()
    find_path(BOOST_HEADER_DIR boost/math/quaternion.hpp
        HINTS
            "$ENV{CONDA_PREFIX}/include"
            "$ENV{HOME}/miniconda3/envs/pdev/include"
            "$ENV{HOME}/miniconda3/include"
            "/usr/include"
            "/usr/local/include"
    )
    if(BOOST_HEADER_DIR)
        message(STATUS "[Nion Tests] Found Boost headers at: ${BOOST_HEADER_DIR}")
        target_include_directories(${NION_TEST_TARGET} PRIVATE ${BOOST_HEADER_DIR})
    endif()
endif()

if(TEST_LARGE_DEGREE)
    file(GLOB_RECURSE number_systems_files ${CMAKE_BINARY_DIR}/_deps/numbersystems/*.cpp)
    add_library(numbersystems_objs OBJECT ${number_systems_files})
    target_include_directories(numbersystems_objs PUBLIC ${CMAKE_BINARY_DIR}/_deps/numbersystems)
    target_link_libraries(${NION_TEST_TARGET} PRIVATE numbersystems_objs)
endif()

# ============================================================================
# Automated CI/CD Accuracy & Benchmark Assertion Executable (nion_ci_test)
# ============================================================================
add_executable(nion_ci_test
    ../niolloc.hpp
    ../nion_simd.hpp
    ../nion.hpp
    ../nion_trig.hpp
    benchmark_ci_test.cpp
)
target_include_directories(nion_ci_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
if(Boost_FOUND)
    target_include_directories(nion_ci_test PRIVATE ${Boost_INCLUDE_DIRS})
elseif(BOOST_HEADER_DIR)
    target_include_directories(nion_ci_test PRIVATE ${BOOST_HEADER_DIR})
endif()

if(TEST_LARGE_DEGREE)
    target_link_libraries(nion_ci_test PRIVATE numbersystems_objs)
endif()

enable_testing()
add_test(NAME NionCI_N2_Complex          COMMAND nion_ci_test --suite=n2)
add_test(NAME NionCI_N4_Quaternion       COMMAND nion_ci_test --suite=n4)
add_test(NAME NionCI_N8_Octonion         COMMAND nion_ci_test --suite=n8)
add_test(NAME NionCI_N16_Sedenion        COMMAND nion_ci_test --suite=n16)
add_test(NAME NionCI_N32_Trigintaduonion COMMAND nion_ci_test --suite=n32)
add_test(NAME NionCI_SIMD_Batches        COMMAND nion_ci_test --suite=simd)
add_test(NAME NionCI_Speed_Benchmarks    COMMAND nion_ci_test --suite=speed)

