cmake_minimum_required(VERSION 3.22)
project(nion)

set(CMAKE_CXX_STANDARD 17)


if ( TEST_LARGE_DEGREE )
    # fetch content from Douglas Wilhelm Harder's NumberSystems library for comparison of our algorithms
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")
    include(FetchContent)
    FetchContent_Declare(
            numbersystems
            URL https://ece.uwaterloo.ca/~dwharder/C++/CQOST/src/NumberSystems.zip
            DOWNLOAD_EXTRACT_TIMESTAMP NEW
    )

    FetchContent_GetProperties(numbersystems)
    if(NOT numbersystems_POPULATED)
        FetchContent_Populate(numbersystems)
    endif()
endif()

add_executable(nion nion.hpp main.cpp)

if ( TEST_LARGE_DEGREE )
    target_include_directories(nion PRIVATE ${numbersystems_SOURCE_DIR})
    file(GLOB_RECURSE number_systems_files ${numbersystems_SOURCE_DIR}/*.cpp)
    target_sources(nion PRIVATE ${number_systems_files})
endif()
