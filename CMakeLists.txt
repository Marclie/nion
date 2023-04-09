cmake_minimum_required(VERSION 3.22)

project(nion VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(nion STATIC nion_lib.cpp)

# Optional: If you want to install the library
install(TARGETS nion
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

install(FILES nion.hpp nion.cpp DESTINATION include)
