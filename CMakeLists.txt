# Copyright 2022-2026 Marcus Dante Liebenthal & Contributors
# Licensed under the Apache License, Version 2.0 (the "License");
#
# CMakeLists.txt - Modern CMake (3.20+) Build System for Nion
#
# Usage Examples:
#   1. Standard CPU Build with Native SIMD (AVX2 / AVX-512 / ARM NEON):
#      cmake -B build -DCMAKE_BUILD_TYPE=Release
#      cmake --build build -j
#      ctest --test-dir build --output-on-failure
#
#   2. NVIDIA CUDA GPU Build:
#      cmake -B build_cuda -DCMAKE_BUILD_TYPE=Release -DNION_ENABLE_CUDA=ON
#      cmake --build build_cuda -j
#
#   3. AMD ROCm / HIP GPU Build:
#      cmake -B build_hip -DCMAKE_BUILD_TYPE=Release -DNION_ENABLE_HIP=ON
#      cmake --build build_hip -j

cmake_minimum_required(VERSION 3.20)

project(nion
    VERSION 2.0.0
    DESCRIPTION "High-Performance CPU/GPU & SIMD C++20 Hypercomplex Arithmetic Library"
    HOMEPAGE_URL "https://github.com/Marclie/nion"
    LANGUAGES CXX
)

# ============================================================================
# 1. Build Options & Configuration Switches
# ============================================================================
# Determine if nion is the top-level project or included via add_subdirectory/FetchContent
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(NION_IS_TOP_LEVEL ON)
else()
    set(NION_IS_TOP_LEVEL OFF)
endif()

option(NION_ENABLE_NATIVE_SIMD "Enable -march=native / /arch:AVX2 for maximum CPU SIMD speed" ON)
option(NION_ENABLE_FAST_MATH   "Enable aggressive floating-point optimizations (-ffast-math / /fp:fast)" OFF)
option(NION_ENABLE_CUDA        "Enable NVIDIA CUDA GPU memory spaces and kernel support" OFF)
option(NION_ENABLE_HIP         "Enable AMD ROCm / HIP GPU memory spaces and kernel support" OFF)
option(NION_BUILD_EXAMPLES     "Build the Nion verification & SIMD/GPU benchmark executable" ${NION_IS_TOP_LEVEL})
option(NION_BUILD_TESTS        "Register verification suite with CTest" ${NION_IS_TOP_LEVEL})
option(NION_INSTALL            "Generate install and package export targets" ${NION_IS_TOP_LEVEL})

if(NION_ENABLE_CUDA AND NION_ENABLE_HIP)
    message(FATAL_ERROR "NION_ENABLE_CUDA and NION_ENABLE_HIP cannot be enabled simultaneously.")
endif()

# Default to Release build if not specified on single-config generators
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "[Nion] Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# ============================================================================
# 2. Core Header-Only Interface Library Target (nion::nion)
# ============================================================================
add_library(nion INTERFACE)
add_library(nion::nion ALIAS nion)

# Require C++20 for Concepts, std::bit, constexpr algorithms, and aligned new
target_compile_features(nion INTERFACE cxx_std_20)

target_include_directories(nion INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

# Provide an IDE-friendly header list for Visual Studio / CLion / Xcode
set(NION_HEADERS
    niolloc.hpp
    nion_simd.hpp
    nion.hpp
    nion_trig.hpp
)

# In CMake >= 3.19, INTERFACE libraries can have PRIVATE sources to populate IDE project trees 
# without improperly propagating them as compilation steps to consumers.
target_sources(nion PRIVATE ${NION_HEADERS})

# ============================================================================
# 3. Compiler Warnings & CPU SIMD Architecture Flags
# ============================================================================
include(CheckCXXCompilerFlag)

# Compiler warning flags (applied only when building examples/tests)
add_library(nion_compiler_flags INTERFACE)

if(MSVC)
    target_compile_options(nion_compiler_flags INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:/W4 /permissive- /Zc:__cplusplus /EHsc>
    )
else()
    target_compile_options(nion_compiler_flags INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -Wno-unknown-pragmas>
    )
endif()

# Native CPU SIMD Dispatch Flags (AVX-512 / AVX2+FMA / SSE2 / ARM NEON)
if(NION_ENABLE_NATIVE_SIMD)
    if(MSVC)
        check_cxx_compiler_flag("/arch:AVX2" NION_COMPILER_SUPPORTS_AVX2)
        if(NION_COMPILER_SUPPORTS_AVX2)
            target_compile_options(nion INTERFACE
                $<$<COMPILE_LANGUAGE:CXX>:/arch:AVX2>
            )
            message(STATUS "[Nion] Enabled MSVC SIMD flag: /arch:AVX2")
        endif()
    else()
        check_cxx_compiler_flag("-march=native" NION_COMPILER_SUPPORTS_MARCH_NATIVE)
        if(NION_COMPILER_SUPPORTS_MARCH_NATIVE)
            target_compile_options(nion INTERFACE
                $<$<COMPILE_LANGUAGE:CXX>:-march=native -mtune=native>
                $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-march=native>
            )
            message(STATUS "[Nion] Enabled native CPU SIMD instructions (-march=native)")
        else()
            # Apple Silicon / ARM64 fallback where -mcpu=native is preferred over -march=native
            check_cxx_compiler_flag("-mcpu=native" NION_COMPILER_SUPPORTS_MCPU_NATIVE)
            if(NION_COMPILER_SUPPORTS_MCPU_NATIVE)
                target_compile_options(nion INTERFACE
                    $<$<COMPILE_LANGUAGE:CXX>:-mcpu=native>
                )
                message(STATUS "[Nion] Enabled native CPU SIMD instructions (-mcpu=native)")
            endif()
        endif()
    endif()
endif()

# Optional Fast-Math FMA & Reciprocal Optimizations
if(NION_ENABLE_FAST_MATH)
    if(MSVC)
        target_compile_options(nion INTERFACE $<$<COMPILE_LANGUAGE:CXX>:/fp:fast>)
    else()
        target_compile_options(nion INTERFACE
            $<$<COMPILE_LANGUAGE:CXX>:-ffast-math>
            $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math>
        )
    endif()
endif()

# ============================================================================
# 4. Optional GPU Backends (NVIDIA CUDA / AMD ROCm HIP)
# ============================================================================

# --- 4a. NVIDIA CUDA Backend ---
if(NION_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(NOT CMAKE_CUDA_COMPILER)
        message(FATAL_ERROR "NION_ENABLE_CUDA is ON, but no CUDA compiler (nvcc) was found.")
    endif()

    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)

    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES "native" CACHE STRING "Target CUDA GPU architectures" FORCE)
    endif()

    target_compile_definitions(nion INTERFACE NION_ENABLE_CUDA=1)
    target_link_libraries(nion INTERFACE CUDA::cudart)

    # Enable constexpr host/device interoperability in NVCC
    target_compile_options(nion INTERFACE
        $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr --extended-lambda>
    )

    message(STATUS "[Nion] CUDA GPU backend enabled (Architectures: ${CMAKE_CUDA_ARCHITECTURES})")
endif()

# --- 4b. AMD ROCm / HIP Backend ---
if(NION_ENABLE_HIP)
    find_package(hip REQUIRED)
    target_compile_definitions(nion INTERFACE NION_ENABLE_HIP=1)
    target_link_libraries(nion INTERFACE hip::host hip::device)
    message(STATUS "[Nion] AMD ROCm / HIP GPU backend enabled")
endif()

# ============================================================================
# 5. Verification Suite, SIMD Benchmark, and GPU Kernel Demo Executable
# ============================================================================
if(NION_BUILD_EXAMPLES)
    add_executable(nion_example nion_example.cpp)
    target_link_libraries(nion_example PRIVATE nion::nion nion_compiler_flags)

    # When CUDA is enabled, compile nion_example.cpp through NVCC so __global__ kernels run on the GPU
    if(NION_ENABLE_CUDA)
        set_source_files_properties(nion_example.cpp PROPERTIES LANGUAGE CUDA)
        set_target_properties(nion_example PROPERTIES
            CUDA_STANDARD 20
            CUDA_STANDARD_REQUIRED ON
            CUDA_SEPARABLE_COMPILATION ON
        )
    else()
        set_target_properties(nion_example PROPERTIES
            CXX_STANDARD 20
            CXX_STANDARD_REQUIRED ON
            CXX_EXTENSIONS OFF
        )
    endif()
endif()

# ============================================================================
# 6. CTest Integration
# ============================================================================
if(NION_BUILD_TESTS AND NION_BUILD_EXAMPLES)
    enable_testing()
    add_test(
        NAME NionVerificationAndSimdSuite
        COMMAND nion_example
    )
endif()

# ============================================================================
# 7. Installation & CMake Package Export (for Downstream Projects like nionix)
# ============================================================================
if(NION_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    # Install headers
    install(FILES ${NION_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    # Install interface target
    install(TARGETS nion
        EXPORT nionTargets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    # Export targets for in-tree and installed find_package(nion) support
    install(EXPORT nionTargets
        FILE nionTargets.cmake
        NAMESPACE nion::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nion
    )

    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/nionConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
        ARCH_INDEPENDENT
    )

    # Generate a minimal self-contained nionConfig.cmake
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/nionConfig.cmake"
        "include(CMakeFindDependencyMacro)\n"
        "if(${NION_ENABLE_CUDA})\n"
        "    find_dependency(CUDAToolkit REQUIRED)\n"
        "endif()\n"
        "if(${NION_ENABLE_HIP})\n"
        "    find_dependency(hip REQUIRED)\n"
        "endif()\n"
        "include(\"\${CMAKE_CURRENT_LIST_DIR}/nionTargets.cmake\")\n"
    )

    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/nionConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/nionConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nion
    )
endif()
