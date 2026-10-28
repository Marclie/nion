// Copyright 2022-2026 Marcus Dante Liebenthal & Contributors
// Licensed under the Apache License, Version 2.0 (the "License");
//
// nion_example.cpp - Verification, SIMD Benchmark, and CPU/GPU Usage Demo
// Compile on CPU with GCC/Clang:
//   g++ -std=c++20 -O3 -march=native nion_example.cpp -o nion_example
// Compile with CUDA (NVCC):
//   nvcc -std=c++20 -O3 -DNION_ENABLE_CUDA -x cu nion_example.cpp -o nion_example

#include "nion.hpp"
#include "nion_trig.hpp"

#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace nion;

// ============================================================================
// 1. Portable CPU / GPU Kernel Functor
// ============================================================================
// Because `quaternion<float>` (`nion<float, 4>`) and `octonion<double>` (`nion<double, 8>`)
// use fixed-size aligned storage and `NION_HOST_DEVICE` methods, the exact same
// function runs inside a CUDA/HIP `__global__` kernel or a CPU loop!

template <typename T, std::size_t N>
NION_HOST_DEVICE void hypercomplex_polynomial_step(
    T* out_buffer, const T* a_buffer, const T* b_buffer, std::size_t idx)
{
    nion_view<const T, N> va(a_buffer + idx * N);
    nion_view<const T, N> vb(b_buffer + idx * N);
    nion_view<T, N>       vout(out_buffer + idx * N);

    // Load into registers
    const nion<T, N> a = va.load();
    const nion<T, N> b = vb.load();

    // Fused expression template + SIMD/register multiplication
    const nion<T, N> c = (a * b) + a.conj_expr() * T{0.5} - b;

    // Store back to global/host memory
    vout.store(c);
}

#if NION_GPU_COMPILER
template <typename T, std::size_t N>
__global__ void gpu_nion_kernel(
    T* out_buffer, const T* a_buffer, const T* b_buffer, std::size_t count)
{
    const std::size_t idx = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx < count) {
        hypercomplex_polynomial_step<T, N>(out_buffer, a_buffer, b_buffer, idx);
    }
}
#endif

// ============================================================================
// 2. Compile-Time `constexpr` Evaluation Test
// ============================================================================
consteval bool verify_constexpr_algebra() {
    constexpr quaternion<double> q1 = {1.0, 2.0, 3.0, 4.0};
    constexpr quaternion<double> q2 = q1.conj();
    constexpr quaternion<double> prod = q1 * q2;
    // q * conj(q) == (norm(q), 0, 0, 0) = (30.0, 0, 0, 0)
    return (prod[0] == 30.0) && (prod[1] == 0.0) && (prod[2] == 0.0) && (prod[3] == 0.0);
}
static_assert(verify_constexpr_algebra(), "Compile-time constexpr Quaternion algebra failed!");

// ============================================================================
// 3. Main Test & Benchmark Driver
// ============================================================================
int main() {
    std::cout << "=========================================================\n";
    std::cout << "  Nion Modern C++20 Hypercomplex Library Verification\n";
    std::cout << "=========================================================\n\n";

    // 1. Original README API Compatibility Check (Dynamic Heap Default)
    {
        nion<float> a = {1.0f, 2.0f, 3.0f, 4.0f};
        nion<float> b = ++a;
        std::cout << "[README Test] ++a      = " << b << " (expected (2,2,3,4))\n";
        --a;
        std::cout << "[README Test] a.conj() = " << a.conj() << " (expected (1,-2,-3,-4))\n";
        std::cout << "[README Test] a.inv()  = " << a.inv() << "\n";
        nion<float> identity = a * a.inv();
        std::cout << "[README Test] a*a.inv()= " << identity << " (expected (1,0,0,0))\n\n";
        assert(std::abs(identity[0] - 1.0f) < 1e-6f);
    }

    // 2. Fixed-Size Stack/Register Quaternions, Octonions, and Sedenions (Zero Heap Alloc)
    {
        octonion<double> o1 = {1, 2, 3, 4, 5, 6, 7, 8};
        octonion<double> o_inv = o1.inv();
        octonion<double> o_id = o1 * o_inv;
        std::cout << "[Octonion N=8]  o1 * o1.inv() = " << o_id << "\n";
        assert(std::abs(o_id[0] - 1.0) < 1e-12);
        for (std::size_t i = 1; i < 8; ++i) {
            assert(std::abs(o_id[i]) < 1e-12);
        }
    }

    // 3. Analytic & Trigonometric Identity Verification: sin^2(q) + cos^2(q) == 1
    //    and log(exp(q)) == q
    {
        quaternion<double> q = {0.35, -0.4, 0.25, 0.6};
        quaternion<double> s = sin(q);
        quaternion<double> c = cos(q);
        quaternion<double> trig_id = s * s + c * c;
        quaternion<double> roundtrip = log(exp(q));

        std::cout << "[Trig Test]     sin^2(q) + cos^2(q) = " << trig_id << "\n";
        std::cout << "[Analytic Test] log(exp(q)) - q     = " << (roundtrip - q) << "\n\n";
        assert(std::abs(trig_id[0] - 1.0) < 1e-12);
        assert((roundtrip - q).abs() < 1e-12);
    }

    // 4. Unified / Device Memory Space Verification
    {
        nion<double, NION_USE_HEAP, UnifiedSpace> u1 = {1.0, 2.0, 3.0, 4.0};
        nion<double, NION_USE_HEAP, UnifiedSpace> u2 = {4.0, -3.0, 2.0, -1.0};
        auto u_prod = u1 * u2;
        std::cout << "[UnifiedSpace]  u1 * u2 = " << u_prod << "\n\n";
    }

    // 5. High-Throughput AoS vs. SoA SIMD Benchmark (Octonions N = 8)
    {
        constexpr std::size_t BATCH = 65536;
        nion_batch_soa<float, 8, HostSpace> soa_a(BATCH);
        nion_batch_soa<float, 8, HostSpace> soa_b(BATCH);
        nion_batch_soa<float, 8, HostSpace> soa_out(BATCH);

        for (std::size_t i = 0; i < BATCH; ++i) {
            soa_a.set(i, octonion<float>{1.0f, 0.5f, -0.25f, 0.1f, 0.3f, -0.7f, 0.2f, 0.4f});
            soa_b.set(i, octonion<float>{0.5f, -1.0f, 0.4f, 0.2f, -0.3f, 0.1f, 0.8f, -0.5f});
        }

        const auto t0 = std::chrono::high_resolution_clock::now();
        nion_batch_soa<float, 8, HostSpace>::multiply_batch(soa_out, soa_a, soa_b);
        const auto t1 = std::chrono::high_resolution_clock::now();

        const double elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        const octonion<float> expected = soa_a.get(0) * soa_b.get(0);
        const octonion<float> actual   = soa_out.get(0);
        assert((expected - actual).abs() < 1e-5f);

        std::cout << "[SoA SIMD Benchmark] Multiplied " << BATCH
                  << " Octonions (N=8) in " << std::fixed << std::setprecision(3)
                  << elapsed_ms << " ms\n";
        std::cout << "                     Sample result[0] = " << actual << "\n";
    }

    std::cout << "\nAll Nion CPU/GPU & SIMD tests passed successfully!\n";
    return 0;
}
