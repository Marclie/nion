// Copyright 2022-2026 Marcus Dante Liebenthal & Contributors
// Licensed under the Apache License, Version 2.0 (the "License");
//
// nion_simd.hpp - Compile-Time Cayley-Dickson Basis Engine & Hardware SIMD Dispatch
// Resolves the Tension Between Heavy Templating and Hardware SIMD Speed:
//   - Any arbitrary algebraic type T uses constexpr / GPU-ready unrolled templates.
//   - Primitive float / double types on CPU automatically dispatch to AVX-512,
//     AVX2+FMA, SSE2, or ARM NEON intrinsics via `if constexpr` & specialization.

#ifndef NION_SIMD_HPP
#define NION_SIMD_HPP

#include "niolloc.hpp"
#include <array>

// Detect CPU SIMD ISA Extensions (Disabled when compiling device code on GPU)
#if !NION_ON_DEVICE
    #if defined(__AVX512F__)
        #include <immintrin.h>
        #define NION_SIMD_AVX512 1
        #define NION_SIMD_AVX2   1
        #define NION_SIMD_SSE2   1
    #elif defined(__AVX2__)
        #include <immintrin.h>
        #define NION_SIMD_AVX512 0
        #define NION_SIMD_AVX2   1
        #define NION_SIMD_SSE2   1
    #elif defined(__SSE2__) || (defined(_MSC_VER) && (defined(_M_X64) || _M_IX86_FP >= 2))
        #include <emmintrin.h>
        #if defined(__SSE3__)
            #include <pmmintrin.h>
        #endif
        #define NION_SIMD_AVX512 0
        #define NION_SIMD_AVX2   0
        #define NION_SIMD_SSE2   1
    #else
        #define NION_SIMD_AVX512 0
        #define NION_SIMD_AVX2   0
        #define NION_SIMD_SSE2   0
    #endif

    #if defined(__ARM_NEON) || defined(__ARM_NEON__)
        #include <arm_neon.h>
        #define NION_SIMD_NEON 1
    #else
        #define NION_SIMD_NEON 0
    #endif
#else
    #define NION_SIMD_AVX512 0
    #define NION_SIMD_AVX2   0
    #define NION_SIMD_SSE2   0
    #define NION_SIMD_NEON   0
#endif

namespace nion_math::simd {

// ============================================================================
// 1. Compile-Time O(log N) Bitwise & O(1) Table-Driven Cayley-Dickson Signs
// ============================================================================
// In the Cayley-Dickson algebra of dimension 2^k, basis elements e_p and e_q
// multiply as:
//     e_p * e_q = sign(p, q) * e_{p ^ q}
// where `p ^ q` is the bitwise XOR of indices p and q, and `sign(p, q)` in {-1, +1}
// is determined recursively by the Cayley-Dickson doubling rule:
//     (a, b) * (c, d) = (a*c - conj(d)*b, d*a + b*conj(c))

/// Bitwise constexpr computation of the Cayley-Dickson multiplication sign:
/// Returns +1 or -1 for basis product e_p * e_q = (+/- 1) * e_{p ^ q}.
/// Works for ANY arbitrary dimension on both Host and GPU Device!
NION_HOST_DEVICE constexpr int cd_sign_bitwise(std::size_t p, std::size_t q) noexcept {
    if (p == 0 || q == 0) return 1;
    if (p == q) return -1;

    // Find highest power of two splitting p and q
    const std::size_t max_pq = (p > q) ? p : q;
    std::size_t half = 1;
    while ((half << 1) <= max_pq) {
        half <<= 1;
    }

    const bool p_hi = (p >= half);
    const bool q_hi = (q >= half);
    const std::size_t p_rem = p & (half - 1);
    const std::size_t q_rem = q & (half - 1);

    if (!p_hi && q_hi) {
        // (a, 0) * (0, d) = (0, d * a) -> sign(q_rem, p_rem)
        return cd_sign_bitwise(q_rem, p_rem);
    } else if (p_hi && !q_hi) {
        // (0, b) * (c, 0) = (0, b * conj(c))
        const int conj_c = (q_rem == 0) ? 1 : -1;
        return conj_c * cd_sign_bitwise(p_rem, q_rem);
    } else {
        // Both p_hi and q_hi: (0, b) * (0, d) = (-conj(d) * b, 0)
        const int conj_d = (q_rem == 0) ? 1 : -1;
        return -conj_d * cd_sign_bitwise(q_rem, p_rem);
    }
}

/// Precomputed Compile-Time Sign Table for fixed dimensions up to N = 64.
/// Eliminates all recursion at runtime!
/// For output index `k = p ^ q`, `e_p * e_{p ^ k} = row_sign(k, p) * e_k`.
template <std::size_t N>
struct CdSignTable {
    // signs[k][p] stores the sign of e_p * e_{p ^ k} contributing to output component e_k
    std::array<std::array<int8_t, N>, N> signs{};

    constexpr CdSignTable() noexcept {
        for (std::size_t k = 0; k < N; ++k) {
            for (std::size_t p = 0; p < N; ++p) {
                const std::size_t q = p ^ k;
                signs[k][p] = static_cast<int8_t>(cd_sign_bitwise(p, q));
            }
        }
    }
};

template <std::size_t N>
inline constexpr CdSignTable<N> cd_table{};

// ============================================================================
// 2. Vectorized Element-Wise Kernels (Add, Sub, Scale, Conj, NormSq)
// ============================================================================

template <typename T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr void vec_add(
    T* NION_RESTRICT dst, const T* NION_RESTRICT a, const T* NION_RESTRICT b, std::size_t n) noexcept
{
#if !NION_ON_DEVICE
    if (!std::is_constant_evaluated()) {
    #if NION_SIMD_AVX2
        if constexpr (std::is_same_v<T, float>) {
            std::size_t i = 0;
            for (; i + 8 <= n; i += 8) {
                __m256 va = _mm256_loadu_ps(a + i);
                __m256 vb = _mm256_loadu_ps(b + i);
                _mm256_storeu_ps(dst + i, _mm256_add_ps(va, vb));
            }
            for (; i < n; ++i) dst[i] = a[i] + b[i];
            return;
        } else if constexpr (std::is_same_v<T, double>) {
            std::size_t i = 0;
            for (; i + 4 <= n; i += 4) {
                __m256d va = _mm256_loadu_pd(a + i);
                __m256d vb = _mm256_loadu_pd(b + i);
                _mm256_storeu_pd(dst + i, _mm256_add_pd(va, vb));
            }
            for (; i < n; ++i) dst[i] = a[i] + b[i];
            return;
        }
    #elif NION_SIMD_SSE2
        if constexpr (std::is_same_v<T, float>) {
            std::size_t i = 0;
            for (; i + 4 <= n; i += 4) {
                __m128 va = _mm_loadu_ps(a + i);
                __m128 vb = _mm_loadu_ps(b + i);
                _mm_storeu_ps(dst + i, _mm_add_ps(va, vb));
            }
            for (; i < n; ++i) dst[i] = a[i] + b[i];
            return;
        } else if constexpr (std::is_same_v<T, double>) {
            std::size_t i = 0;
            for (; i + 2 <= n; i += 2) {
                __m128d va = _mm_loadu_pd(a + i);
                __m128d vb = _mm_loadu_pd(b + i);
                _mm_storeu_pd(dst + i, _mm_add_pd(va, vb));
            }
            for (; i < n; ++i) dst[i] = a[i] + b[i];
            return;
        }
    #elif NION_SIMD_NEON
        if constexpr (std::is_same_v<T, float>) {
            std::size_t i = 0;
            for (; i + 4 <= n; i += 4) {
                float32x4_t va = vld1q_f32(a + i);
                float32x4_t vb = vld1q_f32(b + i);
                vst1q_f32(dst + i, vaddq_f32(va, vb));
            }
            for (; i < n; ++i) dst[i] = a[i] + b[i];
            return;
        }
    #endif
    }
#endif
    for (std::size_t i = 0; i < n; ++i) {
        dst[i] = a[i] + b[i];
    }
}

template <typename T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr void vec_sub(
    T* NION_RESTRICT dst, const T* NION_RESTRICT a, const T* NION_RESTRICT b, std::size_t n) noexcept
{
#if !NION_ON_DEVICE
    if (!std::is_constant_evaluated()) {
    #if NION_SIMD_AVX2
        if constexpr (std::is_same_v<T, float>) {
            std::size_t i = 0;
            for (; i + 8 <= n; i += 8) {
                __m256 va = _mm256_loadu_ps(a + i);
                __m256 vb = _mm256_loadu_ps(b + i);
                _mm256_storeu_ps(dst + i, _mm256_sub_ps(va, vb));
            }
            for (; i < n; ++i) dst[i] = a[i] - b[i];
            return;
        } else if constexpr (std::is_same_v<T, double>) {
            std::size_t i = 0;
            for (; i + 4 <= n; i += 4) {
                __m256d va = _mm256_loadu_pd(a + i);
                __m256d vb = _mm256_loadu_pd(b + i);
                _mm256_storeu_pd(dst + i, _mm256_sub_pd(va, vb));
            }
            for (; i < n; ++i) dst[i] = a[i] - b[i];
            return;
        }
    #elif NION_SIMD_SSE2
        if constexpr (std::is_same_v<T, float>) {
            std::size_t i = 0;
            for (; i + 4 <= n; i += 4) {
                __m128 va = _mm_loadu_ps(a + i);
                __m128 vb = _mm_loadu_ps(b + i);
                _mm_storeu_ps(dst + i, _mm_sub_ps(va, vb));
            }
            for (; i < n; ++i) dst[i] = a[i] - b[i];
            return;
        } else if constexpr (std::is_same_v<T, double>) {
            std::size_t i = 0;
            for (; i + 2 <= n; i += 2) {
                __m128d va = _mm_loadu_pd(a + i);
                __m128d vb = _mm_loadu_pd(b + i);
                _mm_storeu_pd(dst + i, _mm_sub_pd(va, vb));
            }
            for (; i < n; ++i) dst[i] = a[i] - b[i];
            return;
        }
    #endif
    }
#endif
    for (std::size_t i = 0; i < n; ++i) {
        dst[i] = a[i] - b[i];
    }
}

template <typename T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr void vec_scale(
    T* NION_RESTRICT dst, const T* NION_RESTRICT a, const T& s, std::size_t n) noexcept
{
#if !NION_ON_DEVICE
    if (!std::is_constant_evaluated()) {
    #if NION_SIMD_AVX2
        if constexpr (std::is_same_v<T, float>) {
            std::size_t i = 0;
            __m256 vs = _mm256_set1_ps(s);
            for (; i + 8 <= n; i += 8) {
                __m256 va = _mm256_loadu_ps(a + i);
                _mm256_storeu_ps(dst + i, _mm256_mul_ps(va, vs));
            }
            for (; i < n; ++i) dst[i] = a[i] * s;
            return;
        } else if constexpr (std::is_same_v<T, double>) {
            std::size_t i = 0;
            __m256d vs = _mm256_set1_pd(s);
            for (; i + 4 <= n; i += 4) {
                __m256d va = _mm256_loadu_pd(a + i);
                _mm256_storeu_pd(dst + i, _mm256_mul_pd(va, vs));
            }
            for (; i < n; ++i) dst[i] = a[i] * s;
            return;
        }
    #elif NION_SIMD_SSE2
        if constexpr (std::is_same_v<T, float>) {
            std::size_t i = 0;
            __m128 vs = _mm_set1_ps(s);
            for (; i + 4 <= n; i += 4) {
                __m128 va = _mm_loadu_ps(a + i);
                _mm_storeu_ps(dst + i, _mm_mul_ps(va, vs));
            }
            for (; i < n; ++i) dst[i] = a[i] * s;
            return;
        } else if constexpr (std::is_same_v<T, double>) {
            std::size_t i = 0;
            __m128d vs = _mm_set1_pd(s);
            for (; i + 2 <= n; i += 2) {
                __m128d va = _mm_loadu_pd(a + i);
                _mm_storeu_pd(dst + i, _mm_mul_pd(va, vs));
            }
            for (; i < n; ++i) dst[i] = a[i] * s;
            return;
        }
    #endif
    }
#endif
    for (std::size_t i = 0; i < n; ++i) {
        dst[i] = a[i] * s;
    }
}

/// Hypercomplex Conjugation: a_0 remains positive, a_1 ... a_{n-1} are negated
template <typename T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr void vec_conj(
    T* NION_RESTRICT dst, const T* NION_RESTRICT a, std::size_t n) noexcept
{
    if (n == 0) return;
    dst[0] = a[0];
    for (std::size_t i = 1; i < n; ++i) {
        dst[i] = -a[i];
    }
}

/// Squared Euclidean / Cayley-Dickson Norm: sum_{i=0}^{n-1} a_i^2
template <typename T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr T vec_norm_sq(
    const T* NION_RESTRICT a, std::size_t n) noexcept
{
    if (n == 0) return T{};
    if (n == 1) return a[0] * a[0];
    if (n == 2) return a[0] * a[0] + a[1] * a[1];
    if (n == 3) return (a[0] * a[0] + a[1] * a[1]) + a[2] * a[2];
    if (n == 4) return (a[0] * a[0] + a[1] * a[1]) + (a[2] * a[2] + a[3] * a[3]);

#if !NION_ON_DEVICE
    if (!std::is_constant_evaluated()) {
    #if NION_SIMD_AVX2
        if constexpr (std::is_same_v<T, float>) {
            std::size_t i = 0;
            __m256 acc = _mm256_setzero_ps();
            for (; i + 8 <= n; i += 8) {
                __m256 va = _mm256_loadu_ps(a + i);
            #if defined(__FMA__)
                acc = _mm256_fmadd_ps(va, va, acc);
            #else
                acc = _mm256_add_ps(acc, _mm256_mul_ps(va, va));
            #endif
            }
            __m128 lo = _mm256_castps256_ps128(acc);
            __m128 hi = _mm256_extractf128_ps(acc, 1);
            __m128 sum128 = _mm_add_ps(lo, hi);
            __m128 shuf = _mm_movehdup_ps(sum128);
            __m128 sums = _mm_add_ps(sum128, shuf);
            shuf = _mm_movehl_ps(shuf, sums);
            sums = _mm_add_ss(sums, shuf);
            float total = _mm_cvtss_f32(sums);
            for (; i < n; ++i) total += a[i] * a[i];
            return total;
        } else if constexpr (std::is_same_v<T, double>) {
            std::size_t i = 0;
            __m256d acc = _mm256_setzero_pd();
            for (; i + 4 <= n; i += 4) {
                __m256d va = _mm256_loadu_pd(a + i);
            #if defined(__FMA__)
                acc = _mm256_fmadd_pd(va, va, acc);
            #else
                acc = _mm256_add_pd(acc, _mm256_mul_pd(va, va));
            #endif
            }
            __m128d lo = _mm256_castpd256_pd128(acc);
            __m128d hi = _mm256_extractf128_pd(acc, 1);
            __m128d sum128 = _mm_add_pd(lo, hi);
            __m128d hi64 = _mm_unpackhi_pd(sum128, sum128);
            double total = _mm_cvtsd_f64(_mm_add_sd(sum128, hi64));
            for (; i < n; ++i) total += a[i] * a[i];
            return total;
        }
    #endif
    }
#endif
    T sum = a[0] * a[0];
    for (std::size_t i = 1; i < n; ++i) {
        sum = sum + a[i] * a[i];
    }
    return sum;
}

// ============================================================================
// 3. Fast Hypercomplex Multiplication Kernels (Specialized + SIMD + Generic)
// ============================================================================
// Standard Cayley-Dickson recurrence:
//   (a, b) * (c, d) = (a*c - conj(d)*b, d*a + b*conj(c))
// Verify for N=2 (Complex):
//   (a0, a1) * (b0, b1) = (a0*b0 - b1*a1, b1*a0 + a1*b0)
// Verify for N=4 (Quaternions):
//   e0 = a0*b0 - a1*b1 - a2*b2 - a3*b3
//   e1 = a0*b1 + a1*b0 + a2*b3 - a3*b2
//   e2 = a0*b2 - a1*b3 + a2*b0 + a3*b1
//   e3 = a0*b3 + a1*b2 - a2*b1 + a3*b0

/// Unrolled N = 2 (Complex) multiplication
template <typename T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr void mul_fixed_2(
    T* NION_RESTRICT out, const T* NION_RESTRICT a, const T* NION_RESTRICT b) noexcept
{
    const T a0 = a[0], a1 = a[1];
    const T b0 = b[0], b1 = b[1];
    out[0] = a0 * b0 - a1 * b1;
    out[1] = a0 * b1 + a1 * b0;
}

/// Unrolled & SIMD-Accelerated N = 4 (Quaternion) multiplication
template <typename T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr void mul_fixed_4(
    T* NION_RESTRICT out, const T* NION_RESTRICT a, const T* NION_RESTRICT b) noexcept
{
#if !NION_ON_DEVICE
    if (!std::is_constant_evaluated()) {
    #if NION_SIMD_SSE2
        if constexpr (std::is_same_v<T, float>) {
            // 128-bit SSE Quaternion Multiplication
            const __m128 va = _mm_loadu_ps(a);
            const __m128 vb = _mm_loadu_ps(b);

            const __m128 a0 = _mm_shuffle_ps(va, va, _MM_SHUFFLE(0, 0, 0, 0));
            const __m128 a1 = _mm_shuffle_ps(va, va, _MM_SHUFFLE(1, 1, 1, 1));
            const __m128 a2 = _mm_shuffle_ps(va, va, _MM_SHUFFLE(2, 2, 2, 2));
            const __m128 a3 = _mm_shuffle_ps(va, va, _MM_SHUFFLE(3, 3, 3, 3));

            const __m128 b_shuf1 = _mm_shuffle_ps(vb, vb, _MM_SHUFFLE(2, 3, 0, 1)); // (b1, b0, b3, b2)
            const __m128 b_shuf2 = _mm_shuffle_ps(vb, vb, _MM_SHUFFLE(1, 0, 3, 2)); // (b2, b3, b0, b1)
            const __m128 b_shuf3 = _mm_shuffle_ps(vb, vb, _MM_SHUFFLE(0, 1, 2, 3)); // (b3, b2, b1, b0)

            const __m128 s1 = _mm_setr_ps(-1.0f,  1.0f, -1.0f,  1.0f);
            const __m128 s2 = _mm_setr_ps(-1.0f,  1.0f,  1.0f, -1.0f);
            const __m128 s3 = _mm_setr_ps(-1.0f, -1.0f,  1.0f,  1.0f);

            __m128 res = _mm_mul_ps(a0, vb);
        #if defined(__FMA__)
            res = _mm_fmadd_ps(_mm_mul_ps(a1, s1), b_shuf1, res);
            res = _mm_fmadd_ps(_mm_mul_ps(a2, s2), b_shuf2, res);
            res = _mm_fmadd_ps(_mm_mul_ps(a3, s3), b_shuf3, res);
        #else
            res = _mm_add_ps(res, _mm_mul_ps(_mm_mul_ps(a1, s1), b_shuf1));
            res = _mm_add_ps(res, _mm_mul_ps(_mm_mul_ps(a2, s2), b_shuf2));
            res = _mm_add_ps(res, _mm_mul_ps(_mm_mul_ps(a3, s3), b_shuf3));
        #endif

            _mm_storeu_ps(out, res);
            return;
        }
    #endif
    #if NION_SIMD_AVX2
        if constexpr (std::is_same_v<T, double>) {
            // 256-bit AVX2 Double-Precision Quaternion Multiplication
            const __m256d va = _mm256_loadu_pd(a);
            const __m256d vb = _mm256_loadu_pd(b);

            const __m256d a0 = _mm256_permute4x64_pd(va, _MM_SHUFFLE(0, 0, 0, 0));
            const __m256d a1 = _mm256_permute4x64_pd(va, _MM_SHUFFLE(1, 1, 1, 1));
            const __m256d a2 = _mm256_permute4x64_pd(va, _MM_SHUFFLE(2, 2, 2, 2));
            const __m256d a3 = _mm256_permute4x64_pd(va, _MM_SHUFFLE(3, 3, 3, 3));

            const __m256d b_shuf1 = _mm256_permute4x64_pd(vb, _MM_SHUFFLE(2, 3, 0, 1));
            const __m256d b_shuf2 = _mm256_permute4x64_pd(vb, _MM_SHUFFLE(1, 0, 3, 2));
            const __m256d b_shuf3 = _mm256_permute4x64_pd(vb, _MM_SHUFFLE(0, 1, 2, 3));

            const __m256d s1 = _mm256_setr_pd(-1.0,  1.0, -1.0,  1.0);
            const __m256d s2 = _mm256_setr_pd(-1.0,  1.0,  1.0, -1.0);
            const __m256d s3 = _mm256_setr_pd(-1.0, -1.0,  1.0,  1.0);

            __m256d res = _mm256_mul_pd(a0, vb);
        #if defined(__FMA__)
            res = _mm256_fmadd_pd(_mm256_mul_pd(a1, s1), b_shuf1, res);
            res = _mm256_fmadd_pd(_mm256_mul_pd(a2, s2), b_shuf2, res);
            res = _mm256_fmadd_pd(_mm256_mul_pd(a3, s3), b_shuf3, res);
        #else
            res = _mm256_add_pd(res, _mm256_mul_pd(_mm256_mul_pd(a1, s1), b_shuf1));
            res = _mm256_add_pd(res, _mm256_mul_pd(_mm256_mul_pd(a2, s2), b_shuf2));
            res = _mm256_add_pd(res, _mm256_mul_pd(_mm256_mul_pd(a3, s3), b_shuf3));
        #endif

            _mm256_storeu_pd(out, res);
            return;
        }
    #endif
    }
#endif
    const T a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
    const T b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];

    out[0] = a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3;
    out[1] = a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2;
    out[2] = a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1;
    out[3] = a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0;
}

/// Helper: update all output lanes K = 0..N-1 for a fixed input lane p = P + 1 (SIMD-friendly).
template <std::size_t N, std::size_t P, typename T, std::size_t... K>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr void mul_fixed_step_p(
    T* NION_RESTRICT acc, const T* NION_RESTRICT a, const T* NION_RESTRICT b,
    std::index_sequence<K...>) noexcept
{
    constexpr std::size_t p = P + 1;
    const T ap = a[p];
    ((cd_table<N>.signs[K][p] > 0
        ? (acc[K] = acc[K] + ap * b[p ^ K])
        : (acc[K] = acc[K] - ap * b[p ^ K])), ...);
}

/// Helper: compute single output lane K across P = 0..N-2 in one scalar accumulator (x87-friendly).
template <std::size_t N, std::size_t K, typename T, std::size_t... P>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr T mul_fixed_lane_k(
    const T* NION_RESTRICT a, const T* NION_RESTRICT b,
    std::index_sequence<P...>) noexcept
{
    T acc = a[0] * b[K];
    ((cd_table<N>.signs[K][P + 1] > 0
        ? (acc = acc + a[P + 1] * b[(P + 1) ^ K])
        : (acc = acc - a[P + 1] * b[(P + 1) ^ K])), ...);
    return acc;
}

/// Compile-time unrolled multiplication for N = 8, 16, 32.
/// - For SIMD floating-point types (`float`, `double`), iterates `p` in the outer pack and updates
///   all `K = 0..N-1` independent output lanes in parallel in YMM/XMM registers (AVX2 FMA).
/// - For x87 `long double` or composite types, iterates `K` in the outer pack and accumulates a single
///   register `acc` (`st(0)`) across `p = 0..N-1` to avoid x87 register-stack spills!
template <std::size_t N, typename T, std::size_t... K, std::size_t... P>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr void mul_fixed_unrolled_impl(
    T* NION_RESTRICT out, const T* NION_RESTRICT a, const T* NION_RESTRICT b,
    std::index_sequence<K...> ks, std::index_sequence<P...> ps) noexcept
{
    if constexpr (concepts::SimdFloat<T>) {
        T acc[N]{ (a[0] * b[K])... };
        (mul_fixed_step_p<N, P>(acc, a, b, ks), ...);
        ((out[K] = acc[K]), ...);
    } else {
        ((out[K] = mul_fixed_lane_k<N, K>(a, b, ps)), ...);
    }
}

/// Compile-Time Table-Driven Multiplication for Fixed Dimensions N <= 64 (e.g., Octonions N=8, Sedenions N=16)
/// Uses the constexpr `cd_table<N>` so the compiler completely unrolls the dot products without branches!
template <std::size_t N, typename T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr void mul_fixed_table(
    T* NION_RESTRICT out, const T* NION_RESTRICT a, const T* NION_RESTRICT b) noexcept
{
    if constexpr (N == 1) {
        out[0] = a[0] * b[0];
    } else if constexpr (N == 2) {
        mul_fixed_2(out, a, b);
    } else if constexpr (N == 4) {
        mul_fixed_4(out, a, b);
    } else if constexpr (N <= 32) {
        mul_fixed_unrolled_impl<N>(out, a, b, std::make_index_sequence<N>{}, std::make_index_sequence<N - 1>{});
    } else if constexpr (N <= 64) {
        T tmp[N];
        for (std::size_t k = 0; k < N; ++k) {
            tmp[k] = a[0] * b[k];
        }
        for (std::size_t p = 1; p < N; ++p) {
            const T ap = a[p];
            for (std::size_t k = 0; k < N; ++k) {
                const std::size_t q = p ^ k;
                if (cd_table<N>.signs[k][p] > 0) {
                    tmp[k] = tmp[k] + ap * b[q];
                } else {
                    tmp[k] = tmp[k] - ap * b[q];
                }
            }
        }
        for (std::size_t k = 0; k < N; ++k) {
            out[k] = tmp[k];
        }
    } else {
        // For very large fixed N > 64, use two-level blocked Cayley-Dickson doubling
        constexpr std::size_t half = N / 2;
        // (a, b_hi) * (c, d_hi) = (a*c - conj(d_hi)*b_hi, d_hi*a + b_hi*conj(c))
        T ac[half], conj_d[half], conj_d_b[half], da[half], conj_c[half], b_conj_c[half];
        mul_fixed_table<half>(ac, a, b);
        vec_conj(conj_d, b + half, half);
        mul_fixed_table<half>(conj_d_b, conj_d, a + half);
        mul_fixed_table<half>(da, b + half, a);
        vec_conj(conj_c, b, half);
        mul_fixed_table<half>(b_conj_c, a + half, conj_c);

        vec_sub(out, ac, conj_d_b, half);
        vec_add(out + half, da, b_conj_c, half);
    }
}

/// Dynamic Runtime Multiplication for arbitrary power-of-two dimension `n`
/// Dispatches small dynamic sizes (1, 2, 4, 8, 16) to the fast unrolled/SIMD kernels
/// and uses cache-friendly bitwise XOR-table evaluation for arbitrary higher degrees.
template <typename T>
NION_HOST_DEVICE inline void mul_dynamic(
    T* NION_RESTRICT out, const T* NION_RESTRICT a, const T* NION_RESTRICT b, std::size_t n) noexcept
{
    switch (n) {
        case 0: return;
        case 1: out[0] = a[0] * b[0]; return;
        case 2: mul_fixed_2(out, a, b); return;
        case 4: mul_fixed_4(out, a, b); return;
        case 8: mul_fixed_table<8>(out, a, b); return;
        case 16: mul_fixed_table<16>(out, a, b); return;
        case 32: mul_fixed_table<32>(out, a, b); return;
        case 64: mul_fixed_table<64>(out, a, b); return;
        default: break;
    }

    // General power-of-two dimension > 16: zero-allocation bitwise XOR basis product
    for (std::size_t k = 0; k < n; ++k) {
        T acc{};
        for (std::size_t p = 0; p < n; ++p) {
            const std::size_t q = p ^ k;
            const int s = cd_sign_bitwise(p, q);
            if (s > 0) {
                acc = acc + a[p] * b[q];
            } else {
                acc = acc - a[p] * b[q];
            }
        }
        out[k] = acc;
    }
}

} // namespace nion_math::simd

#endif // NION_SIMD_HPP
