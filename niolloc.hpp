// Copyright 2022-2026 Marcus Dante Liebenthal & Contributors
// Licensed under the Apache License, Version 2.0 (the "License");
//
// niolloc.hpp - Modernized Modular Memory & Storage Subsystem for Nion
// Supports:
//   1. Compile-time fixed-size stack/register storage (GPU kernel by-value safe)
//   2. SIMD-aligned CPU heap allocation (64-byte cache-line / AVX-512 aligned)
//   3. GPU Device & Unified (Managed) memory spaces (CUDA / HIP ready)
//   4. Non-owning views/spans for zero-copy GPU kernel execution

#ifndef NIOLLOC_HPP
#define NIOLLOC_HPP

#include <algorithm>
#include <bit>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

// ============================================================================
// 1. Execution Space & Compiler Portability Macros (CPU / CUDA / HIP)
// ============================================================================
#if defined(__CUDACC__) || defined(__HIPCC__)
    #define NION_GPU_COMPILER 1
    #define NION_HOST_DEVICE __host__ __device__
    #define NION_HOST __host__
    #define NION_DEVICE __device__
#else
    #define NION_GPU_COMPILER 0
    #define NION_HOST_DEVICE
    #define NION_HOST
    #define NION_DEVICE
#endif

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    #define NION_ON_DEVICE 1
#else
    #define NION_ON_DEVICE 0
#endif

#if defined(_MSC_VER)
    #define NION_FORCE_INLINE __forceinline
    #define NION_RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define NION_FORCE_INLINE inline __attribute__((always_inline))
    #define NION_RESTRICT __restrict__
#else
    #define NION_FORCE_INLINE inline
    #define NION_RESTRICT
#endif

// Optional CUDA / HIP runtime inclusion for GPU allocators
#if defined(NION_ENABLE_CUDA) && __has_include(<cuda_runtime.h>)
    #include <cuda_runtime.h>
    #define NION_HAS_CUDA_RUNTIME 1
#else
    #define NION_HAS_CUDA_RUNTIME 0
#endif

#if defined(NION_ENABLE_HIP) && __has_include(<hip/hip_runtime.h>)
    #include <hip/hip_runtime.h>
    #define NION_HAS_HIP_RUNTIME 1
#else
    #define NION_HAS_HIP_RUNTIME 0
#endif

namespace nion_math {

/// Dynamic size sentinel (backward-compatible with original NION_USE_HEAP)
inline constexpr std::size_t NION_USE_HEAP = 0;
inline constexpr std::size_t NION_DYNAMIC  = 0;

/// Default SIMD / Cache-line byte alignment (64B covers AVX-512, AVX2, NEON, GPU coalescing)
inline constexpr std::size_t NION_DEFAULT_ALIGNMENT = 64;

// ============================================================================
// 2. C++20 Algebraic & Memory Concepts
// ============================================================================
namespace concepts {

/// A valid hypercomplex scalar coefficient type T must support ring/field arithmetic
template <typename T>
concept Scalar = requires(T a, T b) {
    { a + b } -> std::convertible_to<T>;
    { a - b } -> std::convertible_to<T>;
    { a * b } -> std::convertible_to<T>;
    { a / b } -> std::convertible_to<T>;
    { -a }    -> std::convertible_to<T>;
} && std::is_default_constructible_v<T> && std::is_copy_constructible_v<T>;

/// Primitive floating-point types eligible for hardware SIMD intrinsics
template <typename T>
concept SimdFloat = std::is_same_v<std::remove_cv_t<T>, float> ||
                    std::is_same_v<std::remove_cv_t<T>, double>;

/// Checks whether a compile-time dimension N is either dynamic (0) or a power of two (2^k)
template <std::size_t N>
concept ValidDimension = (N == NION_DYNAMIC) || std::has_single_bit(N);

} // namespace concepts

// ============================================================================
// 3. Utility Functions for Power-of-Two Hypercomplex Dimensions
// ============================================================================
namespace detail {

/// Tag for internal uninitialized construction when the caller immediately overwrites all elements
struct uninit_tag {
    explicit constexpr uninit_tag() = default;
};

NION_HOST_DEVICE constexpr std::size_t next_power_of_two(std::size_t n) noexcept {
    if (n <= 1) return 1;
    // Bit-twiddling ceiling power of two that works in both constexpr and GPU device code
    --n;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    if constexpr (sizeof(std::size_t) > 4) {
        n |= n >> 32;
    }
    return n + 1;
}

NION_HOST_DEVICE constexpr std::size_t log2_pow2(std::size_t n) noexcept {
    std::size_t deg = 0;
    while (n > 1) {
        n >>= 1;
        ++deg;
    }
    return deg;
}

/// Safe element-wise fill (avoids memset misuse on non-trivial or non-zero floats)
template <typename T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr void fill_n(
    T* NION_RESTRICT dst, std::size_t count, const T& val) noexcept
{
    for (std::size_t i = 0; i < count; ++i) {
        dst[i] = val;
    }
}

/// Safe element-wise copy (constexpr and GPU device compatible)
template <typename T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr void copy_n(
    const T* NION_RESTRICT src, std::size_t count, T* NION_RESTRICT dst) noexcept
{
    for (std::size_t i = 0; i < count; ++i) {
        dst[i] = src[i];
    }
}

} // namespace detail

// ============================================================================
// 4. Modular Memory Space Policies (CPU Aligned, GPU Unified, GPU Device)
// ============================================================================

enum class MemoryLocation {
    Host,
    Device,
    Unified
};

/// Policy 1: Aligned CPU Host Memory (Default)
/// Allocates cache-line / SIMD aligned memory on the CPU heap.
struct HostSpace {
    static constexpr MemoryLocation location = MemoryLocation::Host;
    static constexpr bool host_accessible    = true;
    static constexpr bool device_accessible  = false;

    template <typename T>
    static T* allocate(std::size_t count) {
        if (count == 0) return nullptr;
        const std::size_t bytes = count * sizeof(T);
        constexpr std::size_t align = (alignof(T) > NION_DEFAULT_ALIGNMENT)
                                      ? alignof(T) : NION_DEFAULT_ALIGNMENT;
        void* ptr = ::operator new(bytes, std::align_val_t{align});
        if constexpr (!std::is_trivially_default_constructible_v<T>) {
            T* typed_ptr = static_cast<T*>(ptr);
            for (std::size_t i = 0; i < count; ++i) {
                std::construct_at(typed_ptr + i);
            }
            return typed_ptr;
        }
        return static_cast<T*>(ptr);
    }

    template <typename T>
    static void deallocate(T* ptr, std::size_t count) noexcept {
        if (!ptr) return;
        if constexpr (!std::is_trivially_destructible_v<T>) {
            for (std::size_t i = 0; i < count; ++i) {
                std::destroy_at(ptr + i);
            }
        }
        constexpr std::size_t align = (alignof(T) > NION_DEFAULT_ALIGNMENT)
                                      ? alignof(T) : NION_DEFAULT_ALIGNMENT;
        ::operator delete(static_cast<void*>(ptr), std::align_val_t{align});
    }

    template <typename T>
    static void copy(T* dst, const T* src, std::size_t count) {
        detail::copy_n(src, count, dst);
    }
};

/// Policy 2: Unified (Managed) CPU/GPU Memory
/// Accessible on both CPU and GPU via hardware page migration (cudaMallocManaged / hipMallocManaged).
/// Falls back cleanly to aligned HostSpace when compiled without a GPU runtime.
struct UnifiedSpace {
    static constexpr MemoryLocation location = MemoryLocation::Unified;
    static constexpr bool host_accessible    = true;
    static constexpr bool device_accessible  = true;

    template <typename T>
    static T* allocate(std::size_t count) {
        if (count == 0) return nullptr;
#if NION_HAS_CUDA_RUNTIME
        T* ptr = nullptr;
        cudaError_t err = cudaMallocManaged(reinterpret_cast<void**>(&ptr), count * sizeof(T));
        if (err != cudaSuccess) throw std::bad_alloc();
        return ptr;
#elif NION_HAS_HIP_RUNTIME
        T* ptr = nullptr;
        hipError_t err = hipMallocManaged(reinterpret_cast<void**>(&ptr), count * sizeof(T));
        if (err != hipSuccess) throw std::bad_alloc();
        return ptr;
#else
        return HostSpace::allocate<T>(count);
#endif
    }

    template <typename T>
    static void deallocate(T* ptr, std::size_t count) noexcept {
        if (!ptr) return;
#if NION_HAS_CUDA_RUNTIME
        (void)count;
        cudaFree(ptr);
#elif NION_HAS_HIP_RUNTIME
        (void)count;
        hipFree(ptr);
#else
        HostSpace::deallocate(ptr, count);
#endif
    }

    template <typename T>
    static void copy(T* dst, const T* src, std::size_t count) {
        detail::copy_n(src, count, dst);
    }
};

/// Policy 3: Dedicated High-Bandwidth GPU Device Memory (VRAM)
/// Allocated on the host via cudaMalloc / hipMalloc and passed into GPU kernels via views.
struct DeviceSpace {
    static constexpr MemoryLocation location = MemoryLocation::Device;
    static constexpr bool host_accessible    = false;
    static constexpr bool device_accessible  = true;

    template <typename T>
    static T* allocate(std::size_t count) {
        if (count == 0) return nullptr;
        static_assert(std::is_trivially_copyable_v<T>,
                      "DeviceSpace requires trivially copyable scalar types T.");
#if NION_HAS_CUDA_RUNTIME
        T* ptr = nullptr;
        cudaError_t err = cudaMalloc(reinterpret_cast<void**>(&ptr), count * sizeof(T));
        if (err != cudaSuccess) throw std::bad_alloc();
        return ptr;
#elif NION_HAS_HIP_RUNTIME
        T* ptr = nullptr;
        hipError_t err = hipMalloc(reinterpret_cast<void**>(&ptr), count * sizeof(T));
        if (err != hipSuccess) throw std::bad_alloc();
        return ptr;
#else
        // Emulated device memory on CPU-only builds for seamless testing
        return HostSpace::allocate<T>(count);
#endif
    }

    template <typename T>
    static void deallocate(T* ptr, std::size_t count) noexcept {
        if (!ptr) return;
#if NION_HAS_CUDA_RUNTIME
        (void)count;
        cudaFree(ptr);
#elif NION_HAS_HIP_RUNTIME
        (void)count;
        hipFree(ptr);
#else
        HostSpace::deallocate(ptr, count);
#endif
    }

    template <typename T>
    static void copy(T* dst, const T* src, std::size_t count) {
        if (count == 0) return;
#if NION_HAS_CUDA_RUNTIME
        cudaMemcpy(dst, src, count * sizeof(T), cudaMemcpyDefault);
#elif NION_HAS_HIP_RUNTIME
        hipMemcpy(dst, src, count * sizeof(T), hipMemcpyDefault);
#else
        detail::copy_n(src, count, dst);
#endif
    }

    template <typename T>
    static void copy_from_host(T* dst_dev, const T* src_host, std::size_t count) {
        if (count == 0) return;
#if NION_HAS_CUDA_RUNTIME
        cudaMemcpy(dst_dev, src_host, count * sizeof(T), cudaMemcpyHostToDevice);
#elif NION_HAS_HIP_RUNTIME
        hipMemcpy(dst_dev, src_host, count * sizeof(T), hipMemcpyHostToDevice);
#else
        detail::copy_n(src_host, count, dst_dev);
#endif
    }

    template <typename T>
    static void copy_to_host(T* dst_host, const T* src_dev, std::size_t count) {
        if (count == 0) return;
#if NION_HAS_CUDA_RUNTIME
        cudaMemcpy(dst_host, src_dev, count * sizeof(T), cudaMemcpyDeviceToHost);
#elif NION_HAS_HIP_RUNTIME
        hipMemcpy(dst_host, src_dev, count * sizeof(T), hipMemcpyDeviceToHost);
#else
        detail::copy_n(src_dev, count, dst_host);
#endif
    }
};

// ============================================================================
// 5. niolloc_storage: Unified Storage Backend (Fixed Stack vs. Dynamic Space)
// ============================================================================

/// Primary Template: Compile-Time Fixed Storage (N > 0)
/// Zero heap allocations, SIMD aligned, trivially copyable into GPU kernels!
template <concepts::Scalar T, std::size_t N, typename Space = HostSpace>
    requires concepts::ValidDimension<N>
class niolloc_storage {
public:
    static constexpr std::size_t static_size = N;
    static constexpr bool is_dynamic         = false;
    static constexpr std::size_t alignment   =
        concepts::SimdFloat<T> ? (
            (N * sizeof(T) >= 64) ? 64 :
            (N * sizeof(T) >= 32) ? 32 :
            (N * sizeof(T) >= 16) ? 16 : alignof(T)
        ) : alignof(T);

    alignas(alignment) T m_data[N];

    NION_HOST_DEVICE constexpr niolloc_storage() noexcept : m_data{} {}

    NION_HOST_DEVICE constexpr niolloc_storage(detail::uninit_tag, std::size_t count) noexcept {
        assert(count <= N && "Fixed-size niolloc_storage cannot exceed N");
        (void)count;
        if (std::is_constant_evaluated()) {
            detail::fill_n(m_data, N, T{});
        }
    }

    NION_HOST_DEVICE constexpr explicit niolloc_storage(std::size_t count) noexcept : m_data{} {
        assert(count <= N && "Fixed-size niolloc_storage cannot exceed N");
        (void)count;
    }

    NION_HOST_DEVICE constexpr niolloc_storage(std::size_t count, const T& val) noexcept {
        assert(count <= N && "Fixed-size niolloc_storage cannot exceed N");
        (void)count;
        detail::fill_n(m_data, N, val);
    }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T* data() noexcept { return m_data; }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr const T* data() const noexcept { return m_data; }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr std::size_t size() const noexcept { return N; }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T& operator[](std::size_t i) noexcept {
        return m_data[i];
    }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr const T& operator[](std::size_t i) const noexcept {
        return m_data[i];
    }

    NION_HOST_DEVICE constexpr void resize(std::size_t new_size) noexcept {
        assert(new_size == N && "Cannot resize fixed-dimension niolloc_storage");
        (void)new_size;
    }

    NION_HOST_DEVICE constexpr void fill(const T& val) noexcept {
        detail::fill_n(m_data, N, val);
    }
};

/// Partial Specialization: Dynamic Runtime Storage (N == NION_USE_HEAP / 0)
/// Uses RAII with the specified Memory Space Policy (HostSpace, UnifiedSpace, or DeviceSpace).
template <concepts::Scalar T, typename Space>
class niolloc_storage<T, NION_DYNAMIC, Space> {
private:
    T* m_data          = nullptr;
    std::size_t m_size = 0;

public:
    static constexpr std::size_t static_size = NION_DYNAMIC;
    static constexpr bool is_dynamic         = true;
    using space_type                         = Space;

    NION_HOST_DEVICE constexpr niolloc_storage() noexcept = default;

    niolloc_storage(detail::uninit_tag, std::size_t count)
        : m_size(detail::next_power_of_two(count))
    {
        m_data = Space::template allocate<T>(m_size);
    }

    explicit niolloc_storage(std::size_t count)
        : m_size(detail::next_power_of_two(count))
    {
        m_data = Space::template allocate<T>(m_size);
        if constexpr (Space::host_accessible) {
            detail::fill_n(m_data, m_size, T{});
        }
    }

    niolloc_storage(std::size_t count, const T& val)
        : m_size(detail::next_power_of_two(count))
    {
        m_data = Space::template allocate<T>(m_size);
        if constexpr (Space::host_accessible) {
            detail::fill_n(m_data, m_size, val);
        }
    }

    // Copy constructor
    niolloc_storage(const niolloc_storage& other)
        : m_size(other.m_size)
    {
        if (m_size > 0) {
            m_data = Space::template allocate<T>(m_size);
            Space::copy(m_data, other.m_data, m_size);
        }
    }

    // Move constructor
    NION_HOST_DEVICE niolloc_storage(niolloc_storage&& other) noexcept
        : m_data(std::exchange(other.m_data, nullptr)),
          m_size(std::exchange(other.m_size, 0))
    {}

    // Copy assignment
    niolloc_storage& operator=(const niolloc_storage& other) {
        if (this != &other) {
            if (m_size != other.m_size) {
                Space::deallocate(m_data, m_size);
                m_size = other.m_size;
                m_data = (m_size > 0) ? Space::template allocate<T>(m_size) : nullptr;
            }
            if (m_size > 0) {
                Space::copy(m_data, other.m_data, m_size);
            }
        }
        return *this;
    }

    // Move assignment
    NION_HOST_DEVICE niolloc_storage& operator=(niolloc_storage&& other) noexcept {
        if (this != &other) {
#if !NION_ON_DEVICE
            Space::deallocate(m_data, m_size);
#endif
            m_data = std::exchange(other.m_data, nullptr);
            m_size = std::exchange(other.m_size, 0);
        }
        return *this;
    }

    ~niolloc_storage() noexcept {
#if !NION_ON_DEVICE
        Space::deallocate(m_data, m_size);
#endif
    }

    NION_HOST_DEVICE NION_FORCE_INLINE T* data() noexcept { return m_data; }
    NION_HOST_DEVICE NION_FORCE_INLINE const T* data() const noexcept { return m_data; }
    NION_HOST_DEVICE NION_FORCE_INLINE std::size_t size() const noexcept { return m_size; }

    NION_HOST_DEVICE NION_FORCE_INLINE T& operator[](std::size_t i) noexcept {
        return m_data[i];
    }
    NION_HOST_DEVICE NION_FORCE_INLINE const T& operator[](std::size_t i) const noexcept {
        return m_data[i];
    }

    void resize(std::size_t new_count) {
        const std::size_t target = detail::next_power_of_two(new_count);
        if (target == m_size) return;
        T* new_data = Space::template allocate<T>(target);
        if constexpr (Space::host_accessible) {
            const std::size_t copy_len = (m_size < target) ? m_size : target;
            if (m_data && copy_len > 0) {
                detail::copy_n(m_data, copy_len, new_data);
            }
            if (target > copy_len) {
                detail::fill_n(new_data + copy_len, target - copy_len, T{});
            }
        } else if (m_data && m_size > 0) {
            Space::copy(new_data, m_data, (m_size < target) ? m_size : target);
        }
        Space::deallocate(m_data, m_size);
        m_data = new_data;
        m_size = target;
    }

    void fill(const T& val) {
        static_assert(Space::host_accessible,
                      "Direct fill() requires host-accessible memory (HostSpace or UnifiedSpace).");
        detail::fill_n(m_data, m_size, val);
    }
};

} // namespace nion_math

namespace nions = nion_math;

// Bring NION_USE_HEAP and memory spaces into global scope for seamless compatibility
using nion_math::NION_USE_HEAP;
using nion_math::NION_DYNAMIC;
using nion_math::HostSpace;
using nion_math::UnifiedSpace;
using nion_math::DeviceSpace;

#endif // NIOLLOC_HPP
