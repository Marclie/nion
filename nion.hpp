// Copyright 2022-2026 Marcus Dante Liebenthal & Contributors
// Licensed under the Apache License, Version 2.0 (the "License");
//
// nion.hpp - Core Hypercomplex Class, Expression Templates, Views, and SoA Batching
// Features:
//   - Backward-compatible API (`nion<float> a = {1, 2, 3, 4};`, `a.conj()`, `a.inv()`, `++a`, `--a`)
//   - Compile-Time Fixed (`nion<T, 4>`) and Dynamic (`nion<T, NION_USE_HEAP, Space>`) dimensions
//   - Zero-Overhead Expression Templates for fused addition/subtraction/scaling/conjugation
//   - Non-owning `nion_view<T, N>` for zero-copy GPU kernel execution
//   - `nion_batch_soa<T, N, Space>` for coalesced GPU and multi-lane SIMD batch processing

#ifndef NION_HPP
#define NION_HPP

#include "niolloc.hpp"
#include "nion_simd.hpp"

#include <cmath>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <vector>

namespace nion_math {

// ============================================================================
// 1. Expression Template Infrastructure (Zero-Temporary Linear Combinations)
// ============================================================================
// Any linear expression of the form `a + b - c * scalar + d.conj()` is fused
// into a single loop at assignment time, completely eliminating temporary heap
// or stack allocations!

template <typename Derived, typename T>
struct NionExpr {
    using is_nion_expr_tag = void;

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr const Derived& derived() const noexcept {
        return static_cast<const Derived&>(*this);
    }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr std::size_t size() const noexcept {
        return derived().size();
    }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr bool has_uniform_size(std::size_t n) const noexcept {
        return derived().has_uniform_size(n);
    }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T eval(std::size_t i) const noexcept {
        return derived().eval(i);
    }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T eval_unchecked(std::size_t i) const noexcept {
        return derived().eval_unchecked(i);
    }
    NION_HOST_DEVICE constexpr T abs() const noexcept {
        const Derived& d = derived();
        const std::size_t n = d.size();
        if (n == 0) return T{};
        if (d.has_uniform_size(n)) {
            T sum = d.eval_unchecked(0);
            sum = sum * sum;
            for (std::size_t i = 1; i < n; ++i) {
                const T v = d.eval_unchecked(i);
                sum = sum + v * v;
            }
            return sum;
        } else {
            T sum = d.eval(0);
            sum = sum * sum;
            for (std::size_t i = 1; i < n; ++i) {
                const T v = d.eval(i);
                sum = sum + v * v;
            }
            return sum;
        }
    }
    NION_HOST_DEVICE T norm() const noexcept {
        using std::sqrt;
        return static_cast<T>(sqrt(abs()));
    }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T real() const noexcept {
        return derived().eval(0);
    }
};

template <typename S, typename T>
concept ConvertibleScalar =
    std::convertible_to<S, T> &&
    (!requires { typename std::remove_cvref_t<S>::is_nion_expr_tag; });

template <typename Lhs, typename Rhs, typename T>
class NionAddExpr : public NionExpr<NionAddExpr<Lhs, Rhs, T>, T> {
    const Lhs& m_lhs;
    const Rhs& m_rhs;
    std::size_t m_size;
public:
    static constexpr std::size_t static_size =
        (Lhs::static_size != NION_DYNAMIC && Lhs::static_size == Rhs::static_size)
            ? Lhs::static_size : NION_DYNAMIC;

    NION_HOST_DEVICE constexpr NionAddExpr(const Lhs& l, const Rhs& r) noexcept
        : m_lhs(l), m_rhs(r), m_size(l.size() > r.size() ? l.size() : r.size()) {}

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr std::size_t size() const noexcept {
        if constexpr (static_size != NION_DYNAMIC) return static_size;
        else return m_size;
    }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr bool has_uniform_size(std::size_t n) const noexcept {
        return m_lhs.has_uniform_size(n) && m_rhs.has_uniform_size(n);
    }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T eval(std::size_t i) const noexcept {
        if constexpr (static_size != NION_DYNAMIC) {
            return m_lhs.eval_unchecked(i) + m_rhs.eval_unchecked(i);
        } else {
            const T lv = (i < m_lhs.size()) ? m_lhs.eval(i) : T{};
            const T rv = (i < m_rhs.size()) ? m_rhs.eval(i) : T{};
            return lv + rv;
        }
    }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T eval_unchecked(std::size_t i) const noexcept {
        return m_lhs.eval_unchecked(i) + m_rhs.eval_unchecked(i);
    }
};

template <typename Lhs, typename Rhs, typename T>
class NionSubExpr : public NionExpr<NionSubExpr<Lhs, Rhs, T>, T> {
    const Lhs& m_lhs;
    const Rhs& m_rhs;
    std::size_t m_size;
public:
    static constexpr std::size_t static_size =
        (Lhs::static_size != NION_DYNAMIC && Lhs::static_size == Rhs::static_size)
            ? Lhs::static_size : NION_DYNAMIC;

    NION_HOST_DEVICE constexpr NionSubExpr(const Lhs& l, const Rhs& r) noexcept
        : m_lhs(l), m_rhs(r), m_size(l.size() > r.size() ? l.size() : r.size()) {}

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr std::size_t size() const noexcept {
        if constexpr (static_size != NION_DYNAMIC) return static_size;
        else return m_size;
    }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr bool has_uniform_size(std::size_t n) const noexcept {
        return m_lhs.has_uniform_size(n) && m_rhs.has_uniform_size(n);
    }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T eval(std::size_t i) const noexcept {
        if constexpr (static_size != NION_DYNAMIC) {
            return m_lhs.eval_unchecked(i) - m_rhs.eval_unchecked(i);
        } else {
            const T lv = (i < m_lhs.size()) ? m_lhs.eval(i) : T{};
            const T rv = (i < m_rhs.size()) ? m_rhs.eval(i) : T{};
            return lv - rv;
        }
    }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T eval_unchecked(std::size_t i) const noexcept {
        return m_lhs.eval_unchecked(i) - m_rhs.eval_unchecked(i);
    }
};

template <typename Expr, typename T>
class NionScaleExpr : public NionExpr<NionScaleExpr<Expr, T>, T> {
    const Expr& m_expr;
    T m_scalar;
public:
    static constexpr std::size_t static_size = Expr::static_size;

    NION_HOST_DEVICE constexpr NionScaleExpr(const Expr& e, const T& s) noexcept
        : m_expr(e), m_scalar(s) {}

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr std::size_t size() const noexcept { return m_expr.size(); }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr bool has_uniform_size(std::size_t n) const noexcept {
        return m_expr.has_uniform_size(n);
    }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T eval(std::size_t i) const noexcept {
        return m_expr.eval(i) * m_scalar;
    }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T eval_unchecked(std::size_t i) const noexcept {
        return m_expr.eval_unchecked(i) * m_scalar;
    }
};

template <typename Expr, typename T>
class NionConjExpr : public NionExpr<NionConjExpr<Expr, T>, T> {
    const Expr& m_expr;
public:
    static constexpr std::size_t static_size = Expr::static_size;

    NION_HOST_DEVICE constexpr explicit NionConjExpr(const Expr& e) noexcept : m_expr(e) {}

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr std::size_t size() const noexcept { return m_expr.size(); }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr bool has_uniform_size(std::size_t n) const noexcept {
        return m_expr.has_uniform_size(n);
    }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T eval(std::size_t i) const noexcept {
        return (i == 0) ? m_expr.eval(0) : -m_expr.eval(i);
    }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T eval_unchecked(std::size_t i) const noexcept {
        return (i == 0) ? m_expr.eval_unchecked(0) : -m_expr.eval_unchecked(i);
    }
};

template <typename Lhs, typename Rhs, typename T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr auto operator+(
    const NionExpr<Lhs, T>& l, const NionExpr<Rhs, T>& r) noexcept
{
    return NionAddExpr<Lhs, Rhs, T>(l.derived(), r.derived());
}

template <typename Lhs, typename Rhs, typename T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr auto operator-(
    const NionExpr<Lhs, T>& l, const NionExpr<Rhs, T>& r) noexcept
{
    return NionSubExpr<Lhs, Rhs, T>(l.derived(), r.derived());
}

template <typename Expr, typename T, ConvertibleScalar<T> S>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr auto operator*(
    const NionExpr<Expr, T>& e, const S& s) noexcept
{
    return NionScaleExpr<Expr, T>(e.derived(), static_cast<T>(s));
}

template <typename Expr, typename T, ConvertibleScalar<T> S>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr auto operator*(
    const S& s, const NionExpr<Expr, T>& e) noexcept
{
    return NionScaleExpr<Expr, T>(e.derived(), static_cast<T>(s));
}

template <typename Expr, typename T, ConvertibleScalar<T> S>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr auto operator/(
    const NionExpr<Expr, T>& e, const S& s) noexcept
{
    if constexpr (std::is_floating_point_v<T>) {
        return NionScaleExpr<Expr, T>(e.derived(), T{1} / static_cast<T>(s));
    } else {
        return NionScaleExpr<Expr, T>(e.derived(), static_cast<T>(T{1} / static_cast<T>(s)));
    }
}

} // namespace nion_math

// ============================================================================
// 2. Core `nion<T, N, Space>` Class Template
// ============================================================================
// Defined at global scope (matching the original Nion library API) while
// internal utilities live in `namespace nion_math`, completely preventing any
// class-vs-namespace name collision!

template <
    nion_math::concepts::Scalar T,
    std::size_t N = nion_math::NION_USE_HEAP,
    typename Space = nion_math::HostSpace
>
    requires nion_math::concepts::ValidDimension<N>
class nion : public nion_math::NionExpr<nion<T, N, Space>, T> {
private:
    nion_math::niolloc_storage<T, N, Space> m_store;

public:
    using value_type                         = T;
    using space_type                         = Space;
    static constexpr std::size_t static_size = N;
    static constexpr bool is_dynamic         = (N == nion_math::NION_DYNAMIC);

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    /// Default constructor: creates a 1-D real scalar (if dynamic) or N-D zero nion (if fixed)
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr nion() : m_store(is_dynamic ? 1 : N, T{}) {}

    /// Internal uninitialized constructor (for operations that immediately overwrite all elements)
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr nion(nion_math::detail::uninit_tag tag, std::size_t dim)
        : m_store(tag, is_dynamic ? dim : N)
    {}

    /// Construct with explicit integral dimension `dim` (zero-initialized), matching master API: `nion<double> n(2);`
    template <std::integral Z>
        requires (!std::is_integral_v<T>)
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr explicit nion(Z dim)
        : m_store(is_dynamic ? static_cast<std::size_t>(dim) : N, T{})
    {}

    /// Construct from scalar value `re` (e_0 = re, all imaginary components = 0)
    template <typename S = T>
        requires (nion_math::ConvertibleScalar<S, T> && (!std::is_integral_v<S> || std::is_integral_v<T>))
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr explicit nion(const S& re)
        : m_store(is_dynamic ? 1 : N, T{})
    {
        if constexpr (Space::host_accessible || !is_dynamic) {
            m_store[0] = static_cast<T>(re);
        }
    }

    /// Construct dynamic or fixed nion with explicit dimension `dim` and initial fill `val`
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr nion(std::size_t dim, const T& val)
        : m_store(dim, val)
    {}

    /// Initializer-list constructor: `nion<float> a = {1, 2, 3, 4};`
    /// Automatically pads to the next power of two if needed.
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr nion(std::initializer_list<T> init)
        : m_store(nion_math::detail::uninit_tag{},
                  is_dynamic ? nion_math::detail::next_power_of_two(init.size()) : N)
    {
        const std::size_t sz = size();
        const std::size_t limit = (init.size() < sz) ? init.size() : sz;
        nion_math::detail::copy_n(init.begin(), limit, m_store.data());
        for (std::size_t i = limit; i < sz; ++i) {
            m_store[i] = T{};
        }
    }

    /// Construct from std::vector<T>
    explicit nion(const std::vector<T>& vec)
        : m_store(nion_math::detail::uninit_tag{},
                  is_dynamic ? nion_math::detail::next_power_of_two(vec.size()) : N)
    {
        const std::size_t sz = size();
        const std::size_t limit = (vec.size() < sz) ? vec.size() : sz;
        for (std::size_t i = 0; i < limit; ++i) {
            m_store[i] = vec[i];
        }
        for (std::size_t i = limit; i < sz; ++i) {
            m_store[i] = T{};
        }
    }

    /// Construct from raw pointer and count
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr nion(const T* ptr, std::size_t count)
        : m_store(nion_math::detail::uninit_tag{},
                  is_dynamic ? nion_math::detail::next_power_of_two(count) : N)
    {
#if !NION_ON_DEVICE && (defined(__GNUC__) || defined(__clang__))
        if (!std::is_constant_evaluated()) {
            asm("" : "+r"(ptr));
        }
#endif
        const std::size_t sz = size();
        const std::size_t limit = (count < sz) ? count : sz;
        nion_math::detail::copy_n(ptr, limit, m_store.data());
        for (std::size_t i = limit; i < sz; ++i) {
            m_store[i] = T{};
        }
    }

    /// Construct & evaluate from any Expression Template node (Single Fused Loop!)
    template <typename Expr>
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr nion(const nion_math::NionExpr<Expr, T>& expr)
        : m_store(nion_math::detail::uninit_tag{}, is_dynamic ? expr.size() : N)
    {
        const Expr& d = expr.derived();
        T* NION_RESTRICT dst = m_store.data();
        if constexpr (!is_dynamic && Expr::static_size == N) {
            for (std::size_t i = 0; i < N; ++i) {
                dst[i] = d.eval_unchecked(i);
            }
        } else {
            const std::size_t len = size();
            if (d.has_uniform_size(len)) {
                for (std::size_t i = 0; i < len; ++i) {
                    dst[i] = d.eval_unchecked(i);
                }
            } else {
                const std::size_t d_sz = d.size();
                for (std::size_t i = 0; i < len; ++i) {
                    dst[i] = (i < d_sz) ? d.eval(i) : T{};
                }
            }
        }
    }

    template <typename Expr>
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr nion& operator=(const nion_math::NionExpr<Expr, T>& expr) {
        const auto& d = expr.derived();
        if constexpr (!is_dynamic && Expr::static_size == N) {
            T* NION_RESTRICT dst = m_store.data();
            for (std::size_t i = 0; i < N; ++i) {
                dst[i] = d.eval_unchecked(i);
            }
        } else {
            if constexpr (is_dynamic) {
                if (size() != d.size()) {
                    m_store.resize(d.size());
                }
            }
            const std::size_t len = size();
            T* NION_RESTRICT dst = m_store.data();
            if (d.has_uniform_size(len)) {
                for (std::size_t i = 0; i < len; ++i) {
                    dst[i] = d.eval_unchecked(i);
                }
            } else {
                const std::size_t d_sz = d.size();
                for (std::size_t i = 0; i < len; ++i) {
                    dst[i] = (i < d_sz) ? d.eval(i) : T{};
                }
            }
        }
        return *this;
    }

    // Standard copy/move constructors & assignments
    // NOTE: no NION_HOST_DEVICE here. NVCC ignores execution-space annotations on functions
    // explicitly defaulted on their first declaration (warning #20012-D) and instead infers
    // __host__ __device__ from the (already NION_HOST_DEVICE) niolloc_storage members.
    constexpr nion(const nion&) = default;
    constexpr nion(nion&&) noexcept = default;
    constexpr nion& operator=(const nion&) = default;
    constexpr nion& operator=(nion&&) noexcept = default;

    /// Cross-space copy helper (e.g., copying HostSpace nion <-> DeviceSpace / UnifiedSpace nion)
    template <std::size_t OtherN, typename OtherSpace>
    explicit nion(const nion<T, OtherN, OtherSpace>& other)
        : m_store(is_dynamic ? other.size() : N)
    {
        const std::size_t copy_len = (size() < other.size()) ? size() : other.size();
        if constexpr (Space::location == nion_math::MemoryLocation::Device &&
                      OtherSpace::location == nion_math::MemoryLocation::Host) {
            nion_math::DeviceSpace::copy_from_host(m_store.data(), other.data(), copy_len);
        } else if constexpr (Space::location == nion_math::MemoryLocation::Host &&
                             OtherSpace::location == nion_math::MemoryLocation::Device) {
            nion_math::DeviceSpace::copy_to_host(m_store.data(), other.data(), copy_len);
        } else {
            Space::copy(m_store.data(), other.data(), copy_len);
        }
    }

    // ------------------------------------------------------------------------
    // Accessors & Dimension Queries
    // ------------------------------------------------------------------------

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T* data() noexcept { return m_store.data(); }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr const T* data() const noexcept { return m_store.data(); }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr std::size_t size() const noexcept {
        return m_store.size();
    }

    /// Returns the Cayley-Dickson algebra degree D where size == 2^D
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr std::size_t degree() const noexcept {
        return nion_math::detail::log2_pow2(size());
    }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T& operator[](std::size_t i) noexcept {
        return m_store[i];
    }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr const T& operator[](std::size_t i) const noexcept {
        return m_store[i];
    }

    /// Expression Template evaluation hooks
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr bool has_uniform_size(std::size_t n) const noexcept {
        return size() == n;
    }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T eval(std::size_t i) const noexcept {
        return m_store[i];
    }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T eval_unchecked(std::size_t i) const noexcept {
        return m_store[i];
    }

    /// Real part (scalar coefficient e_0)
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T real() const noexcept {
        return m_store[0];
    }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T& real() noexcept {
        return m_store[0];
    }

    /// Pure imaginary part (returns a nion with e_0 = 0 and e_1...e_{n-1} preserved)
    NION_HOST_DEVICE constexpr nion imag() const {
        nion res(*this);
        res[0] = T{};
        return res;
    }

    /// Resize dynamic nion to the next power-of-two >= new_dim and return *this
    nion& resize(std::size_t new_dim) {
        m_store.resize(new_dim);
        return *this;
    }

    /// Checks whether all imaginary components (e_1 ... e_{n-1}) are zero
    NION_HOST_DEVICE constexpr bool is_real(const T& tol = T{}) const noexcept {
        for (std::size_t i = 1; i < size(); ++i) {
            const T v = (m_store[i] < T{}) ? -m_store[i] : m_store[i];
            if (v > tol) return false;
        }
        return true;
    }

    // ------------------------------------------------------------------------
    // Core Algebraic Member Functions (conj, norm, abs, inv)
    // ------------------------------------------------------------------------

    /// Returns a materialized conjugate nion
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr nion conj() const {
        nion res(nion_math::detail::uninit_tag{}, size());
        nion_math::simd::vec_conj(res.data(), data(), size());
        return res;
    }

    /// Lazy conjugate expression (zero temporaries when chained in expressions)
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr auto conj_expr() const noexcept {
        return nion_math::NionConjExpr<nion, T>(*this);
    }

    /// Squared norm: ||q||^2 = sum_{i=0}^{n-1} q_i^2 (SIMD-accelerated)
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T norm_sq() const noexcept {
        return nion_math::simd::vec_norm_sq(data(), size());
    }

    /// Squared magnitude / Cayley-Dickson norm: sum_{i=0}^{n-1} q_i^2
    /// Matches README.md, original Nion API, std::norm(std::complex), and Boost/Harder norm()
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T abs() const noexcept {
        return norm_sq();
    }

    /// Euclidean magnitude: ||q|| = sqrt(sum_{i=0}^{n-1} q_i^2)
    /// Matches README.md and original Nion API
    NION_HOST_DEVICE NION_FORCE_INLINE T norm() const noexcept {
        using std::sqrt;
        return static_cast<T>(sqrt(norm_sq()));
    }

    /// Multiplicative inverse: q^{-1} = conj(q) / abs(q)
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr nion inv() const {
        const T n_sq = abs();
        const std::size_t sz = size();
        nion res(nion_math::detail::uninit_tag{}, sz);
        if constexpr (std::is_floating_point_v<T>) {
            const T inv_norm = T{1} / n_sq;
            res[0] = m_store[0] * inv_norm;
            for (std::size_t i = 1; i < sz; ++i) {
                res[i] = -m_store[i] * inv_norm;
            }
        } else {
            res[0] = m_store[0] / n_sq;
            for (std::size_t i = 1; i < sz; ++i) {
                res[i] = (-m_store[i]) / n_sq;
            }
        }
        return res;
    }

    // ------------------------------------------------------------------------
    // Unary & Increment/Decrement Operators (Matches Original Nion Semantics)
    // ------------------------------------------------------------------------

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr nion operator+() const { return *this; }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr nion operator-() const {
        const std::size_t sz = size();
        nion res(nion_math::detail::uninit_tag{}, sz);
        for (std::size_t i = 0; i < sz; ++i) {
            res[i] = -m_store[i];
        }
        return res;
    }

    /// Pre-increment: increments real component e_0 (`++a`)
    NION_HOST_DEVICE constexpr nion& operator++() noexcept {
        m_store[0] = m_store[0] + T{1};
        return *this;
    }

    /// Post-increment: increments real component e_0 (`a++`)
    NION_HOST_DEVICE constexpr nion operator++(int) {
        nion tmp(*this);
        m_store[0] = m_store[0] + T{1};
        return tmp;
    }

    /// Pre-decrement: decrements real component e_0 (`--a`)
    NION_HOST_DEVICE constexpr nion& operator--() noexcept {
        m_store[0] = m_store[0] - T{1};
        return *this;
    }

    /// Post-decrement: decrements real component e_0 (`a--`)
    NION_HOST_DEVICE constexpr nion operator--(int) {
        nion tmp(*this);
        m_store[0] = m_store[0] - T{1};
        return tmp;
    }

    // ------------------------------------------------------------------------
    // Compound Assignment Operators (SIMD-Accelerated)
    // ------------------------------------------------------------------------

    NION_HOST_DEVICE constexpr nion& operator+=(const nion& rhs) {
        if constexpr (is_dynamic) {
            if (rhs.size() > size()) m_store.resize(rhs.size());
        }
        if (size() == rhs.size()) {
            nion_math::simd::vec_add(data(), data(), rhs.data(), size());
        } else {
            const std::size_t min_sz = (size() < rhs.size()) ? size() : rhs.size();
            nion_math::simd::vec_add(data(), data(), rhs.data(), min_sz);
        }
        return *this;
    }

    NION_HOST_DEVICE constexpr nion& operator-=(const nion& rhs) {
        if constexpr (is_dynamic) {
            if (rhs.size() > size()) m_store.resize(rhs.size());
        }
        if (size() == rhs.size()) {
            nion_math::simd::vec_sub(data(), data(), rhs.data(), size());
        } else {
            const std::size_t min_sz = (size() < rhs.size()) ? size() : rhs.size();
            nion_math::simd::vec_sub(data(), data(), rhs.data(), min_sz);
        }
        return *this;
    }

    NION_HOST_DEVICE constexpr nion& operator+=(const T& scalar) noexcept {
        m_store[0] = m_store[0] + scalar;
        return *this;
    }

    NION_HOST_DEVICE constexpr nion& operator-=(const T& scalar) noexcept {
        m_store[0] = m_store[0] - scalar;
        return *this;
    }

    NION_HOST_DEVICE constexpr nion& operator*=(const T& scalar) noexcept {
        nion_math::simd::vec_scale(data(), data(), scalar, size());
        return *this;
    }

    NION_HOST_DEVICE constexpr nion& operator/=(const T& scalar) noexcept {
        nion_math::simd::vec_scale(data(), data(), T{1} / scalar, size());
        return *this;
    }

    NION_HOST_DEVICE constexpr nion& operator*=(const nion& rhs) {
        *this = (*this) * rhs;
        return *this;
    }

    NION_HOST_DEVICE constexpr nion& operator/=(const nion& rhs) {
        *this = (*this) / rhs;
        return *this;
    }

    // ------------------------------------------------------------------------
    // Comparison Operators
    // ------------------------------------------------------------------------

    NION_HOST_DEVICE constexpr bool operator==(const nion& rhs) const noexcept {
        const std::size_t max_sz = (size() > rhs.size()) ? size() : rhs.size();
        for (std::size_t i = 0; i < max_sz; ++i) {
            const T lv = (i < size()) ? m_store[i] : T{};
            const T rv = (i < rhs.size()) ? rhs[i] : T{};
            if (lv != rv) return false;
        }
        return true;
    }

    NION_HOST_DEVICE constexpr bool operator!=(const nion& rhs) const noexcept {
        return !(*this == rhs);
    }

    /// Stream formatting: `(a0,a1,a2,a3)` matching original Nion output
    friend std::ostream& operator<<(std::ostream& os, const nion& q) {
        os << '(';
        for (std::size_t i = 0; i < q.size(); ++i) {
            os << q[i];
            if (i + 1 < q.size()) os << ',';
        }
        os << ')';
        return os;
    }
};

// ============================================================================
// 3. Binary Arithmetic Operators (Scalar & Hypercomplex Multiplication/Division)
// ============================================================================

/// nion + scalar
template <nion_math::concepts::Scalar T, std::size_t N, typename Space, nion_math::ConvertibleScalar<T> S>
NION_HOST_DEVICE constexpr nion<T, N, Space> operator+(
    const nion<T, N, Space>& q, const S& s)
{
    nion<T, N, Space> res(q);
    res += static_cast<T>(s);
    return res;
}

/// scalar + nion
template <nion_math::concepts::Scalar T, std::size_t N, typename Space, nion_math::ConvertibleScalar<T> S>
NION_HOST_DEVICE constexpr nion<T, N, Space> operator+(
    const S& s, const nion<T, N, Space>& q)
{
    nion<T, N, Space> res(q);
    res += static_cast<T>(s);
    return res;
}

/// nion - scalar
template <nion_math::concepts::Scalar T, std::size_t N, typename Space, nion_math::ConvertibleScalar<T> S>
NION_HOST_DEVICE constexpr nion<T, N, Space> operator-(
    const nion<T, N, Space>& q, const S& s)
{
    nion<T, N, Space> res(q);
    res -= static_cast<T>(s);
    return res;
}

/// scalar - nion
template <nion_math::concepts::Scalar T, std::size_t N, typename Space, nion_math::ConvertibleScalar<T> S>
NION_HOST_DEVICE constexpr nion<T, N, Space> operator-(
    const S& s, const nion<T, N, Space>& q)
{
    nion<T, N, Space> res = -q;
    res += static_cast<T>(s);
    return res;
}

/// Hypercomplex Multiplication: `nion * nion`
/// Dispatches at compile time to SIMD / unrolled table kernels when `N > 0`,
/// and to dynamic SIMD/table kernels when `N == NION_USE_HEAP`.
template <nion_math::concepts::Scalar T, std::size_t N1, std::size_t N2, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr auto operator*(
    const nion<T, N1, Space>& lhs, const nion<T, N2, Space>& rhs)
{
    if constexpr (N1 != nion_math::NION_DYNAMIC && N2 != nion_math::NION_DYNAMIC) {
        constexpr std::size_t OutN = (N1 > N2) ? N1 : N2;
        nion<T, OutN, Space> res(nion_math::detail::uninit_tag{}, OutN);
        if constexpr (N1 == N2) {
            nion_math::simd::mul_fixed_table<OutN>(res.data(), lhs.data(), rhs.data());
        } else {
            nion<T, OutN, Space> l_pad(lhs.data(), N1);
            nion<T, OutN, Space> r_pad(rhs.data(), N2);
            nion_math::simd::mul_fixed_table<OutN>(res.data(), l_pad.data(), r_pad.data());
        }
        return res;
    } else {
        const std::size_t out_sz = (lhs.size() > rhs.size()) ? lhs.size() : rhs.size();
        nion<T, nion_math::NION_DYNAMIC, Space> res(nion_math::detail::uninit_tag{}, out_sz);
        if (lhs.size() == rhs.size()) {
            nion_math::simd::mul_dynamic(res.data(), lhs.data(), rhs.data(), out_sz);
        } else {
            nion<T, nion_math::NION_DYNAMIC, Space> l_pad(out_sz, T{});
            nion_math::detail::copy_n(lhs.data(), lhs.size(), l_pad.data());
            nion<T, nion_math::NION_DYNAMIC, Space> r_pad(out_sz, T{});
            nion_math::detail::copy_n(rhs.data(), rhs.size(), r_pad.data());
            nion_math::simd::mul_dynamic(res.data(), l_pad.data(), r_pad.data(), out_sz);
        }
        return res;
    }
}

/// Hypercomplex Right Division: `lhs / rhs = lhs * rhs.inv()`
template <nion_math::concepts::Scalar T, std::size_t N1, std::size_t N2, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr auto operator/(
    const nion<T, N1, Space>& lhs, const nion<T, N2, Space>& rhs)
{
    return lhs * rhs.inv();
}

/// scalar / nion: `s * q.inv()`
template <nion_math::concepts::Scalar T, std::size_t N, typename Space, nion_math::ConvertibleScalar<T> S>
NION_HOST_DEVICE constexpr nion<T, N, Space> operator/(
    const S& s, const nion<T, N, Space>& q)
{
    nion<T, N, Space> res = q.inv();
    res *= static_cast<T>(s);
    return res;
}

namespace nion_math {

// ============================================================================
// 3b. Free Function Equivalents (Matching Original Nion API)
// ============================================================================
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr T real(const ::nion<T, N, Space>& q) noexcept {
    return q.real();
}

template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr ::nion<T, N, Space> imag(const ::nion<T, N, Space>& q) {
    return q.imag();
}

template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr ::nion<T, N, Space> conj(const ::nion<T, N, Space>& q) {
    return q.conj();
}

template <typename Derived, concepts::Scalar T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr T abs(const NionExpr<Derived, T>& q) noexcept {
    return q.derived().abs();
}

template <typename Derived, concepts::Scalar T>
NION_HOST_DEVICE NION_FORCE_INLINE T norm(const NionExpr<Derived, T>& q) noexcept {
    return q.derived().norm();
}

template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr ::nion<T, N, Space> inv(const ::nion<T, N, Space>& q) {
    return q.inv();
}

template <concepts::Scalar T, std::size_t N1, std::size_t N2, typename Space>
NION_HOST_DEVICE constexpr T dot(
    const ::nion<T, N1, Space>& lhs, const ::nion<T, N2, Space>& rhs) noexcept
{
    const std::size_t min_sz = (lhs.size() < rhs.size()) ? lhs.size() : rhs.size();
    T acc{};
    for (std::size_t i = 0; i < min_sz; ++i) {
        acc = acc + lhs[i] * rhs[i];
    }
    return acc;
}

// ============================================================================
// 4. Convenient Fixed-Size Hypercomplex Type Aliases
// ============================================================================
template <concepts::Scalar T, typename Space = HostSpace> using real_nion       = ::nion<T, 1, Space>;
template <concepts::Scalar T, typename Space = HostSpace> using complex_t       = ::nion<T, 2, Space>;
template <concepts::Scalar T, typename Space = HostSpace> using quaternion      = ::nion<T, 4, Space>;
template <concepts::Scalar T, typename Space = HostSpace> using octonion        = ::nion<T, 8, Space>;
template <concepts::Scalar T, typename Space = HostSpace> using sedenion        = ::nion<T, 16, Space>;
template <concepts::Scalar T, typename Space = HostSpace> using trigintaduonion = ::nion<T, 32, Space>;

// ============================================================================
// 5. Non-Owning `nion_view<T, N>`: Zero-Copy GPU Kernel & Buffer View
// ============================================================================
/// Allows raw GPU device pointers (`T*`) or slices of tensors (`nionix`) to be
/// treated as first-class `nion` objects inside CUDA/HIP kernels with zero copies!
template <concepts::Scalar T, std::size_t N>
    requires (N > 0 && concepts::ValidDimension<N>)
class nion_view {
private:
    T* m_ptr;

public:
    NION_HOST_DEVICE constexpr explicit nion_view(T* ptr) noexcept : m_ptr(ptr) {}

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T* data() noexcept { return m_ptr; }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr const T* data() const noexcept { return m_ptr; }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr std::size_t size() const noexcept { return N; }

    NION_HOST_DEVICE NION_FORCE_INLINE constexpr T& operator[](std::size_t i) noexcept { return m_ptr[i]; }
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr const T& operator[](std::size_t i) const noexcept { return m_ptr[i]; }

    /// Load into a register-backed fixed `nion<T, N>`
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr ::nion<std::remove_cv_t<T>, N> load() const noexcept {
        return ::nion<std::remove_cv_t<T>, N>(m_ptr, N);
    }

    /// Store a register-backed fixed `nion<T, N>` back into the underlying buffer
    NION_HOST_DEVICE NION_FORCE_INLINE constexpr void store(const ::nion<std::remove_cv_t<T>, N>& val) noexcept
        requires (!std::is_const_v<T>)
    {
        detail::copy_n(val.data(), N, m_ptr);
    }
};

// ============================================================================
// 6. Structure-of-Arrays (`nion_batch_soa<T, N, Space>`) for Peak SIMD & GPU Speed
// ============================================================================
/// Stores a batch of `M` hypercomplex numbers of dimension `N` in Structure-of-Arrays
/// (SoA) layout:
///   Component e_0: [q_{0,0}, q_{1,0}, ..., q_{M-1,0}]
///   Component e_1: [q_{0,1}, q_{1,1}, ..., q_{M-1,1}]
///   ...
/// Benefits:
///   1. On CPU (AVX2 / AVX-512 / NEON): Every basis operation uses 100% full-width
///      vector FMA instructions across 8 floats / 4 doubles simultaneously with ZERO
///      cross-lane permutation overhead!
///   2. On GPU (CUDA / HIP): Warp threads `tid = 0..31` access contiguous addresses
///      `component_ptr(k)[tid]`, achieving 100% coalesced global memory bandwidth!
template <concepts::Scalar T, std::size_t N, typename Space = HostSpace>
    requires (N > 0 && concepts::ValidDimension<N>)
class nion_batch_soa {
private:
    std::size_t m_count  = 0; // Logical number of Nions in batch
    std::size_t m_stride = 0; // Padded stride (multiple of 16 elements for 64B alignment)
    T* m_buffer          = nullptr;

    static constexpr std::size_t pad_stride(std::size_t count) noexcept {
        constexpr std::size_t lane_align = 64 / sizeof(T) > 0 ? (64 / sizeof(T)) : 1;
        return ((count + lane_align - 1) / lane_align) * lane_align;
    }

public:
    using value_type                         = T;
    using space_type                         = Space;
    static constexpr std::size_t dimension   = N;

    constexpr nion_batch_soa() noexcept = default;

    explicit nion_batch_soa(std::size_t batch_count)
        : m_count(batch_count),
          m_stride(pad_stride(batch_count))
    {
        if (m_stride > 0) {
            m_buffer = Space::template allocate<T>(N * m_stride);
            if constexpr (Space::host_accessible) {
                detail::fill_n(m_buffer, N * m_stride, T{});
            }
        }
    }

    ~nion_batch_soa() noexcept {
        Space::deallocate(m_buffer, N * m_stride);
    }

    nion_batch_soa(const nion_batch_soa& other)
        : m_count(other.m_count), m_stride(other.m_stride)
    {
        if (m_stride > 0) {
            m_buffer = Space::template allocate<T>(N * m_stride);
            Space::copy(m_buffer, other.m_buffer, N * m_stride);
        }
    }

    nion_batch_soa(nion_batch_soa&& other) noexcept
        : m_count(std::exchange(other.m_count, 0)),
          m_stride(std::exchange(other.m_stride, 0)),
          m_buffer(std::exchange(other.m_buffer, nullptr))
    {}

    NION_HOST_DEVICE NION_FORCE_INLINE std::size_t count() const noexcept { return m_count; }
    NION_HOST_DEVICE NION_FORCE_INLINE std::size_t stride() const noexcept { return m_stride; }

    /// Pointer to the contiguous array of the k-th hypercomplex basis component (e_k)
    NION_HOST_DEVICE NION_FORCE_INLINE T* component(std::size_t k) noexcept {
        return m_buffer + k * m_stride;
    }
    NION_HOST_DEVICE NION_FORCE_INLINE const T* component(std::size_t k) const noexcept {
        return m_buffer + k * m_stride;
    }

    /// Set the idx-th Nion in the batch from an AoS `nion<T, N>`
    NION_HOST_DEVICE NION_FORCE_INLINE void set(std::size_t idx, const ::nion<T, N>& q) noexcept {
        for (std::size_t k = 0; k < N; ++k) {
            component(k)[idx] = q[k];
        }
    }

    /// Gather the idx-th Nion in the batch into an AoS `nion<T, N>`
    NION_HOST_DEVICE NION_FORCE_INLINE ::nion<T, N> get(std::size_t idx) const noexcept {
        ::nion<T, N> q;
        for (std::size_t k = 0; k < N; ++k) {
            q[k] = component(k)[idx];
        }
        return q;
    }

    /// High-Throughput SIMD Batch Multiplication (`out = a * b`) across all `M` elements!
    /// Zero cross-lane shuffles: compilers vectorize the inner `i` loop with pure FMA instructions.
    static void multiply_batch(
        nion_batch_soa& out, const nion_batch_soa& a, const nion_batch_soa& b) noexcept
    {
        assert(out.count() == a.count() && a.count() == b.count());
        const std::size_t M = out.count();

        for (std::size_t k = 0; k < N; ++k) {
            T* NION_RESTRICT out_k = out.component(k);
            detail::fill_n(out_k, M, T{});

            for (std::size_t p = 0; p < N; ++p) {
                const std::size_t q = p ^ k;
                const int8_t sign = simd::cd_table<N>.signs[k][p];
                const T* NION_RESTRICT a_p = a.component(p);
                const T* NION_RESTRICT b_q = b.component(q);

                if (sign > 0) {
                    for (std::size_t i = 0; i < M; ++i) {
                        out_k[i] = out_k[i] + a_p[i] * b_q[i];
                    }
                } else {
                    for (std::size_t i = 0; i < M; ++i) {
                        out_k[i] = out_k[i] - a_p[i] * b_q[i];
                    }
                }
            }
        }
    }
};

} // namespace nion_math

// Expose convenient fixed-dimension aliases, free functions, and batch/view types at global scope
using nion_math::real;
using nion_math::imag;
using nion_math::conj;
using nion_math::abs;
using nion_math::norm;
using nion_math::inv;
using nion_math::dot;
using nion_math::real_nion;
using nion_math::complex_t;
using nion_math::quaternion;
using nion_math::octonion;
using nion_math::sedenion;
using nion_math::trigintaduonion;
using nion_math::nion_view;
using nion_math::nion_batch_soa;

#include "nion_trig.hpp"

#endif // NION_HPP
