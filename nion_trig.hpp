// Copyright 2022-2026 Marcus Dante Liebenthal & Contributors
// Licensed under the Apache License, Version 2.0 (the "License");
//
// nion_trig.hpp - Hypercomplex Polar Decomposition, Exponential, Power, and Trigonometry
// Header-only, Host & GPU Device (`__host__ __device__`) compatible.
// Every hypercomplex number q = a_0 + v (where v = imag(q)) spans a 2D complex-like
// subplane {1, u} with unit pure-imaginary versor u = v / ||v|| and u^2 = -1.
// Analytic functions f(z) extend naturally to any 2^k-nion via:
//   f(a_0 + ||v|| u) = Re(f(a_0 + i ||v||)) + Im(f(a_0 + i ||v||)) * u

#ifndef NION_TRIG_HPP
#define NION_TRIG_HPP

#include "nion.hpp"
#include <cmath>
#include <limits>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

#ifndef M_PI_2
#define M_PI_2 1.570796326794896619231321691639751442L
#endif

namespace nion_math {

/// Holds the polar / complex-plane projection of a hypercomplex number:
///   q = r * (cos(theta) + u * sin(theta)) = a_0 + v_norm * u
template <concepts::Scalar T, std::size_t N, typename Space>
struct PolarDecomposition {
    T real_part;                 // a_0 = Re(q)
    T imag_norm;                 // ||v|| = sqrt(sum_{i=1}^{n-1} a_i^2)
    T magnitude;                 // |q| = hypot(a_0, ||v||)
    T angle;                     // theta = atan2(||v||, a_0) in [0, pi]
    nion<T, N, Space> versor;    // Unit pure-imaginary axis u = v / ||v|| (or e_1 if ||v|| == 0)
};

namespace detail {

template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> from_complex_projection(
    const nion<T, N, Space>& q, T re_out, T im_over_vnorm)
{
    const std::size_t len = q.size();
    nion<T, N, Space> res(uninit_tag{}, len);
    res[0] = re_out;
    for (std::size_t i = 1; i < len; ++i) {
        res[i] = q[i] * im_over_vnorm;
    }
    return res;
}

template <concepts::Scalar T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr T imag_norm_sq(const T* data, std::size_t len) noexcept {
    if (len <= 1) return T{};
    if (len == 2) return data[1] * data[1];
    return simd::vec_norm_sq(data + 1, len - 1);
}

template <concepts::Scalar T>
NION_HOST_DEVICE NION_FORCE_INLINE T imag_norm(const T* data, std::size_t len) noexcept {
    using std::abs;
    using std::sqrt;
    if (len <= 1) return T{};
    if constexpr (std::is_arithmetic_v<T>) {
        if (len == 2) return static_cast<T>(abs(data[1]));
    }
    return static_cast<T>(sqrt(simd::vec_norm_sq(data + 1, len - 1)));
}

template <concepts::Scalar T>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr T safe_denorm_min() noexcept {
    if constexpr (std::is_arithmetic_v<T>) {
        return std::numeric_limits<T>::denorm_min();
    } else {
        return T{};
    }
}

} // namespace detail

/// Computes the general angle and unit imaginary axis in the hypercomplex plane (`atrig`)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE PolarDecomposition<T, N, Space> atrig(const nion<T, N, Space>& q) {
    using std::atan2;
    using std::sqrt;

    PolarDecomposition<T, N, Space> pd{};
    pd.real_part = q.real();

    const std::size_t len = q.size();
    const T v_norm_sq = detail::imag_norm_sq(q.data(), len);
    pd.imag_norm = (len == 2 && std::is_arithmetic_v<T>)
        ? detail::imag_norm(q.data(), len)
        : static_cast<T>(sqrt(v_norm_sq));
    pd.magnitude = static_cast<T>(sqrt(pd.real_part * pd.real_part + v_norm_sq));
    pd.angle     = static_cast<T>(atan2(pd.imag_norm, pd.real_part));

    pd.versor = nion<T, N, Space>(len, T{});
    const T eps = detail::safe_denorm_min<T>();
    if (pd.imag_norm > eps) {
        const T inv_v = T{1} / pd.imag_norm;
        for (std::size_t i = 1; i < len; ++i) {
            pd.versor[i] = q[i] * inv_v;
        }
    } else if (len > 1) {
        pd.versor[1] = T{1};
    }
    return pd;
}

// ============================================================================
// 1. Exponential, Logarithm, Power, and Roots (exp, log, pow, sqr, sqrt, cbrt)
// ============================================================================

/// Hypercomplex Exponential: exp(a_0 + v) = exp(a_0) * (cos(||v||) + v * (sin(||v||) / ||v||))
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> exp(const nion<T, N, Space>& q) {
    using std::cos;
    using std::exp;
    using std::sin;

    const T a0 = q.real();
    const T exp_a0 = static_cast<T>(exp(a0));

    if constexpr (N == 2 && std::is_arithmetic_v<T>) {
        const T b = q[1];
        const T sb = static_cast<T>(sin(b));
        const T cb = static_cast<T>(cos(b));
        nion<T, N, Space> res(detail::uninit_tag{}, 2);
        res[0] = exp_a0 * cb;
        res[1] = exp_a0 * sb;
        return res;
    } else {
        const T v_norm = detail::imag_norm(q.data(), q.size());
        const T sv = static_cast<T>(sin(v_norm));
        const T cv = static_cast<T>(cos(v_norm));
        const T eps = detail::safe_denorm_min<T>();
        const T scale = (v_norm > eps) ? (exp_a0 * (sv / v_norm)) : exp_a0;
        return detail::from_complex_projection(q, exp_a0 * cv, scale);
    }
}

/// Hypercomplex Natural Logarithm: log(q) = 0.5*log(|q|^2) + v * (atan2(||v||, a_0) / ||v||)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> log(const nion<T, N, Space>& q) {
    using std::atan2;
    using std::log;
    using std::sqrt;

    const T a0 = q.real();
    if constexpr (N == 2 && std::is_arithmetic_v<T>) {
        const T b = q[1];
        nion<T, N, Space> res(detail::uninit_tag{}, 2);
        res[0] = static_cast<T>(T{0.5} * log(a0 * a0 + b * b));
        res[1] = static_cast<T>(atan2(b, a0));
        return res;
    } else {
        const std::size_t len = q.size();
        const T v_norm_sq = detail::imag_norm_sq(q.data(), len);
        const T v_norm = static_cast<T>(sqrt(v_norm_sq));
        const T z_abs = a0 * a0 + v_norm_sq;
        const T re_out = static_cast<T>(T{0.5} * log(z_abs));
        const T theta = static_cast<T>(atan2(v_norm, a0));

        const T eps = detail::safe_denorm_min<T>();
        const T scale = (v_norm > eps) ? (theta / v_norm) : theta;
        return detail::from_complex_projection(q, re_out, scale);
    }
}

/// Hypercomplex Square: sqr(a_0 + v) = (a_0^2 - ||v||^2) + (2 * a_0) * v
/// Evaluates in O(N) instead of O(N^2)!
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE constexpr nion<T, N, Space> sqr(const nion<T, N, Space>& q) {
    if constexpr (std::is_arithmetic_v<T>) {
        const T a0 = q.real();
        const T v_sq = detail::imag_norm_sq(q.data(), q.size());
        return detail::from_complex_projection(q, a0 * a0 - v_sq, a0 + a0);
    } else {
        return q * q;
    }
}

/// Hypercomplex Principal Square Root: sqrt(q)
/// Uses algebraic 2D complex square root requiring only 2 scalar sqrts (no atan2/sin/cos)!
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> sqrt(const nion<T, N, Space>& q) {
    using std::sqrt;

    const T x = q.real();
    const std::size_t len = q.size();
    const T y2 = detail::imag_norm_sq(q.data(), len);
    const T eps = detail::safe_denorm_min<T>();

    if (y2 <= eps) {
        if (x >= T{0}) {
            return detail::from_complex_projection(q, static_cast<T>(sqrt(x)), T{0});
        } else {
            nion<T, N, Space> res(len, T{});
            if (len > 1) res[1] = static_cast<T>(sqrt(-x));
            return res;
        }
    }

    const T r = static_cast<T>(sqrt(x * x + y2));
    if (x >= T{0}) {
        const T w = static_cast<T>(sqrt(T{0.5} * (r + x)));
        return detail::from_complex_projection(q, w, T{0.5} / w);
    } else {
        const T u = static_cast<T>(sqrt(T{0.5} * (r - x)));
        const T y = (len == 2 && std::is_arithmetic_v<T>)
            ? detail::imag_norm(q.data(), len)
            : static_cast<T>(sqrt(y2));
        return detail::from_complex_projection(q, T{0.5} * y / u, u / y);
    }
}

/// Hypercomplex Principal Cube Root: cbrt(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> cbrt(const nion<T, N, Space>& q) {
    using std::atan2;
    using std::cbrt;
    using std::cos;
    using std::sin;
    using std::sqrt;

    const T x = q.real();
    const std::size_t len = q.size();
    const T y2 = detail::imag_norm_sq(q.data(), len);
    const T y = (len == 2 && std::is_arithmetic_v<T>)
        ? detail::imag_norm(q.data(), len)
        : static_cast<T>(sqrt(y2));
    const T r = static_cast<T>(sqrt(x * x + y2));
    const T root_r = static_cast<T>(cbrt(r));
    const T third_theta = static_cast<T>(atan2(y, x)) / T{3};

    const T eps = detail::safe_denorm_min<T>();
    const T scale = (y > eps)
        ? (root_r * static_cast<T>(sin(third_theta)) / y)
        : T{0};
    return detail::from_complex_projection(
        q, root_r * static_cast<T>(cos(third_theta)), scale);
}

/// Hypercomplex Power (nion^scalar): pow(q, p)
/// Includes exact O(N) polynomial fast paths for small integer powers p in {-2, -1, 0, 1, 2, 3, 4}
/// and algebraic fast path for p = 0.5.
template <concepts::Scalar T, std::size_t N, typename Space, ConvertibleScalar<T> S>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> pow(const nion<T, N, Space>& q, const S& exponent) {
    using std::atan2;
    using std::cos;
    using std::exp;
    using std::log;
    using std::sin;
    using std::sqrt;

    const T p = static_cast<T>(exponent);

    if (p == T{0}) {
        nion<T, N, Space> one(q.size(), T{});
        one[0] = T{1};
        return one;
    }
    if (p == T{1}) return q;
    if (p == T{2}) return sqr(q);
    if (p == T{-1}) return q.inv();
    if (p == T{-2}) return sqr(q.inv());

    if constexpr (std::is_arithmetic_v<T>) {
        if (p == T{3}) {
            const T x = q.real();
            const T x2 = x * x;
            const T y2 = detail::imag_norm_sq(q.data(), q.size());
            return detail::from_complex_projection(
                q, x * (x2 - T{3} * y2), T{3} * x2 - y2);
        }
        if (p == T{4}) {
            const T x = q.real();
            const T x2 = x * x;
            const T y2 = detail::imag_norm_sq(q.data(), q.size());
            const T d = x2 - y2;
            return detail::from_complex_projection(
                q, d * d - T{4} * x2 * y2, T{4} * x * d);
        }
        if (p == T{0.5}) {
            return sqrt(q);
        }
        if (p == T{1.5}) {
            const nion<T, N, Space> s = sqrt(q);
            const T x = q.real();
            const T y2 = detail::imag_norm_sq(q.data(), q.size());
            const T w = s.real();
            // Since s = w + scale_s * v where scale_s = 0.5 / w (for x >= 0),
            // q * s = (x + v)(w + scale_s * v) = (x*w - scale_s*y2) + (w + x*scale_s)*v
            const T eps = detail::safe_denorm_min<T>();
            if (y2 <= eps) {
                return q * s;
            }
            const T scale_s = (x >= T{0}) ? (T{0.5} / w) : (w / (T{2} * (q.size() > 1 ? T{1} : T{1})));
            (void)scale_s;
            // Direct 2D complex multiplication of (x, y) and (s.real(), ||imag(s)||)
            const T y = detail::imag_norm(q.data(), q.size());
            const T su = (x >= T{0}) ? (T{0.5} * y / w) : static_cast<T>(sqrt(T{0.5} * (static_cast<T>(sqrt(x * x + y2)) - x)));
            return detail::from_complex_projection(
                q, x * w - y * su, w + x * (su / y));
        }
    }

    const T x = q.real();
    const std::size_t len = q.size();
    const T y2 = detail::imag_norm_sq(q.data(), len);
    const T r2 = x * x + y2;
    const T eps = detail::safe_denorm_min<T>();
    if (r2 <= eps) {
        return nion<T, N, Space>(len, T{});
    }

    const T y = (len == 2 && std::is_arithmetic_v<T>)
        ? detail::imag_norm(q.data(), len)
        : static_cast<T>(sqrt(y2));
    const T mag_p = static_cast<T>(exp(T{0.5} * p * static_cast<T>(log(r2))));
    const T theta_p = p * static_cast<T>(atan2(y, x));

    if (y <= eps) {
        nion<T, N, Space> res(len, T{});
        res[0] = mag_p * static_cast<T>(cos(theta_p));
        if (len > 1) res[1] = mag_p * static_cast<T>(sin(theta_p));
        return res;
    }

    return detail::from_complex_projection(
        q,
        mag_p * static_cast<T>(cos(theta_p)),
        mag_p * static_cast<T>(sin(theta_p)) / y);
}

/// Hypercomplex Power (nion^nion): pow(base, exp_nion) = exp(log(base) * exp_nion)
template <concepts::Scalar T, std::size_t N1, std::size_t N2, typename Space>
NION_HOST_DEVICE auto pow(const nion<T, N1, Space>& base, const nion<T, N2, Space>& exp_nion) {
    return exp(log(base) * exp_nion);
}

/// Hypercomplex Power (scalar^nion): pow(base, exp_nion) = exp(log(base) * exp_nion)
template <concepts::Scalar T, std::size_t N, typename Space, ConvertibleScalar<T> S>
NION_HOST_DEVICE nion<T, N, Space> pow(const S& base, const nion<T, N, Space>& exp_nion) {
    using std::log;
    const T b = static_cast<T>(base);
    if constexpr (std::is_arithmetic_v<T>) {
        if (b > T{0}) {
            return exp(static_cast<T>(log(b)) * exp_nion);
        }
    }
    nion<T, N, Space> b_nion(exp_nion.size(), T{});
    b_nion[0] = b;
    return exp(log(b_nion) * exp_nion);
}

// ============================================================================
// 2. Hyperbolic Functions (sinh, cosh, tanh, coth, sech, csch & inverses)
// ============================================================================

/// Hypercomplex Hyperbolic Sine:
///   sinh(a_0 + v) = sinh(a_0) cos(||v||) + cosh(a_0) (sin(||v||) / ||v||) v
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> sinh(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T sh = static_cast<T>(sinh(a0));
    const T ch = static_cast<T>(cosh(a0));

    if constexpr (N == 2 && std::is_arithmetic_v<T>) {
        const T b = q[1];
        const T sb = static_cast<T>(sin(b));
        const T cb = static_cast<T>(cos(b));
        nion<T, N, Space> res(detail::uninit_tag{}, 2);
        res[0] = sh * cb;
        res[1] = ch * sb;
        return res;
    } else {
        const T v_norm = detail::imag_norm(q.data(), q.size());
        const T sv = static_cast<T>(sin(v_norm));
        const T cv = static_cast<T>(cos(v_norm));
        const T eps = detail::safe_denorm_min<T>();
        const T sinc = (v_norm > eps) ? (sv / v_norm) : T{1};
        return detail::from_complex_projection(q, sh * cv, ch * sinc);
    }
}

/// Hypercomplex Hyperbolic Cosine:
///   cosh(a_0 + v) = cosh(a_0) cos(||v||) + sinh(a_0) (sin(||v||) / ||v||) v
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> cosh(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T sh = static_cast<T>(sinh(a0));
    const T ch = static_cast<T>(cosh(a0));

    if constexpr (N == 2 && std::is_arithmetic_v<T>) {
        const T b = q[1];
        const T sb = static_cast<T>(sin(b));
        const T cb = static_cast<T>(cos(b));
        nion<T, N, Space> res(detail::uninit_tag{}, 2);
        res[0] = ch * cb;
        res[1] = sh * sb;
        return res;
    } else {
        const T v_norm = detail::imag_norm(q.data(), q.size());
        const T sv = static_cast<T>(sin(v_norm));
        const T cv = static_cast<T>(cos(v_norm));
        const T eps = detail::safe_denorm_min<T>();
        const T sinc = (v_norm > eps) ? (sv / v_norm) : T{1};
        return detail::from_complex_projection(q, ch * cv, sh * sinc);
    }
}

/// Hypercomplex Hyperbolic Tangent:
/// Direct 2D complex formula (zero hypercomplex division!):
///   denom = sinh^2(a_0) + cos^2(||v||)
///   tanh(a_0 + v) = (sinh(a_0) cosh(a_0) + (sin(||v||)/||v||) cos(||v||) v) / denom
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> tanh(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T sinhr = static_cast<T>(sinh(a0));
    const T coshr = static_cast<T>(cosh(a0));

    if constexpr (N == 2 && std::is_arithmetic_v<T>) {
        const T b = q[1];
        const T sinb = static_cast<T>(sin(b));
        const T cosb = static_cast<T>(cos(b));
        const T inv_denom = T{1} / (sinhr * sinhr + cosb * cosb);
        nion<T, N, Space> res(detail::uninit_tag{}, 2);
        res[0] = (sinhr * coshr) * inv_denom;
        res[1] = (sinb * cosb) * inv_denom;
        return res;
    } else {
        const T y = detail::imag_norm(q.data(), q.size());
        const T siny = static_cast<T>(sin(y));
        const T cosy = static_cast<T>(cos(y));
        const T inv_denom = T{1} / (sinhr * sinhr + cosy * cosy);
        const T eps = detail::safe_denorm_min<T>();
        const T sinc = (y > eps) ? (siny / y) : T{1};
        return detail::from_complex_projection(
            q, (sinhr * coshr) * inv_denom, (sinc * cosy) * inv_denom);
    }
}

/// Hypercomplex Hyperbolic Cotangent: coth(q) = cosh(q) / sinh(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> coth(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T y = detail::imag_norm(q.data(), q.size());
    const T sinhr = static_cast<T>(sinh(a0));
    const T coshr = static_cast<T>(cosh(a0));
    const T siny  = static_cast<T>(sin(y));
    const T cosy  = static_cast<T>(cos(y));
    const T denom = sinhr * sinhr + siny * siny;

    const T eps = detail::safe_denorm_min<T>();
    if (y <= eps) {
        return detail::from_complex_projection(q, coshr / sinhr, T{0});
    }
    const T inv_denom = T{1} / denom;
    return detail::from_complex_projection(
        q,
        (sinhr * coshr) * inv_denom,
        -(cosy * (siny / y)) * inv_denom);
}

/// Hypercomplex Hyperbolic Secant: sech(q) = 1 / cosh(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> sech(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T y = detail::imag_norm(q.data(), q.size());
    const T sinhr = static_cast<T>(sinh(a0));
    const T coshr = static_cast<T>(cosh(a0));
    const T cosy  = static_cast<T>(cos(y));
    const T denom = sinhr * sinhr + cosy * cosy;

    const T eps = detail::safe_denorm_min<T>();
    const T sinc = (y > eps) ? (static_cast<T>(sin(y)) / y) : T{1};
    const T inv_denom = T{1} / denom;

    return detail::from_complex_projection(
        q,
        (coshr * cosy) * inv_denom,
        -(sinhr * sinc) * inv_denom);
}

/// Hypercomplex Hyperbolic Cosecant: csch(q) = 1 / sinh(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> csch(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T y = detail::imag_norm(q.data(), q.size());
    const T sinhr = static_cast<T>(sinh(a0));
    const T coshr = static_cast<T>(cosh(a0));
    const T siny  = static_cast<T>(sin(y));
    const T cosy  = static_cast<T>(cos(y));
    const T denom = sinhr * sinhr + siny * siny;

    const T eps = detail::safe_denorm_min<T>();
    if (y <= eps) {
        return detail::from_complex_projection(q, T{1} / sinhr, T{0});
    }
    const T inv_denom = T{1} / denom;
    return detail::from_complex_projection(
        q,
        (sinhr * cosy) * inv_denom,
        -(coshr * (siny / y)) * inv_denom);
}

/// Inverse Hyperbolic Sine: asinh(q)
/// Direct 2D Hull-Fairgrieve-Tang / Harder complex projection (zero O(N^2) multiplications!)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> asinh(const nion<T, N, Space>& q) {
    using std::asin;
    using std::log;
    using std::sqrt;

    const T r = q.real();
    const std::size_t len = q.size();
    const T y2 = detail::imag_norm_sq(q.data(), len);
    const T y = (len == 2 && std::is_arithmetic_v<T>)
        ? detail::imag_norm(q.data(), len)
        : static_cast<T>(sqrt(y2));
    const T eps = detail::safe_denorm_min<T>();

    if (y <= eps) {
        return detail::from_complex_projection(
            q, static_cast<T>(log(r + static_cast<T>(sqrt(r * r + T{1})))), T{1});
    }

    const T ss = r * r + y2 + T{1};
    const T ssp2y = static_cast<T>(sqrt(ss + T{2} * y));
    const T ssm2y = static_cast<T>(sqrt(ss - T{2} * y));
    const T sum = T{0.5} * (ssp2y + ssm2y);
    const T sign_r = (r < T{0}) ? T{-1} : T{1};

    const T re_out = sign_r * static_cast<T>(log(sum + static_cast<T>(sqrt(sum * sum - T{1}))));
    const T im_out = static_cast<T>(asin(T{0.5} * (ssp2y - ssm2y)));

    return detail::from_complex_projection(q, re_out, im_out / y);
}

/// Inverse Hyperbolic Cosine: acosh(q)
/// Direct 2D Hull-Fairgrieve-Tang / Harder complex projection with principal branch cut (-oo, 1)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> acosh(const nion<T, N, Space>& q) {
    using std::acos;
    using std::log;
    using std::sqrt;

    const T r = q.real();
    const std::size_t len = q.size();
    const T y2 = detail::imag_norm_sq(q.data(), len);
    const T y = (len == 2 && std::is_arithmetic_v<T>)
        ? detail::imag_norm(q.data(), len)
        : static_cast<T>(sqrt(y2));
    const T eps = detail::safe_denorm_min<T>();

    if (y <= eps) {
        if (r >= T{1}) {
            return detail::from_complex_projection(
                q, static_cast<T>(log(r + static_cast<T>(sqrt(r * r - T{1})))), T{0});
        } else if (r >= T{-1}) {
            nion<T, N, Space> res(len, T{});
            if (len > 1) res[1] = static_cast<T>(acos(r));
            return res;
        } else {
            nion<T, N, Space> res(len, T{});
            res[0] = static_cast<T>(log(-r + static_cast<T>(sqrt(r * r - T{1}))));
            if (len > 1) res[1] = static_cast<T>(M_PI);
            return res;
        }
    }

    const T ss = r * r + y2 + T{1};
    const T ssp2r = static_cast<T>(sqrt(ss + T{2} * r));
    const T ssm2r = static_cast<T>(sqrt(ss - T{2} * r));
    const T sum = T{0.5} * (ssp2r + ssm2r);

    const T re_out = static_cast<T>(log(sum + static_cast<T>(sqrt(sum * sum - T{1}))));
    const T im_out = static_cast<T>(acos(T{0.5} * (ssp2r - ssm2r)));

    return detail::from_complex_projection(q, re_out, im_out / y);
}

/// Inverse Hyperbolic Tangent: atanh(q) = 0.5 * (log(1 + q) - log(1 - q))
/// Direct 2D complex projection using a single atan2 call (zero O(N) temporaries!)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> atanh(const nion<T, N, Space>& q) {
    using std::atan2;
    using std::log;
    using std::sqrt;

    const T r = q.real();
    if constexpr (N == 2 && std::is_arithmetic_v<T>) {
        const T b = q[1];
        const T b2 = b * b;
        const T opr = T{1} + r;
        const T omr = T{1} - r;
        nion<T, N, Space> res(detail::uninit_tag{}, 2);
        res[0] = T{0.25} * static_cast<T>(log((b2 + opr * opr) / (b2 + omr * omr)));
        res[1] = T{0.5} * static_cast<T>(atan2(b + b, T{1} - r * r - b2));
        return res;
    } else {
        const std::size_t len = q.size();
        const T y2 = detail::imag_norm_sq(q.data(), len);
        const T y = static_cast<T>(sqrt(y2));
        const T eps = detail::safe_denorm_min<T>();

        const T opr = T{1} + r;
        const T omr = T{1} - r;
        const T re_out = T{0.25} * static_cast<T>(log((y2 + opr * opr) / (y2 + omr * omr)));

        if (y <= eps) {
            if (r > T{-1} && r < T{1}) {
                return detail::from_complex_projection(q, re_out, T{1});
            }
            nion<T, N, Space> res(len, T{});
            res[0] = re_out;
            if (len > 1) res[1] = static_cast<T>(M_PI_2);
            return res;
        }

        const T im_out = T{0.5} * static_cast<T>(atan2(y + y, T{1} - r * r - y2));
        return detail::from_complex_projection(q, re_out, im_out / y);
    }
}

/// Inverse Hyperbolic Cotangent: acoth(q) = atanh(q.inv())
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> acoth(const nion<T, N, Space>& q) {
    return atanh(q.inv());
}

/// Inverse Hyperbolic Secant: asech(q) = acosh(q.inv())
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> asech(const nion<T, N, Space>& q) {
    return acosh(q.inv());
}

/// Inverse Hyperbolic Cosecant: acsch(q) = asinh(q.inv())
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> acsch(const nion<T, N, Space>& q) {
    return asinh(q.inv());
}

// ============================================================================
// 3. Circular Trigonometric Functions (sin, cos, tan, cot, sec, csc & inverses)
// ============================================================================

/// Hypercomplex Sine:
///   sin(a_0 + v) = sin(a_0) cosh(||v||) + cos(a_0) (sinh(||v||) / ||v||) v
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> sin(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T sa = static_cast<T>(sin(a0));
    const T ca = static_cast<T>(cos(a0));

    if constexpr (N == 2 && std::is_arithmetic_v<T>) {
        const T b = q[1];
        const T shb = static_cast<T>(sinh(b));
        const T chb = static_cast<T>(cosh(b));
        nion<T, N, Space> res(detail::uninit_tag{}, 2);
        res[0] = sa * chb;
        res[1] = ca * shb;
        return res;
    } else {
        const T v_norm = detail::imag_norm(q.data(), q.size());
        const T shv = static_cast<T>(sinh(v_norm));
        const T chv = static_cast<T>(cosh(v_norm));
        const T eps = detail::safe_denorm_min<T>();
        const T sinhc = (v_norm > eps) ? (shv / v_norm) : T{1};
        return detail::from_complex_projection(q, sa * chv, ca * sinhc);
    }
}

/// Hypercomplex Cosine:
///   cos(a_0 + v) = cos(a_0) cosh(||v||) - sin(a_0) (sinh(||v||) / ||v||) v
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> cos(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T sa = static_cast<T>(sin(a0));
    const T ca = static_cast<T>(cos(a0));

    if constexpr (N == 2 && std::is_arithmetic_v<T>) {
        const T b = q[1];
        const T shb = static_cast<T>(sinh(b));
        const T chb = static_cast<T>(cosh(b));
        nion<T, N, Space> res(detail::uninit_tag{}, 2);
        res[0] = ca * chb;
        res[1] = -sa * shb;
        return res;
    } else {
        const T v_norm = detail::imag_norm(q.data(), q.size());
        const T shv = static_cast<T>(sinh(v_norm));
        const T chv = static_cast<T>(cosh(v_norm));
        const T eps = detail::safe_denorm_min<T>();
        const T sinhc = (v_norm > eps) ? (shv / v_norm) : T{1};
        return detail::from_complex_projection(q, ca * chv, -sa * sinhc);
    }
}

/// Hypercomplex Tangent:
/// Direct 2D complex formula (zero hypercomplex division!):
///   denom = cos^2(a_0) + sinh^2(||v||)
///   tan(a_0 + v) = (sin(a_0) cos(a_0) + (sinh(||v||)/||v||) cosh(||v||) v) / denom
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> tan(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T sinr = static_cast<T>(sin(a0));
    const T cosr = static_cast<T>(cos(a0));

    if constexpr (N == 2 && std::is_arithmetic_v<T>) {
        const T b = q[1];
        const T sinhb = static_cast<T>(sinh(b));
        const T coshb = static_cast<T>(cosh(b));
        const T inv_denom = T{1} / (cosr * cosr + sinhb * sinhb);
        nion<T, N, Space> res(detail::uninit_tag{}, 2);
        res[0] = (sinr * cosr) * inv_denom;
        res[1] = (sinhb * coshb) * inv_denom;
        return res;
    } else {
        const T y = detail::imag_norm(q.data(), q.size());
        const T sinhy = static_cast<T>(sinh(y));
        const T coshy = static_cast<T>(cosh(y));
        const T inv_denom = T{1} / (cosr * cosr + sinhy * sinhy);
        const T eps = detail::safe_denorm_min<T>();
        const T sinhc = (y > eps) ? (sinhy / y) : T{1};
        return detail::from_complex_projection(
            q, (sinr * cosr) * inv_denom, (sinhc * coshy) * inv_denom);
    }
}

/// Hypercomplex Cotangent: cot(q) = cos(q) / sin(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> cot(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T y = detail::imag_norm(q.data(), q.size());
    const T sinr  = static_cast<T>(sin(a0));
    const T cosr  = static_cast<T>(cos(a0));
    const T sinhy = static_cast<T>(sinh(y));
    const T coshy = static_cast<T>(cosh(y));
    const T denom = sinr * sinr + sinhy * sinhy;

    const T eps = detail::safe_denorm_min<T>();
    const T sinhc = (y > eps) ? (sinhy / y) : T{1};
    const T inv_denom = T{1} / denom;

    return detail::from_complex_projection(
        q,
        (sinr * cosr) * inv_denom,
        -(sinhc * coshy) * inv_denom);
}

/// Hypercomplex Secant: sec(q) = 1 / cos(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> sec(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T y = detail::imag_norm(q.data(), q.size());
    const T sinr  = static_cast<T>(sin(a0));
    const T cosr  = static_cast<T>(cos(a0));
    const T sinhy = static_cast<T>(sinh(y));
    const T coshy = static_cast<T>(cosh(y));
    const T denom = cosr * cosr + sinhy * sinhy;

    const T eps = detail::safe_denorm_min<T>();
    const T sinhc = (y > eps) ? (sinhy / y) : T{1};
    const T inv_denom = T{1} / denom;

    return detail::from_complex_projection(
        q,
        (cosr * coshy) * inv_denom,
        (sinr * sinhc) * inv_denom);
}

/// Hypercomplex Cosecant: csc(q) = 1 / sin(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> csc(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T y = detail::imag_norm(q.data(), q.size());
    const T sinr  = static_cast<T>(sin(a0));
    const T cosr  = static_cast<T>(cos(a0));
    const T sinhy = static_cast<T>(sinh(y));
    const T coshy = static_cast<T>(cosh(y));
    const T denom = sinr * sinr + sinhy * sinhy;

    const T eps = detail::safe_denorm_min<T>();
    const T sinhc = (y > eps) ? (sinhy / y) : T{1};
    const T inv_denom = T{1} / denom;

    return detail::from_complex_projection(
        q,
        (sinr * coshy) * inv_denom,
        -(cosr * sinhc) * inv_denom);
}

/// Inverse Sine: asin(q)
/// Direct 2D Hull-Fairgrieve-Tang / Harder complex projection (zero O(N^2) multiplications!)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> asin(const nion<T, N, Space>& q) {
    using std::asin;
    using std::log;
    using std::sqrt;

    const T r = q.real();
    const std::size_t len = q.size();
    const T y2 = detail::imag_norm_sq(q.data(), len);
    const T y = (len == 2 && std::is_arithmetic_v<T>)
        ? detail::imag_norm(q.data(), len)
        : static_cast<T>(sqrt(y2));
    const T eps = detail::safe_denorm_min<T>();

    if (y <= eps) {
        if (r >= T{-1} && r <= T{1}) {
            return detail::from_complex_projection(q, static_cast<T>(asin(r)), T{1});
        }
        nion<T, N, Space> res(len, T{});
        res[0] = (r > T{1}) ? static_cast<T>(M_PI_2) : -static_cast<T>(M_PI_2);
        const T abs_r = (r < T{0}) ? -r : r;
        if (len > 1) res[1] = static_cast<T>(log(abs_r + static_cast<T>(sqrt(r * r - T{1}))));
        return res;
    }

    const T ss = r * r + y2 + T{1};
    const T ssp2r = static_cast<T>(sqrt(ss + T{2} * r));
    const T ssm2r = static_cast<T>(sqrt(ss - T{2} * r));
    const T sum = T{0.5} * (ssp2r + ssm2r);

    const T re_out = static_cast<T>(asin(T{0.5} * (ssp2r - ssm2r)));
    const T im_out = static_cast<T>(log(sum + static_cast<T>(sqrt(sum * sum - T{1}))));

    return detail::from_complex_projection(q, re_out, im_out / y);
}

/// Inverse Cosine: acos(q)
/// Direct 2D Hull-Fairgrieve-Tang / Harder complex projection (zero O(N^2) multiplications!)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> acos(const nion<T, N, Space>& q) {
    using std::acos;
    using std::log;
    using std::sqrt;

    const T r = q.real();
    const std::size_t len = q.size();
    const T y2 = detail::imag_norm_sq(q.data(), len);
    const T y = (len == 2 && std::is_arithmetic_v<T>)
        ? detail::imag_norm(q.data(), len)
        : static_cast<T>(sqrt(y2));
    const T eps = detail::safe_denorm_min<T>();

    if (y <= eps) {
        if (r >= T{-1} && r <= T{1}) {
            return detail::from_complex_projection(q, static_cast<T>(acos(r)), T{-1});
        }
        nion<T, N, Space> res(len, T{});
        res[0] = (r > T{1}) ? T{0} : static_cast<T>(M_PI);
        const T abs_r = (r < T{0}) ? -r : r;
        if (len > 1) res[1] = -static_cast<T>(log(abs_r + static_cast<T>(sqrt(r * r - T{1}))));
        return res;
    }

    const T ss = r * r + y2 + T{1};
    const T ssp2r = static_cast<T>(sqrt(ss + T{2} * r));
    const T ssm2r = static_cast<T>(sqrt(ss - T{2} * r));
    const T sum = T{0.5} * (ssp2r + ssm2r);

    const T re_out = static_cast<T>(acos(T{0.5} * (ssp2r - ssm2r)));
    const T im_out = -static_cast<T>(log(sum + static_cast<T>(sqrt(sum * sum - T{1}))));

    return detail::from_complex_projection(q, re_out, im_out / y);
}

/// Inverse Tangent: atan(q)
/// Direct 2D complex projection using a single atan2 call (zero O(N^2) multiplications!)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE NION_FORCE_INLINE nion<T, N, Space> atan(const nion<T, N, Space>& q) {
    using std::atan;
    using std::atan2;
    using std::log;
    using std::sqrt;

    const T r = q.real();
    const T rr = r * r;

    if constexpr (N == 2 && std::is_arithmetic_v<T>) {
        const T b = q[1];
        const T opb = T{1} + b;
        const T omb = T{1} - b;
        nion<T, N, Space> res(detail::uninit_tag{}, 2);
        res[0] = T{0.5} * static_cast<T>(atan2(r + r, T{1} - rr - b * b));
        res[1] = T{0.25} * static_cast<T>(log((rr + opb * opb) / (rr + omb * omb)));
        return res;
    } else {
        const std::size_t len = q.size();
        const T y2 = detail::imag_norm_sq(q.data(), len);
        const T y = static_cast<T>(sqrt(y2));
        const T eps = detail::safe_denorm_min<T>();

        if (y <= eps) {
            return detail::from_complex_projection(q, static_cast<T>(atan(r)), T{1});
        }

        const T opy = T{1} + y;
        const T omy = T{1} - y;

        const T re_out = T{0.5} * static_cast<T>(atan2(r + r, T{1} - rr - y2));
        const T im_out = T{0.25} * static_cast<T>(log((rr + opy * opy) / (rr + omy * omy)));

        return detail::from_complex_projection(q, re_out, im_out / y);
    }
}

/// Inverse Cotangent: acot(q) = pi/2 - atan(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> acot(const nion<T, N, Space>& q) {
    return static_cast<T>(M_PI_2) - atan(q);
}

/// Inverse Secant: asec(q) = acos(q.inv())
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> asec(const nion<T, N, Space>& q) {
    return acos(q.inv());
}

/// Inverse Cosecant: acsc(q) = asin(q.inv())
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> acsc(const nion<T, N, Space>& q) {
    return asin(q.inv());
}

/// Two-argument Inverse Tangent: atan2(y, x) = atan(y / x)
template <concepts::Scalar T, std::size_t N1, std::size_t N2, typename Space>
NION_HOST_DEVICE auto atan2(const nion<T, N1, Space>& y, const nion<T, N2, Space>& x) {
    return atan(y / x);
}

/// General Inverse Trigonometric Helper (`atri`, matching original Nion API)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> atri(
    const nion<T, N, Space>& z,
    const nion<T, N, Space>& a,
    const nion<T, N, Space>& b,
    const nion<T, N, Space>& c)
{
    const auto pd = atrig(z);
    return (-pd.versor) * log((a + b * pd.versor) / c);
}

/// Hypercomplex Gamma Function (Nemes / Stirling approximation matching original Nion API)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> gamma(const nion<T, N, Space>& z) {
    using std::sqrt;
    const T root_2pi = static_cast<T>(sqrt(T{2} * static_cast<T>(M_PI)));
    const nion<T, N, Space> z10 = T{10} * z;
    const nion<T, N, Space> inner_denom = T{12} * z - z10.inv();
    const nion<T, N, Space> inner = inner_denom.inv() + z;
    return exp(-z) * sqrt(z.inv()) * pow(inner, z) * root_2pi;
}

} // namespace nion_math

// Expose analytic & trigonometric functions at global scope so unqualified
// calls (`sin(q)`, `exp(q)`, `pow(q, 2)`, `gamma(q)`) work seamlessly out of the box.
using nion_math::PolarDecomposition;
using nion_math::atrig;
using nion_math::atri;
using nion_math::exp;
using nion_math::log;
using nion_math::pow;
using nion_math::sqr;
using nion_math::sqrt;
using nion_math::cbrt;
using nion_math::sinh;
using nion_math::cosh;
using nion_math::tanh;
using nion_math::coth;
using nion_math::sech;
using nion_math::csch;
using nion_math::asinh;
using nion_math::acosh;
using nion_math::atanh;
using nion_math::acoth;
using nion_math::asech;
using nion_math::acsch;
using nion_math::sin;
using nion_math::cos;
using nion_math::tan;
using nion_math::cot;
using nion_math::sec;
using nion_math::csc;
using nion_math::asin;
using nion_math::acos;
using nion_math::atan;
using nion_math::atan2;
using nion_math::acot;
using nion_math::asec;
using nion_math::acsc;
using nion_math::gamma;

#endif // NION_TRIG_HPP
