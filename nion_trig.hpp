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

/// Computes the general angle and unit imaginary axis in the hypercomplex plane (`atrig`)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE PolarDecomposition<T, N, Space> atrig(const nion<T, N, Space>& q) {
    using std::atan2;
    using std::hypot;
    using std::sqrt;

    PolarDecomposition<T, N, Space> pd{};
    pd.real_part = q.real();

    const std::size_t len = q.size();
    T v_norm_sq{};
    if (len > 1) {
        v_norm_sq = simd::vec_norm_sq(q.data() + 1, len - 1);
    }
    pd.imag_norm = static_cast<T>(sqrt(v_norm_sq));
    pd.magnitude = static_cast<T>(hypot(pd.real_part, pd.imag_norm));
    pd.angle     = static_cast<T>(atan2(pd.imag_norm, pd.real_part));

    pd.versor = nion<T, N, Space>(len, T{});
    const T eps = std::numeric_limits<T>::epsilon() * T{16};
    if (pd.imag_norm > eps) {
        const T inv_v = T{1} / pd.imag_norm;
        for (std::size_t i = 1; i < len; ++i) {
            pd.versor[i] = q[i] * inv_v;
        }
    } else if (len > 1) {
        // Default canonical imaginary axis e_1 when q is purely real
        pd.versor[1] = T{1};
    }
    return pd;
}

/// Helper: constructs `re + im * u` where `v = q.imag()` and `u = v / ||v||`.
/// When `||v|| -> 0`, uses `sinc_scale = im / ||v||` directly to avoid 0/0 NaN singularities!
namespace detail {

template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> from_complex_projection(
    const nion<T, N, Space>& q, T re_out, T im_over_vnorm)
{
    const std::size_t len = q.size();
    nion<T, N, Space> res(len, T{});
    res[0] = re_out;
    for (std::size_t i = 1; i < len; ++i) {
        res[i] = q[i] * im_over_vnorm;
    }
    return res;
}

template <concepts::Scalar T>
NION_HOST_DEVICE T imag_norm(const T* data, std::size_t len) noexcept {
    using std::sqrt;
    if (len <= 1) return T{};
    return static_cast<T>(sqrt(simd::vec_norm_sq(data + 1, len - 1)));
}

} // namespace detail

// ============================================================================
// 1. Exponential, Logarithm, Power, and Roots (exp, log, pow, sqr, sqrt, cbrt)
// ============================================================================

/// Hypercomplex Exponential: exp(a_0 + v) = exp(a_0) * (cos(||v||) + v * (sin(||v||) / ||v||))
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> exp(const nion<T, N, Space>& q) {
    using std::cos;
    using std::exp;
    using std::sin;

    const T a0 = q.real();
    const T v_norm = detail::imag_norm(q.data(), q.size());
    const T exp_a0 = static_cast<T>(exp(a0));

    const T eps = std::numeric_limits<T>::epsilon() * T{16};
    // Numerically stable sinc(v_norm) = sin(v_norm) / v_norm
    const T sinc = (v_norm > eps)
        ? (static_cast<T>(sin(v_norm)) / v_norm)
        : (T{1} - (v_norm * v_norm) / T{6});

    return detail::from_complex_projection(
        q, exp_a0 * static_cast<T>(cos(v_norm)), exp_a0 * sinc);
}

/// Hypercomplex Natural Logarithm: log(q) = log(|q|) + u * atan2(||v||, a_0)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> log(const nion<T, N, Space>& q) {
    using std::log;
    const auto pd = atrig(q);
    nion<T, N, Space> res = pd.versor * pd.angle;
    res[0] = static_cast<T>(log(pd.magnitude));
    return res;
}

/// Hypercomplex Power (nion^scalar): pow(q, p) = exp(p * log(q))
template <concepts::Scalar T, std::size_t N, typename Space, std::convertible_to<T> S>
NION_HOST_DEVICE nion<T, N, Space> pow(const nion<T, N, Space>& q, const S& exponent) {
    const T p = static_cast<T>(exponent);
    if (p == T{0}) {
        nion<T, N, Space> one(q.size(), T{});
        one[0] = T{1};
        return one;
    }
    return exp(log(q) * p);
}

/// Hypercomplex Power (nion^nion): pow(base, exp_nion) = exp(log(base) * exp_nion)
template <concepts::Scalar T, std::size_t N1, std::size_t N2, typename Space>
NION_HOST_DEVICE auto pow(const nion<T, N1, Space>& base, const nion<T, N2, Space>& exp_nion) {
    return exp(log(base) * exp_nion);
}

/// Hypercomplex Square: sqr(q) = q * q
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE constexpr nion<T, N, Space> sqr(const nion<T, N, Space>& q) {
    return q * q;
}

/// Hypercomplex Principal Square Root: sqrt(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> sqrt(const nion<T, N, Space>& q) {
    using std::cos;
    using std::sin;
    using std::sqrt;

    const auto pd = atrig(q);
    const T root_r = static_cast<T>(sqrt(pd.magnitude));
    const T half_theta = pd.angle * T{0.5};

    nion<T, N, Space> res = pd.versor * (root_r * static_cast<T>(sin(half_theta)));
    res[0] = root_r * static_cast<T>(cos(half_theta));
    return res;
}

/// Hypercomplex Principal Cube Root: cbrt(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> cbrt(const nion<T, N, Space>& q) {
    using std::cbrt;
    using std::cos;
    using std::sin;

    const auto pd = atrig(q);
    const T root_r = static_cast<T>(cbrt(pd.magnitude));
    const T third_theta = pd.angle / T{3};

    nion<T, N, Space> res = pd.versor * (root_r * static_cast<T>(sin(third_theta)));
    res[0] = root_r * static_cast<T>(cos(third_theta));
    return res;
}

// ============================================================================
// 2. Hyperbolic Functions (sinh, cosh, tanh, coth, sech, csch & inverses)
// ============================================================================

/// Hypercomplex Hyperbolic Sine:
///   sinh(a_0 + v) = sinh(a_0) cos(||v||) + cosh(a_0) (sin(||v||) / ||v||) v
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> sinh(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T v_norm = detail::imag_norm(q.data(), q.size());
    const T eps = std::numeric_limits<T>::epsilon() * T{16};
    const T sinc = (v_norm > eps)
        ? (static_cast<T>(sin(v_norm)) / v_norm)
        : (T{1} - (v_norm * v_norm) / T{6});

    return detail::from_complex_projection(
        q,
        static_cast<T>(sinh(a0)) * static_cast<T>(cos(v_norm)),
        static_cast<T>(cosh(a0)) * sinc);
}

/// Hypercomplex Hyperbolic Cosine:
///   cosh(a_0 + v) = cosh(a_0) cos(||v||) + sinh(a_0) (sin(||v||) / ||v||) v
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> cosh(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T v_norm = detail::imag_norm(q.data(), q.size());
    const T eps = std::numeric_limits<T>::epsilon() * T{16};
    const T sinc = (v_norm > eps)
        ? (static_cast<T>(sin(v_norm)) / v_norm)
        : (T{1} - (v_norm * v_norm) / T{6});

    return detail::from_complex_projection(
        q,
        static_cast<T>(cosh(a0)) * static_cast<T>(cos(v_norm)),
        static_cast<T>(sinh(a0)) * sinc);
}

/// Hypercomplex Hyperbolic Tangent: tanh(q) = sinh(q) / cosh(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> tanh(const nion<T, N, Space>& q) {
    return sinh(q) / cosh(q);
}

/// Hypercomplex Hyperbolic Cotangent: coth(q) = cosh(q) / sinh(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> coth(const nion<T, N, Space>& q) {
    return cosh(q) / sinh(q);
}

/// Hypercomplex Hyperbolic Secant: sech(q) = 1 / cosh(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> sech(const nion<T, N, Space>& q) {
    return cosh(q).inv();
}

/// Hypercomplex Hyperbolic Cosecant: csch(q) = 1 / sinh(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> csch(const nion<T, N, Space>& q) {
    return sinh(q).inv();
}

/// Inverse Hyperbolic Sine: asinh(q) = log(q + sqrt(q * q + 1))
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> asinh(const nion<T, N, Space>& q) {
    return log(q + sqrt(q * q + T{1}));
}

/// Inverse Hyperbolic Cosine: acosh(q) = log(q + sqrt(q * q - 1))
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> acosh(const nion<T, N, Space>& q) {
    return log(q + sqrt(q * q - T{1}));
}

/// Inverse Hyperbolic Tangent: atanh(q) = 0.5 * (log(1 + q) - log(1 - q))
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> atanh(const nion<T, N, Space>& q) {
    return (log(T{1} + q) - log(T{1} - q)) * T{0.5};
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
NION_HOST_DEVICE nion<T, N, Space> sin(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T v_norm = detail::imag_norm(q.data(), q.size());
    const T eps = std::numeric_limits<T>::epsilon() * T{16};
    const T sinhc = (v_norm > eps)
        ? (static_cast<T>(sinh(v_norm)) / v_norm)
        : (T{1} + (v_norm * v_norm) / T{6});

    return detail::from_complex_projection(
        q,
        static_cast<T>(sin(a0)) * static_cast<T>(cosh(v_norm)),
        static_cast<T>(cos(a0)) * sinhc);
}

/// Hypercomplex Cosine:
///   cos(a_0 + v) = cos(a_0) cosh(||v||) - sin(a_0) (sinh(||v||) / ||v||) v
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> cos(const nion<T, N, Space>& q) {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;

    const T a0 = q.real();
    const T v_norm = detail::imag_norm(q.data(), q.size());
    const T eps = std::numeric_limits<T>::epsilon() * T{16};
    const T sinhc = (v_norm > eps)
        ? (static_cast<T>(sinh(v_norm)) / v_norm)
        : (T{1} + (v_norm * v_norm) / T{6});

    return detail::from_complex_projection(
        q,
        static_cast<T>(cos(a0)) * static_cast<T>(cosh(v_norm)),
        -static_cast<T>(sin(a0)) * sinhc);
}

/// Hypercomplex Tangent: tan(q) = sin(q) / cos(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> tan(const nion<T, N, Space>& q) {
    return sin(q) / cos(q);
}

/// Hypercomplex Cotangent: cot(q) = cos(q) / sin(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> cot(const nion<T, N, Space>& q) {
    return cos(q) / sin(q);
}

/// Hypercomplex Secant: sec(q) = 1 / cos(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> sec(const nion<T, N, Space>& q) {
    return cos(q).inv();
}

/// Hypercomplex Cosecant: csc(q) = 1 / sin(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> csc(const nion<T, N, Space>& q) {
    return sin(q).inv();
}

/// Inverse Sine: asin(q) = -u * asinh(u * q) where u is the unit imaginary versor of q
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> asin(const nion<T, N, Space>& q) {
    const auto pd = atrig(q);
    return (-pd.versor) * asinh(pd.versor * q);
}

/// Inverse Cosine: acos(q) = -u * acosh(q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> acos(const nion<T, N, Space>& q) {
    const auto pd = atrig(q);
    return (-pd.versor) * acosh(q);
}

/// Inverse Tangent: atan(q) = -u * atanh(u * q)
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> atan(const nion<T, N, Space>& q) {
    const auto pd = atrig(q);
    return (-pd.versor) * atanh(pd.versor * q);
}

/// Inverse Cotangent: acot(q) = atan(q.inv())
template <concepts::Scalar T, std::size_t N, typename Space>
NION_HOST_DEVICE nion<T, N, Space> acot(const nion<T, N, Space>& q) {
    return atan(q.inv());
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

} // namespace nion_math

// Expose analytic & trigonometric functions at global scope so unqualified
// calls (`sin(q)`, `exp(q)`, `pow(q, 2)`) work seamlessly out of the box.
using nion_math::PolarDecomposition;
using nion_math::atrig;
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
using nion_math::acot;
using nion_math::asec;
using nion_math::acsc;

#endif // NION_TRIG_HPP
