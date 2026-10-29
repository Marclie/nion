/*
 * Copyright 2023-2026 Marcus Dante Liebenthal
 * Licensed under the Apache License, Version 2.0
 *
 * Tests/benchmark_ci_test.cpp
 * Automated CI/CD test suite with randomized multi-trial validation of
 * numerical accuracy (against std::complex, Boost.Math, Harder's NumberSystems,
 * and mathematical identities) and benchmark speed across N=2, 4, 8, 16, 32.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <complex>
#include <random>
#include <algorithm>
#include <array>
#include <vector>
#include <string>
#include <sstream>

#include "../nion.hpp"
#include "../nion_trig.hpp"

#if __has_include(<boost/math/quaternion.hpp>) && __has_include(<boost/math/octonion.hpp>)
#include <boost/math/quaternion.hpp>
#include <boost/math/octonion.hpp>
#define NION_HAS_BOOST_MATH 1
#else
#define NION_HAS_BOOST_MATH 0
#endif

#if defined(TEST_LARGE_DEGREE) && __has_include("Sedenion.h") && __has_include("Trigintaduonion.h")
#include "Sedenion.h"
#include "Trigintaduonion.h"
#define NION_HAS_HARDER_LIBS 1
#else
#define NION_HAS_HARDER_LIBS 0
#endif

// ============================================================================
// 1. Test Harness & Multi-Trial Statistical Accumulator
// ============================================================================
static int g_tests_run    = 0;
static int g_tests_failed = 0;

/// Prevents the compiler from dead-code eliminating benchmark results
template <typename T>
inline void do_not_optimize(const T& val) noexcept {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(&val) : "memory");
#else
    volatile const char* p = reinterpret_cast<volatile const char*>(&val);
    (void)*p;
#endif
}

/// Euclidean distance ||a - b||_2 between two hypercomplex numbers
template <typename T, std::size_t N>
double l2_diff(const nion<T, N>& a, const nion<T, N>& b) noexcept {
    double sum_sq = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        const double d = static_cast<double>(a[i]) - static_cast<double>(b[i]);
        if (std::isnan(d) || std::isinf(d)) return std::numeric_limits<double>::infinity();
        sum_sq += d * d;
    }
    return std::sqrt(sum_sq);
}

template <typename T, std::size_t N>
double l2_mag(const nion<T, N>& a) noexcept {
    double sum_sq = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        const double v = static_cast<double>(a[i]);
        sum_sq += v * v;
    }
    return std::sqrt(sum_sq);
}

template <typename T, std::size_t N>
std::string format_nion_short(const nion<T, N>& q) {
    std::ostringstream oss;
    oss.precision(4);
    oss << "(";
    const std::size_t limit = (N <= 4) ? N : 4;
    for (std::size_t i = 0; i < limit; ++i) {
        if (i > 0) oss << ", ";
        oss << static_cast<double>(q[i]);
    }
    if (N > limit) oss << ", ...";
    oss << ")";
    return oss.str();
}

/// Tracks MAE, MRE, Max Error, worst-case input, and Nion vs Reference timing across random trials
struct TrialStats {
    int         trials       = 0;
    double      sum_abs_err  = 0.0;
    double      sum_rel_err  = 0.0;
    double      max_abs_err  = 0.0;
    double      max_rel_err  = 0.0;
    double      nion_ns_sum  = 0.0;
    double      ref_ns_sum   = 0.0;
    bool        has_ref_time = false;
    std::string worst_detail;

    void record(double abs_err, double ref_mag, double nion_ns, double ref_ns = -1.0,
                const std::string& detail_on_worst = std::string()) {
        ++trials;
        const double rel_err = (ref_mag > 1e-12) ? (abs_err / ref_mag) : abs_err;
        sum_abs_err += abs_err;
        sum_rel_err += rel_err;
        nion_ns_sum += nion_ns;
        if (ref_ns >= 0.0) {
            ref_ns_sum += ref_ns;
            has_ref_time = true;
        }
        if (!std::isfinite(abs_err) || abs_err > max_abs_err) {
            max_abs_err = abs_err;
        }
        if (!std::isfinite(rel_err) || rel_err > max_rel_err) {
            max_rel_err = rel_err;
            if (!detail_on_worst.empty()) {
                worst_detail = detail_on_worst;
            }
        }
    }

    void assert_and_report(std::size_t degree, const char* op_name, double rel_tol) const {
        ++g_tests_run;
        const double inv_n   = (trials > 0) ? (1.0 / static_cast<double>(trials)) : 0.0;
        const double mae     = sum_abs_err * inv_n;
        const double mre     = sum_rel_err * inv_n;
        const double nion_ns = nion_ns_sum * inv_n;
        const bool   passed  = std::isfinite(max_rel_err) && (max_rel_err <= rel_tol);

        if (!passed) {
            ++g_tests_failed;
        }

        const char* status = passed ? "PASS" : "FAIL";
        if (has_ref_time && ref_ns_sum > 0.0) {
            const double ref_ns  = ref_ns_sum * inv_n;
            const double speedup = ref_ns / (nion_ns > 0.0 ? nion_ns : 1.0);
            std::printf(
                "  [%s] [N=%-2zu] %-28s | trials=%-4d | MAE=%9.2e | MRE=%9.2e | MaxRel=%9.2e (tol %.0e) | nion=%6.1f ns vs ref=%6.1f ns (%5.2fx)\n",
                status, degree, op_name, trials, mae, mre, max_rel_err, rel_tol, nion_ns, ref_ns, speedup);
        } else {
            std::printf(
                "  [%s] [N=%-2zu] %-28s | trials=%-4d | MAE=%9.2e | MRE=%9.2e | MaxRel=%9.2e (tol %.0e) | nion=%6.1f ns\n",
                status, degree, op_name, trials, mae, mre, max_rel_err, rel_tol, nion_ns);
        }

        if (!passed) {
            std::fprintf(stderr,
                "         └── FAILURE DETAIL [N=%zu %s]: MaxRel=%.3e > tol=%.3e (%s)\n",
                degree, op_name, max_rel_err, rel_tol,
                worst_detail.empty() ? "check inputs" : worst_detail.c_str());
        }
    }
};

#define CI_ASSERT_SPEED(elapsed_ns, max_allowed_ns, label)                            \
    do {                                                                              \
        ++g_tests_run;                                                                \
        if ((elapsed_ns) > (max_allowed_ns)) {                                        \
            ++g_tests_failed;                                                         \
            std::fprintf(stderr,                                                      \
                "  [FAIL] %-46s : %7.2f ns/op exceeds limit %.2f ns/op\n",            \
                (label), (elapsed_ns), (max_allowed_ns));                             \
        } else {                                                                      \
            std::printf("  [PASS] %-46s : %7.2f ns/op (limit: %6.2f ns/op)\n",        \
                        (label), (elapsed_ns), (max_allowed_ns));                     \
        }                                                                             \
    } while (0)

/// Measure median nanoseconds per operation across `trials` runs of `iterations` ops
template <typename Fn>
double measure_ns_per_op(Fn&& fn, std::size_t iterations, std::size_t trials = 5) {
    fn(iterations / 4 > 0 ? iterations / 4 : 1);
    std::vector<double> samples(trials);
    for (std::size_t t = 0; t < trials; ++t) {
        const auto t0 = std::chrono::steady_clock::now();
        fn(iterations);
        const auto t1 = std::chrono::steady_clock::now();
        const double total_ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
        samples[t] = total_ns / static_cast<double>(iterations);
    }
    std::sort(samples.begin(), samples.end());
    return samples[trials / 2];
}

// ============================================================================
// 2. Random Hypercomplex Number Generator
// ============================================================================
/// Generates a random nion<T, N> with components in [min_val, max_val] and
/// guaranteed minimum norm >= 0.25 so divisions and inversions are well-conditioned.
template <typename T, std::size_t N>
nion<T, N> random_nion(std::mt19937_64& rng, double min_val = -2.5, double max_val = 2.5) {
    std::uniform_real_distribution<double> dist(min_val, max_val);
    nion<T, N> q;
    double norm_sq = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        double v = dist(rng);
        if (std::abs(v) < 0.05) v = (v >= 0.0 ? 0.15 : -0.15);
        q[i] = static_cast<T>(v);
        norm_sq += v * v;
    }
    if (norm_sq < 0.25) {
        q[0] = static_cast<T>(0.75);
    }
    return q;
}

// ============================================================================
// 3. Generic Mathematical Identity Subsuite (Works for ANY Degree N=2..32)
// ============================================================================
template <typename T, std::size_t N>
void run_algebraic_and_analytic_identities(std::mt19937_64& rng, int num_trials, double tol) {
    TrialStats st_add_sub, st_mul_inv, st_div_mul, st_conj_norm;
    TrialStats st_pythag_trig, st_pythag_hyp, st_sinh_exp, st_cosh_exp;
    TrialStats st_asin_rt, st_acos_rt, st_atan_rt;
    TrialStats st_asinh_rt, st_acosh_rt, st_atanh_rt;

    nion<T, N> unit{};
    unit[0] = T{1};
    nion<T, N> zero{};

    for (int t = 0; t < num_trials; ++t) {
        // Use moderate radius [-1.5, 1.5] so exp/cosh/sinh do not overflow and round-trips stay in principal branch
        const nion<T, N> q1 = random_nion<T, N>(rng, -1.5, 1.5);
        const nion<T, N> q2 = random_nion<T, N>(rng, -1.5, 1.5);
        const double q1_mag = l2_mag(q1);

        // 1. (q1 + q2) - q2 == q1
        {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<T, N> res = (q1 + q2) - q2;
            const auto t1 = std::chrono::steady_clock::now();
            do_not_optimize(res);
            st_add_sub.record(l2_diff(res, q1), q1_mag,
                std::chrono::duration<double, std::nano>(t1 - t0).count());
        }

        // 2. q1 * q1.inv() == 1
        {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<T, N> res = q1 * q1.inv();
            const auto t1 = std::chrono::steady_clock::now();
            do_not_optimize(res);
            st_mul_inv.record(l2_diff(res, unit), 1.0,
                std::chrono::duration<double, std::nano>(t1 - t0).count());
        }

        // 3. (q1 / q2) == q1 * q2.inv()
        {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<T, N> lhs = q1 / q2;
            const auto t1 = std::chrono::steady_clock::now();
            const nion<T, N> rhs = q1 * q2.inv();
            st_div_mul.record(l2_diff(lhs, rhs), l2_mag(rhs),
                std::chrono::duration<double, std::nano>(t1 - t0).count());
        }

        // 4. q1 * conj(q1) == abs(q1) [squared norm ||q1||^2]
        {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<T, N> prod = q1 * q1.conj();
            const T nrm_sq = q1.abs();
            const auto t1 = std::chrono::steady_clock::now();
            nion<T, N> expected_nrm{};
            expected_nrm[0] = nrm_sq;
            st_conj_norm.record(l2_diff(prod, expected_nrm), static_cast<double>(nrm_sq),
                std::chrono::duration<double, std::nano>(t1 - t0).count());
        }

        // 5. sin^2(q1) + cos^2(q1) == 1
        {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<T, N> s = nion_math::sin(q1);
            const nion<T, N> c = nion_math::cos(q1);
            const nion<T, N> res = s * s + c * c;
            const auto t1 = std::chrono::steady_clock::now();
            st_pythag_trig.record(l2_diff(res, unit), 1.0,
                std::chrono::duration<double, std::nano>(t1 - t0).count());
        }

        // 6. cosh^2(q1) - sinh^2(q1) == 1
        {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<T, N> sh = nion_math::sinh(q1);
            const nion<T, N> ch = nion_math::cosh(q1);
            const nion<T, N> res = ch * ch - sh * sh;
            const auto t1 = std::chrono::steady_clock::now();
            st_pythag_hyp.record(l2_diff(res, unit), 1.0,
                std::chrono::duration<double, std::nano>(t1 - t0).count());
        }

        // 7. sinh(q1) == (exp(q1) - exp(-q1)) / 2
        {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<T, N> sh = nion_math::sinh(q1);
            const auto t1 = std::chrono::steady_clock::now();
            const nion<T, N> sh_exp = (nion_math::exp(q1) - nion_math::exp(-q1)) / T{2};
            st_sinh_exp.record(l2_diff(sh, sh_exp), l2_mag(sh_exp),
                std::chrono::duration<double, std::nano>(t1 - t0).count());
        }

        // 8. cosh(q1) == (exp(q1) + exp(-q1)) / 2
        {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<T, N> ch = nion_math::cosh(q1);
            const auto t1 = std::chrono::steady_clock::now();
            const nion<T, N> ch_exp = (nion_math::exp(q1) + nion_math::exp(-q1)) / T{2};
            st_cosh_exp.record(l2_diff(ch, ch_exp), l2_mag(ch_exp),
                std::chrono::duration<double, std::nano>(t1 - t0).count());
        }

        // 9. Round-trip sin(asin(q1)) == q1
        {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<T, N> rt = nion_math::sin(nion_math::asin(q1));
            const auto t1 = std::chrono::steady_clock::now();
            st_asin_rt.record(l2_diff(rt, q1), q1_mag,
                std::chrono::duration<double, std::nano>(t1 - t0).count());
        }

        // 10. Round-trip cos(acos(q1)) == q1
        {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<T, N> rt = nion_math::cos(nion_math::acos(q1));
            const auto t1 = std::chrono::steady_clock::now();
            st_acos_rt.record(l2_diff(rt, q1), q1_mag,
                std::chrono::duration<double, std::nano>(t1 - t0).count());
        }

        // 11. Round-trip tan(atan(q1)) == q1
        {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<T, N> rt = nion_math::tan(nion_math::atan(q1));
            const auto t1 = std::chrono::steady_clock::now();
            st_atan_rt.record(l2_diff(rt, q1), q1_mag,
                std::chrono::duration<double, std::nano>(t1 - t0).count());
        }

        // 12. Round-trip sinh(asinh(q1)) == q1
        {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<T, N> rt = nion_math::sinh(nion_math::asinh(q1));
            const auto t1 = std::chrono::steady_clock::now();
            st_asinh_rt.record(l2_diff(rt, q1), q1_mag,
                std::chrono::duration<double, std::nano>(t1 - t0).count());
        }

        // 13. Round-trip cosh(acosh(q1)) == q1
        {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<T, N> rt = nion_math::cosh(nion_math::acosh(q1));
            const auto t1 = std::chrono::steady_clock::now();
            st_acosh_rt.record(l2_diff(rt, q1), q1_mag,
                std::chrono::duration<double, std::nano>(t1 - t0).count());
        }

        // 14. Round-trip tanh(atanh(q1)) == q1
        {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<T, N> rt = nion_math::tanh(nion_math::atanh(q1));
            const auto t1 = std::chrono::steady_clock::now();
            st_atanh_rt.record(l2_diff(rt, q1), q1_mag,
                std::chrono::duration<double, std::nano>(t1 - t0).count());
        }
    }

    st_add_sub.assert_and_report(N,    "(a + b) - b == a",     tol);
    st_mul_inv.assert_and_report(N,    "a * a.inv() == 1",     tol);
    st_div_mul.assert_and_report(N,    "a / b == a * b.inv()", tol);
    st_conj_norm.assert_and_report(N,  "a * conj(a) == abs(a)",tol);
    st_pythag_trig.assert_and_report(N,"sin^2 + cos^2 == 1",   tol * 10.0);
    st_pythag_hyp.assert_and_report(N, "cosh^2 - sinh^2 == 1", tol * 10.0);
    st_sinh_exp.assert_and_report(N,   "sinh == (e^q-e^-q)/2", tol);
    st_cosh_exp.assert_and_report(N,   "cosh == (e^q+e^-q)/2", tol);
    st_asin_rt.assert_and_report(N,    "sin(asin(q))   == q",  tol);
    st_acos_rt.assert_and_report(N,    "cos(acos(q))   == q",  tol);
    st_atan_rt.assert_and_report(N,    "tan(atan(q))   == q",  tol);
    st_asinh_rt.assert_and_report(N,   "sinh(asinh(q)) == q",  tol);
    st_acosh_rt.assert_and_report(N,   "cosh(acosh(q)) == q",  tol);
    st_atanh_rt.assert_and_report(N,   "tanh(atanh(q)) == q",  tol);
}

// ============================================================================
// 4. Degree N=2 (Complex) vs std::complex<long double>
// ============================================================================
static void test_degree_2_complex(std::mt19937_64& rng, int num_trials) {
    using Real = long double;
    constexpr double kTol = 1e-12;

    std::printf("\n-----------------------------------------------------------------\n");
    std::printf("  Degree N=2 (Complex) — %d Random Trials vs std::complex\n", num_trials);
    std::printf("-----------------------------------------------------------------\n");

    TrialStats st_norm, st_conj, st_add, st_sub, st_mul, st_div, st_pow;
    TrialStats st_exp, st_log, st_sin, st_cos, st_tan, st_asin, st_acos, st_atan;
    TrialStats st_sinh, st_cosh, st_tanh, st_asinh, st_acosh, st_atanh;

    for (int t = 0; t < num_trials; ++t) {
        const nion<Real, 2> z1 = random_nion<Real, 2>(rng, -3.0, 3.0);
        const nion<Real, 2> z2 = random_nion<Real, 2>(rng, -3.0, 3.0);
        const std::complex<Real> c1(z1[0], z1[1]);
        const std::complex<Real> c2(z2[0], z2[1]);

        auto cmp_c = [&](TrialStats& st, auto&& nion_fn, auto&& std_fn) {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<Real, 2> nr = nion_fn();
            const auto t1 = std::chrono::steady_clock::now();
            const auto c_res = std_fn();
            const auto t2 = std::chrono::steady_clock::now();
            do_not_optimize(nr);
            do_not_optimize(c_res);
            const nion<Real, 2> cr({c_res.real(), c_res.imag()});
            st.record(
                l2_diff(nr, cr), l2_mag(cr),
                std::chrono::duration<double, std::nano>(t1 - t0).count(),
                std::chrono::duration<double, std::nano>(t2 - t1).count()
            );
        };

        // Norm (Euclidean magnitude ||z||_2 vs std::abs, and squared norm abs(z) vs std::norm)
        {
            const auto t0 = std::chrono::steady_clock::now();
            const Real nn = z1.norm();
            const auto t1 = std::chrono::steady_clock::now();
            const Real cn = std::abs(c1);
            const auto t2 = std::chrono::steady_clock::now();
            do_not_optimize(nn);
            do_not_optimize(cn);
            st_norm.record(
                static_cast<double>(std::abs(nn - cn)), static_cast<double>(std::abs(cn)),
                std::chrono::duration<double, std::nano>(t1 - t0).count(),
                std::chrono::duration<double, std::nano>(t2 - t1).count()
            );
        }

        cmp_c(st_conj,  [&]() -> nion<Real, 2> { return z1.conj(); },          [&]() { return std::conj(c1); });
        cmp_c(st_add,   [&]() -> nion<Real, 2> { return z1 + z2; },            [&]() { return c1 + c2; });
        cmp_c(st_sub,   [&]() -> nion<Real, 2> { return z1 - z2; },            [&]() { return c1 - c2; });
        cmp_c(st_mul,   [&]() -> nion<Real, 2> { return z1 * z2; },            [&]() { return c1 * c2; });
        cmp_c(st_div,   [&]() -> nion<Real, 2> { return z1 / z2; },            [&]() { return c1 / c2; });
        cmp_c(st_pow,   [&]() -> nion<Real, 2> { return nion_math::pow(z1, z2); }, [&]() { return std::pow(c1, c2); });
        cmp_c(st_exp,   [&]() -> nion<Real, 2> { return nion_math::exp(z1); }, [&]() { return std::exp(c1); });
        cmp_c(st_log,   [&]() -> nion<Real, 2> { return nion_math::log(z1); }, [&]() { return std::log(c1); });
        cmp_c(st_sin,   [&]() -> nion<Real, 2> { return nion_math::sin(z1); }, [&]() { return std::sin(c1); });
        cmp_c(st_cos,   [&]() -> nion<Real, 2> { return nion_math::cos(z1); }, [&]() { return std::cos(c1); });
        cmp_c(st_tan,   [&]() -> nion<Real, 2> { return nion_math::tan(z1); }, [&]() { return std::tan(c1); });
        cmp_c(st_asin,  [&]() -> nion<Real, 2> { return nion_math::asin(z1); },[&]() { return std::asin(c1); });
        cmp_c(st_acos,  [&]() -> nion<Real, 2> { return nion_math::acos(z1); },[&]() { return std::acos(c1); });
        cmp_c(st_atan,  [&]() -> nion<Real, 2> { return nion_math::atan(z1); },[&]() { return std::atan(c1); });
        cmp_c(st_sinh,  [&]() -> nion<Real, 2> { return nion_math::sinh(z1); },[&]() { return std::sinh(c1); });
        cmp_c(st_cosh,  [&]() -> nion<Real, 2> { return nion_math::cosh(z1); },[&]() { return std::cosh(c1); });
        cmp_c(st_tanh,  [&]() -> nion<Real, 2> { return nion_math::tanh(z1); },[&]() { return std::tanh(c1); });
        cmp_c(st_asinh, [&]() -> nion<Real, 2> { return nion_math::asinh(z1);},[&]() { return std::asinh(c1); });
        cmp_c(st_acosh, [&]() -> nion<Real, 2> { return nion_math::acosh(z1);},[&]() { return std::acosh(c1); });
        cmp_c(st_atanh, [&]() -> nion<Real, 2> { return nion_math::atanh(z1);},[&]() { return std::atanh(c1); });
    }

    st_norm.assert_and_report(2,  "norm(z) vs std::norm",       kTol);
    st_conj.assert_and_report(2,  "conj(z) vs std::conj",       kTol);
    st_add.assert_and_report(2,   "add (z1 + z2) vs std",       kTol);
    st_sub.assert_and_report(2,   "sub (z1 - z2) vs std",       kTol);
    st_mul.assert_and_report(2,   "mul (z1 * z2) vs std",       kTol);
    st_div.assert_and_report(2,   "div (z1 / z2) vs std",       kTol);
    st_pow.assert_and_report(2,   "pow(z1, z2) vs std::pow",    kTol * 10.0);
    st_exp.assert_and_report(2,   "exp(z) vs std::exp",         kTol);
    st_log.assert_and_report(2,   "log(z) vs std::log",         kTol);
    st_sin.assert_and_report(2,   "sin(z) vs std::sin",         kTol);
    st_cos.assert_and_report(2,   "cos(z) vs std::cos",         kTol);
    st_tan.assert_and_report(2,   "tan(z) vs std::tan",         kTol);
    st_asin.assert_and_report(2,  "asin(z) vs std::asin",       kTol);
    st_acos.assert_and_report(2,  "acos(z) vs std::acos",       kTol);
    st_atan.assert_and_report(2,  "atan(z) vs std::atan",       kTol);
    st_sinh.assert_and_report(2,  "sinh(z) vs std::sinh",       kTol);
    st_cosh.assert_and_report(2,  "cosh(z) vs std::cosh",       kTol);
    st_tanh.assert_and_report(2,  "tanh(z) vs std::tanh",       kTol);
    st_asinh.assert_and_report(2, "asinh(z) vs std::asinh",     kTol);
    st_acosh.assert_and_report(2, "acosh(z) vs std::acosh",     kTol);
    st_atanh.assert_and_report(2, "atanh(z) vs std::atanh",     kTol);

    run_algebraic_and_analytic_identities<Real, 2>(rng, num_trials, kTol);
}

// ============================================================================
// 5. Degree N=4 (Quaternion) vs boost::math::quaternion & Identities
// ============================================================================
static void test_degree_4_quaternion(std::mt19937_64& rng, int num_trials) {
    using Real = long double;
    constexpr double kTol = 1e-12;

    std::printf("\n-----------------------------------------------------------------\n");
    std::printf("  Degree N=4 (Quaternion) — %d Random Trials vs Boost.Math & Identities\n", num_trials);
    std::printf("-----------------------------------------------------------------\n");

#if NION_HAS_BOOST_MATH
    TrialStats st_norm, st_conj, st_add, st_sub, st_mul, st_div, st_pow;
    TrialStats st_exp, st_sin, st_cos, st_tan, st_sinh, st_cosh, st_tanh;

    for (int t = 0; t < num_trials; ++t) {
        const nion<Real, 4> q1 = random_nion<Real, 4>(rng, -2.5, 2.5);
        const nion<Real, 4> q2 = random_nion<Real, 4>(rng, -2.5, 2.5);
        const int pow_exp = (t % 5) + 1;
        const boost::math::quaternion<Real> b1(q1[0], q1[1], q1[2], q1[3]);
        const boost::math::quaternion<Real> b2(q2[0], q2[1], q2[2], q2[3]);

        auto cmp_q = [&](TrialStats& st, auto&& nion_fn, auto&& boost_fn) {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<Real, 4> nr = nion_fn();
            const auto t1 = std::chrono::steady_clock::now();
            const auto b_res = boost_fn();
            const auto t2 = std::chrono::steady_clock::now();
            do_not_optimize(nr);
            do_not_optimize(b_res);
            const nion<Real, 4> br({
                b_res.R_component_1(), b_res.R_component_2(),
                b_res.R_component_3(), b_res.R_component_4()
            });
            st.record(
                l2_diff(nr, br), l2_mag(br),
                std::chrono::duration<double, std::nano>(t1 - t0).count(),
                std::chrono::duration<double, std::nano>(t2 - t1).count()
            );
        };

        {
            const auto t0 = std::chrono::steady_clock::now();
            const Real nn = q1.abs();
            const auto t1 = std::chrono::steady_clock::now();
            const Real bn = boost::math::norm(b1);
            const auto t2 = std::chrono::steady_clock::now();
            do_not_optimize(nn);
            do_not_optimize(bn);
            st_norm.record(
                static_cast<double>(std::abs(nn - bn)), static_cast<double>(std::abs(bn)),
                std::chrono::duration<double, std::nano>(t1 - t0).count(),
                std::chrono::duration<double, std::nano>(t2 - t1).count()
            );
        }

        cmp_q(st_conj, [&]() -> nion<Real, 4> { return q1.conj(); },                [&]() { return boost::math::conj(b1); });
        cmp_q(st_add,  [&]() -> nion<Real, 4> { return q1 + q2; },                  [&]() { return b1 + b2; });
        cmp_q(st_sub,  [&]() -> nion<Real, 4> { return q1 - q2; },                  [&]() { return b1 - b2; });
        cmp_q(st_mul,  [&]() -> nion<Real, 4> { return q1 * q2; },                  [&]() { return b1 * b2; });
        cmp_q(st_div,  [&]() -> nion<Real, 4> { return q1 / q2; },                  [&]() { return b1 / b2; });
        cmp_q(st_pow,  [&]() -> nion<Real, 4> { return nion_math::pow(q1, pow_exp);},[&]() { return boost::math::pow(b1, pow_exp); });
        cmp_q(st_exp,  [&]() -> nion<Real, 4> { return nion_math::exp(q1); },       [&]() { return boost::math::exp(b1); });
        cmp_q(st_sin,  [&]() -> nion<Real, 4> { return nion_math::sin(q1); },       [&]() { return boost::math::sin(b1); });
        cmp_q(st_cos,  [&]() -> nion<Real, 4> { return nion_math::cos(q1); },       [&]() { return boost::math::cos(b1); });
        cmp_q(st_tan,  [&]() -> nion<Real, 4> { return nion_math::tan(q1); },       [&]() { return boost::math::tan(b1); });
        cmp_q(st_sinh, [&]() -> nion<Real, 4> { return nion_math::sinh(q1); },      [&]() { return boost::math::sinh(b1); });
        cmp_q(st_cosh, [&]() -> nion<Real, 4> { return nion_math::cosh(q1); },      [&]() { return boost::math::cosh(b1); });
        cmp_q(st_tanh, [&]() -> nion<Real, 4> { return nion_math::tanh(q1); },      [&]() { return boost::math::tanh(b1); });
    }

    st_norm.assert_and_report(4, "norm(q) vs Boost",        kTol);
    st_conj.assert_and_report(4, "conj(q) vs Boost",        kTol);
    st_add.assert_and_report(4,  "add (q1 + q2) vs Boost",  kTol);
    st_sub.assert_and_report(4,  "sub (q1 - q2) vs Boost",  kTol);
    st_mul.assert_and_report(4,  "mul (q1 * q2) vs Boost",  kTol);
    st_div.assert_and_report(4,  "div (q1 / q2) vs Boost",  kTol);
    st_pow.assert_and_report(4,  "pow(q, n) vs Boost",      kTol);
    st_exp.assert_and_report(4,  "exp(q) vs Boost",         kTol);
    st_sin.assert_and_report(4,  "sin(q) vs Boost",         kTol);
    st_cos.assert_and_report(4,  "cos(q) vs Boost",         kTol);
    st_tan.assert_and_report(4,  "tan(q) vs Boost",         kTol);
    st_sinh.assert_and_report(4, "sinh(q) vs Boost",        kTol);
    st_cosh.assert_and_report(4, "cosh(q) vs Boost",        kTol);
    st_tanh.assert_and_report(4, "tanh(q) vs Boost",        kTol);
#endif

    run_algebraic_and_analytic_identities<Real, 4>(rng, num_trials, kTol);
}

// ============================================================================
// 6. Degree N=8 (Octonion) vs boost::math::octonion & Identities
// ============================================================================
static void test_degree_8_octonion(std::mt19937_64& rng, int num_trials) {
    using Real = long double;
    constexpr double kTol = 1e-12;

    std::printf("\n-----------------------------------------------------------------\n");
    std::printf("  Degree N=8 (Octonion) — %d Random Trials vs Boost.Math & Identities\n", num_trials);
    std::printf("-----------------------------------------------------------------\n");

#if NION_HAS_BOOST_MATH
    TrialStats st_norm, st_conj, st_add, st_sub, st_mul, st_div, st_pow;
    TrialStats st_exp, st_sin, st_cos, st_tan;

    for (int t = 0; t < num_trials; ++t) {
        const nion<Real, 8> o1 = random_nion<Real, 8>(rng, -2.5, 2.5);
        const nion<Real, 8> o2 = random_nion<Real, 8>(rng, -2.5, 2.5);
        const int pow_exp = (t % 5) + 1;
        const boost::math::octonion<Real> b1(o1[0], o1[1], o1[2], o1[3], o1[4], o1[5], o1[6], o1[7]);
        const boost::math::octonion<Real> b2(o2[0], o2[1], o2[2], o2[3], o2[4], o2[5], o2[6], o2[7]);

        auto cmp_o = [&](TrialStats& st, auto&& nion_fn, auto&& boost_fn) {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<Real, 8> nr = nion_fn();
            const auto t1 = std::chrono::steady_clock::now();
            const auto b_res = boost_fn();
            const auto t2 = std::chrono::steady_clock::now();
            do_not_optimize(nr);
            do_not_optimize(b_res);
            const nion<Real, 8> br({
                b_res.R_component_1(), b_res.R_component_2(), b_res.R_component_3(), b_res.R_component_4(),
                b_res.R_component_5(), b_res.R_component_6(), b_res.R_component_7(), b_res.R_component_8()
            });
            st.record(
                l2_diff(nr, br), l2_mag(br),
                std::chrono::duration<double, std::nano>(t1 - t0).count(),
                std::chrono::duration<double, std::nano>(t2 - t1).count()
            );
        };

        {
            const auto t0 = std::chrono::steady_clock::now();
            const Real nn = o1.abs();
            const auto t1 = std::chrono::steady_clock::now();
            const Real bn = boost::math::norm(b1);
            const auto t2 = std::chrono::steady_clock::now();
            do_not_optimize(nn);
            do_not_optimize(bn);
            st_norm.record(
                static_cast<double>(std::abs(nn - bn)), static_cast<double>(std::abs(bn)),
                std::chrono::duration<double, std::nano>(t1 - t0).count(),
                std::chrono::duration<double, std::nano>(t2 - t1).count()
            );
        }

        cmp_o(st_conj, [&]() -> nion<Real, 8> { return o1.conj(); },                 [&]() { return boost::math::conj(b1); });
        cmp_o(st_add,  [&]() -> nion<Real, 8> { return o1 + o2; },                   [&]() { return b1 + b2; });
        cmp_o(st_sub,  [&]() -> nion<Real, 8> { return o1 - o2; },                   [&]() { return b1 - b2; });
        cmp_o(st_mul,  [&]() -> nion<Real, 8> { return o1 * o2; },                   [&]() { return b1 * b2; });
        cmp_o(st_div,  [&]() -> nion<Real, 8> { return o1 / o2; },                   [&]() { return b1 / b2; });
        cmp_o(st_pow,  [&]() -> nion<Real, 8> { return nion_math::pow(o1, pow_exp); },[&]() { return boost::math::pow(b1, pow_exp); });
        cmp_o(st_exp,  [&]() -> nion<Real, 8> { return nion_math::exp(o1); },        [&]() { return boost::math::exp(b1); });
        cmp_o(st_sin,  [&]() -> nion<Real, 8> { return nion_math::sin(o1); },        [&]() { return boost::math::sin(b1); });
        cmp_o(st_cos,  [&]() -> nion<Real, 8> { return nion_math::cos(o1); },        [&]() { return boost::math::cos(b1); });
        cmp_o(st_tan,  [&]() -> nion<Real, 8> { return nion_math::tan(o1); },        [&]() { return boost::math::tan(b1); });
    }

    st_norm.assert_and_report(8, "norm(o) vs Boost",        kTol);
    st_conj.assert_and_report(8, "conj(o) vs Boost",        kTol);
    st_add.assert_and_report(8,  "add (o1 + o2) vs Boost",  kTol);
    st_sub.assert_and_report(8,  "sub (o1 - o2) vs Boost",  kTol);
    st_mul.assert_and_report(8,  "mul (o1 * o2) vs Boost",  kTol);
    st_div.assert_and_report(8,  "div (o1 / o2) vs Boost",  kTol);
    st_pow.assert_and_report(8,  "pow(o, n) vs Boost",      kTol);
    st_exp.assert_and_report(8,  "exp(o) vs Boost",         kTol);
    st_sin.assert_and_report(8,  "sin(o) vs Boost",         kTol);
    st_cos.assert_and_report(8,  "cos(o) vs Boost",         kTol);
    st_tan.assert_and_report(8,  "tan(o) vs Boost",         kTol);
#endif

    run_algebraic_and_analytic_identities<Real, 8>(rng, num_trials, kTol);
}

// ============================================================================
// 7. Degree N=16 (Sedenion) & N=32 (Trigintaduonion)
// ============================================================================
static void test_degree_16_sedenion(std::mt19937_64& rng, int num_trials) {
    using Real = double;
    constexpr double kTol = 1e-11;

    std::printf("\n-----------------------------------------------------------------\n");
#if NION_HAS_HARDER_LIBS
    std::printf("  Degree N=16 (Sedenion) — %d Random Trials vs Harder's Sedenion & Identities\n", num_trials);
#else
    std::printf("  Degree N=16 (Sedenion) — %d Random Trials (Algebraic & Analytic Identities)\n", num_trials);
#endif
    std::printf("-----------------------------------------------------------------\n");

#if NION_HAS_HARDER_LIBS
    TrialStats st_abs, st_conj, st_add, st_sub, st_mul, st_div, st_inv;
    TrialStats st_exp, st_log, st_sin, st_cos, st_tan, st_asin, st_acos, st_atan;
    TrialStats st_sinh, st_cosh, st_tanh, st_asinh, st_acosh, st_atanh;

    for (int t = 0; t < num_trials; ++t) {
        const nion<Real, 16> s1 = random_nion<Real, 16>(rng, -2.0, 2.0);
        const nion<Real, 16> s2 = random_nion<Real, 16>(rng, -2.0, 2.0);
        const Sedenion<Real> h1(
            s1[0], s1[1], s1[2], s1[3], s1[4], s1[5], s1[6], s1[7],
            s1[8], s1[9], s1[10], s1[11], s1[12], s1[13], s1[14], s1[15]
        );
        const Sedenion<Real> h2(
            s2[0], s2[1], s2[2], s2[3], s2[4], s2[5], s2[6], s2[7],
            s2[8], s2[9], s2[10], s2[11], s2[12], s2[13], s2[14], s2[15]
        );

        auto to_nion16 = [](const Sedenion<Real>& h) -> nion<Real, 16> {
            return nion<Real, 16>({
                h.real(), h.imag_i(), h.imag_j(), h.imag_k(),
                h.imag_u1(), h.imag_i1(), h.imag_j1(), h.imag_k1(),
                h.imag_u2(), h.imag_i2(), h.imag_j2(), h.imag_k2(),
                h.imag_u3(), h.imag_i3(), h.imag_j3(), h.imag_k3()
            });
        };

        auto cmp_s = [&](TrialStats& st, auto&& nion_fn, auto&& harder_fn) {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<Real, 16> nr = nion_fn();
            const auto t1 = std::chrono::steady_clock::now();
            const Sedenion<Real> hr_raw = harder_fn();
            const auto t2 = std::chrono::steady_clock::now();
            do_not_optimize(nr);
            do_not_optimize(hr_raw);
            const nion<Real, 16> hr = to_nion16(hr_raw);
            st.record(
                l2_diff(nr, hr), l2_mag(hr),
                std::chrono::duration<double, std::nano>(t1 - t0).count(),
                std::chrono::duration<double, std::nano>(t2 - t1).count()
            );
        };

        {
            const auto t0 = std::chrono::steady_clock::now();
            const Real na = s1.norm();
            const auto t1 = std::chrono::steady_clock::now();
            const Real ha = h1.abs();
            const auto t2 = std::chrono::steady_clock::now();
            do_not_optimize(na);
            do_not_optimize(ha);
            st_abs.record(
                std::abs(na - ha), std::abs(ha),
                std::chrono::duration<double, std::nano>(t1 - t0).count(),
                std::chrono::duration<double, std::nano>(t2 - t1).count()
            );
        }

        cmp_s(st_conj,  [&]() -> nion<Real, 16> { return s1.conj(); },            [&]() { return h1.conj(); });
        cmp_s(st_add,   [&]() -> nion<Real, 16> { return s1 + s2; },              [&]() { return h1 + h2; });
        cmp_s(st_sub,   [&]() -> nion<Real, 16> { return s1 - s2; },              [&]() { return h1 - h2; });
        cmp_s(st_mul,   [&]() -> nion<Real, 16> { return s1 * s2; },              [&]() { return h1 * h2; });
        cmp_s(st_div,   [&]() -> nion<Real, 16> { return s1 / s2; },              [&]() { return h1 / h2; });
        cmp_s(st_inv,   [&]() -> nion<Real, 16> { return s1.inv(); },             [&]() { return 1.0 / h1; });
        cmp_s(st_exp,   [&]() -> nion<Real, 16> { return nion_math::exp(s1); },   [&]() { return h1.exp(); });
        cmp_s(st_log,   [&]() -> nion<Real, 16> { return nion_math::log(s1); },   [&]() { return h1.log(); });
        cmp_s(st_sin,   [&]() -> nion<Real, 16> { return nion_math::sin(s1); },   [&]() { return h1.sin(); });
        cmp_s(st_cos,   [&]() -> nion<Real, 16> { return nion_math::cos(s1); },   [&]() { return h1.cos(); });
        cmp_s(st_tan,   [&]() -> nion<Real, 16> { return nion_math::tan(s1); },   [&]() { return h1.tan(); });
        cmp_s(st_asin,  [&]() -> nion<Real, 16> { return nion_math::asin(s1); },  [&]() { return h1.asin(); });
        cmp_s(st_acos,  [&]() -> nion<Real, 16> { return nion_math::acos(s1); },  [&]() { return h1.acos(); });
        cmp_s(st_atan,  [&]() -> nion<Real, 16> { return nion_math::atan(s1); },  [&]() { return h1.atan(); });
        cmp_s(st_sinh,  [&]() -> nion<Real, 16> { return nion_math::sinh(s1); },  [&]() { return h1.sinh(); });
        cmp_s(st_cosh,  [&]() -> nion<Real, 16> { return nion_math::cosh(s1); },  [&]() { return h1.cosh(); });
        cmp_s(st_tanh,  [&]() -> nion<Real, 16> { return nion_math::tanh(s1); },  [&]() { return h1.tanh(); });
        cmp_s(st_asinh, [&]() -> nion<Real, 16> { return nion_math::asinh(s1); }, [&]() { return h1.asinh(); });
        cmp_s(st_acosh, [&]() -> nion<Real, 16> { return nion_math::acosh(s1); }, [&]() { return h1.acosh(); });
        cmp_s(st_atanh, [&]() -> nion<Real, 16> { return nion_math::atanh(s1); }, [&]() { return h1.atanh(); });
    }

    st_abs.assert_and_report(16,   "abs(s) vs Sedenion",        kTol);
    st_conj.assert_and_report(16,  "conj(s) vs Sedenion",       kTol);
    st_add.assert_and_report(16,   "add (s1 + s2) vs Sedenion", kTol);
    st_sub.assert_and_report(16,   "sub (s1 - s2) vs Sedenion", kTol);
    st_mul.assert_and_report(16,   "mul (s1 * s2) vs Sedenion", kTol);
    st_div.assert_and_report(16,   "div (s1 / s2) vs Sedenion", kTol);
    st_inv.assert_and_report(16,   "inv(s) vs Sedenion",        kTol);
    st_exp.assert_and_report(16,   "exp(s) vs Sedenion",        kTol);
    st_log.assert_and_report(16,   "log(s) vs Sedenion",        kTol);
    st_sin.assert_and_report(16,   "sin(s) vs Sedenion",        kTol);
    st_cos.assert_and_report(16,   "cos(s) vs Sedenion",        kTol);
    st_tan.assert_and_report(16,   "tan(s) vs Sedenion",        kTol);
    st_asin.assert_and_report(16,  "asin(s) vs Sedenion",       kTol);
    st_acos.assert_and_report(16,  "acos(s) vs Sedenion",       kTol);
    st_atan.assert_and_report(16,  "atan(s) vs Sedenion",       kTol);
    st_sinh.assert_and_report(16,  "sinh(s) vs Sedenion",       kTol);
    st_cosh.assert_and_report(16,  "cosh(s) vs Sedenion",       kTol);
    st_tanh.assert_and_report(16,  "tanh(s) vs Sedenion",       kTol);
    st_asinh.assert_and_report(16, "asinh(s) vs Sedenion",      kTol);
    st_acosh.assert_and_report(16, "acosh(s) vs Sedenion",      kTol);
    st_atanh.assert_and_report(16, "atanh(s) vs Sedenion",      kTol);
#endif

    run_algebraic_and_analytic_identities<Real, 16>(rng, num_trials, kTol);
}

static void test_degree_32_trigintaduonion(std::mt19937_64& rng, int num_trials) {
    using Real = double;
    constexpr double kTol = 1e-11;

    std::printf("\n-----------------------------------------------------------------\n");
#if NION_HAS_HARDER_LIBS
    std::printf("  Degree N=32 (Trigintaduonion) — %d Random Trials vs Harder's Library & Identities\n", num_trials);
#else
    std::printf("  Degree N=32 (Trigintaduonion) — %d Random Trials (Algebraic & Analytic Identities)\n", num_trials);
#endif
    std::printf("-----------------------------------------------------------------\n");

#if NION_HAS_HARDER_LIBS
    TrialStats st_abs, st_conj, st_add, st_sub, st_mul, st_div, st_inv;
    TrialStats st_exp, st_log, st_sin, st_cos, st_tan, st_asin, st_acos, st_atan;
    TrialStats st_sinh, st_cosh, st_tanh, st_asinh, st_acosh, st_atanh;

    for (int t = 0; t < num_trials; ++t) {
        const nion<Real, 32> q1 = random_nion<Real, 32>(rng, -1.5, 1.5);
        const nion<Real, 32> q2 = random_nion<Real, 32>(rng, -1.5, 1.5);
        const Trigintaduonion<Real> h1(
            q1[0],  q1[1],  q1[2],  q1[3],  q1[4],  q1[5],  q1[6],  q1[7],
            q1[8],  q1[9],  q1[10], q1[11], q1[12], q1[13], q1[14], q1[15],
            q1[16], q1[17], q1[18], q1[19], q1[20], q1[21], q1[22], q1[23],
            q1[24], q1[25], q1[26], q1[27], q1[28], q1[29], q1[30], q1[31]
        );
        const Trigintaduonion<Real> h2(
            q2[0],  q2[1],  q2[2],  q2[3],  q2[4],  q2[5],  q2[6],  q2[7],
            q2[8],  q2[9],  q2[10], q2[11], q2[12], q2[13], q2[14], q2[15],
            q2[16], q2[17], q2[18], q2[19], q2[20], q2[21], q2[22], q2[23],
            q2[24], q2[25], q2[26], q2[27], q2[28], q2[29], q2[30], q2[31]
        );

        auto to_nion32 = [](const Trigintaduonion<Real>& h) -> nion<Real, 32> {
            return nion<Real, 32>({
                h.real(),    h.imag_i(),  h.imag_j(),  h.imag_k(),
                h.imag_u1(), h.imag_i1(), h.imag_j1(), h.imag_k1(),
                h.imag_u2(), h.imag_i2(), h.imag_j2(), h.imag_k2(),
                h.imag_u3(), h.imag_i3(), h.imag_j3(), h.imag_k3(),
                h.imag_u4(), h.imag_i4(), h.imag_j4(), h.imag_k4(),
                h.imag_u5(), h.imag_i5(), h.imag_j5(), h.imag_k5(),
                h.imag_u6(), h.imag_i6(), h.imag_j6(), h.imag_k6(),
                h.imag_u7(), h.imag_i7(), h.imag_j7(), h.imag_k7()
            });
        };

        auto cmp_t = [&](TrialStats& st, auto&& nion_fn, auto&& harder_fn) {
            const auto t0 = std::chrono::steady_clock::now();
            const nion<Real, 32> nr = nion_fn();
            const auto t1 = std::chrono::steady_clock::now();
            const Trigintaduonion<Real> hr_raw = harder_fn();
            const auto t2 = std::chrono::steady_clock::now();
            do_not_optimize(nr);
            do_not_optimize(hr_raw);
            const nion<Real, 32> hr = to_nion32(hr_raw);
            st.record(
                l2_diff(nr, hr), l2_mag(hr),
                std::chrono::duration<double, std::nano>(t1 - t0).count(),
                std::chrono::duration<double, std::nano>(t2 - t1).count()
            );
        };

        {
            const auto t0 = std::chrono::steady_clock::now();
            const Real na = q1.norm();
            const auto t1 = std::chrono::steady_clock::now();
            const Real ha = h1.abs();
            const auto t2 = std::chrono::steady_clock::now();
            do_not_optimize(na);
            do_not_optimize(ha);
            st_abs.record(
                std::abs(na - ha), std::abs(ha),
                std::chrono::duration<double, std::nano>(t1 - t0).count(),
                std::chrono::duration<double, std::nano>(t2 - t1).count()
            );
        }

        cmp_t(st_conj,  [&]() -> nion<Real, 32> { return q1.conj(); },            [&]() { return h1.conj(); });
        cmp_t(st_add,   [&]() -> nion<Real, 32> { return q1 + q2; },              [&]() { return h1 + h2; });
        cmp_t(st_sub,   [&]() -> nion<Real, 32> { return q1 - q2; },              [&]() { return h1 - h2; });
        cmp_t(st_mul,   [&]() -> nion<Real, 32> { return q1 * q2; },              [&]() { return h1 * h2; });
        cmp_t(st_div,   [&]() -> nion<Real, 32> { return q1 / q2; },              [&]() { return h1 / h2; });
        cmp_t(st_inv,   [&]() -> nion<Real, 32> { return q1.inv(); },             [&]() { return 1.0 / h1; });
        cmp_t(st_exp,   [&]() -> nion<Real, 32> { return nion_math::exp(q1); },   [&]() { return h1.exp(); });
        cmp_t(st_log,   [&]() -> nion<Real, 32> { return nion_math::log(q1); },   [&]() { return h1.log(); });
        cmp_t(st_sin,   [&]() -> nion<Real, 32> { return nion_math::sin(q1); },   [&]() { return h1.sin(); });
        cmp_t(st_cos,   [&]() -> nion<Real, 32> { return nion_math::cos(q1); },   [&]() { return h1.cos(); });
        cmp_t(st_tan,   [&]() -> nion<Real, 32> { return nion_math::tan(q1); },   [&]() { return h1.tan(); });
        cmp_t(st_asin,  [&]() -> nion<Real, 32> { return nion_math::asin(q1); },  [&]() { return h1.asin(); });
        cmp_t(st_acos,  [&]() -> nion<Real, 32> { return nion_math::acos(q1); },  [&]() { return h1.acos(); });
        cmp_t(st_atan,  [&]() -> nion<Real, 32> { return nion_math::atan(q1); },  [&]() { return h1.atan(); });
        cmp_t(st_sinh,  [&]() -> nion<Real, 32> { return nion_math::sinh(q1); },  [&]() { return h1.sinh(); });
        cmp_t(st_cosh,  [&]() -> nion<Real, 32> { return nion_math::cosh(q1); },  [&]() { return h1.cosh(); });
        cmp_t(st_tanh,  [&]() -> nion<Real, 32> { return nion_math::tanh(q1); },  [&]() { return h1.tanh(); });
        cmp_t(st_asinh, [&]() -> nion<Real, 32> { return nion_math::asinh(q1); }, [&]() { return h1.asinh(); });
        cmp_t(st_acosh, [&]() -> nion<Real, 32> { return nion_math::acosh(q1); }, [&]() { return h1.acosh(); });
        cmp_t(st_atanh, [&]() -> nion<Real, 32> { return nion_math::atanh(q1); }, [&]() { return h1.atanh(); });
    }

    st_abs.assert_and_report(32,   "abs(t) vs Trigintaduonion",  kTol);
    st_conj.assert_and_report(32,  "conj(t) vs Trigintaduonion", kTol);
    st_add.assert_and_report(32,   "add vs Trigintaduonion",     kTol);
    st_sub.assert_and_report(32,   "sub vs Trigintaduonion",     kTol);
    st_mul.assert_and_report(32,   "mul vs Trigintaduonion",     kTol);
    st_div.assert_and_report(32,   "div vs Trigintaduonion",     kTol);
    st_inv.assert_and_report(32,   "inv vs Trigintaduonion",     kTol);
    st_exp.assert_and_report(32,   "exp vs Trigintaduonion",     kTol);
    st_log.assert_and_report(32,   "log vs Trigintaduonion",     kTol);
    st_sin.assert_and_report(32,   "sin vs Trigintaduonion",     kTol);
    st_cos.assert_and_report(32,   "cos vs Trigintaduonion",     kTol);
    st_tan.assert_and_report(32,   "tan vs Trigintaduonion",     kTol);
    st_asin.assert_and_report(32,  "asin vs Trigintaduonion",    kTol);
    st_acos.assert_and_report(32,  "acos vs Trigintaduonion",    kTol);
    st_atan.assert_and_report(32,  "atan vs Trigintaduonion",    kTol);
    st_sinh.assert_and_report(32,  "sinh vs Trigintaduonion",    kTol);
    st_cosh.assert_and_report(32,  "cosh vs Trigintaduonion",    kTol);
    st_tanh.assert_and_report(32,  "tanh vs Trigintaduonion",    kTol);
    st_asinh.assert_and_report(32, "asinh vs Trigintaduonion",   kTol);
    st_acosh.assert_and_report(32, "acosh vs Trigintaduonion",   kTol);
    st_atanh.assert_and_report(32, "atanh vs Trigintaduonion",   kTol);
#endif

    run_algebraic_and_analytic_identities<Real, 32>(rng, num_trials, kTol);
}

// ============================================================================
// 8. SoA SIMD Batch Multiplication Accuracy (Random Trials Across N=4, 8, 16)
// ============================================================================
template <std::size_t N>
static void run_soa_batch_random_trial(std::mt19937_64& rng, std::size_t batch_size) {
    TrialStats st_soa;
    nion_math::nion_batch_soa<double, N> batch_a(batch_size);
    nion_math::nion_batch_soa<double, N> batch_b(batch_size);
    nion_math::nion_batch_soa<double, N> batch_out(batch_size);

    for (std::size_t idx = 0; idx < batch_size; ++idx) {
        batch_a.set(idx, random_nion<double, N>(rng, -3.0, 3.0));
        batch_b.set(idx, random_nion<double, N>(rng, -3.0, 3.0));
    }

    const auto t0 = std::chrono::steady_clock::now();
    nion_math::nion_batch_soa<double, N>::multiply_batch(batch_out, batch_a, batch_b);
    const auto t1 = std::chrono::steady_clock::now();
    const double soa_ns_per_elem =
        std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(batch_size);

    for (std::size_t idx = 0; idx < batch_size; ++idx) {
        const auto s0 = std::chrono::steady_clock::now();
        const nion<double, N> expected = batch_a.get(idx) * batch_b.get(idx);
        const auto s1 = std::chrono::steady_clock::now();
        const double aos_ns = std::chrono::duration<double, std::nano>(s1 - s0).count();
        const nion<double, N> actual = batch_out.get(idx);
        st_soa.record(l2_diff(actual, expected), l2_mag(expected), soa_ns_per_elem, aos_ns);
    }

    st_soa.assert_and_report(N, "SoA SIMD batch mul vs AoS", 1e-12);
}

static void test_soa_simd_batches(std::mt19937_64& rng, int num_trials) {
    const std::size_t batch_size = static_cast<std::size_t>(std::max(64, num_trials));
    std::printf("\n-----------------------------------------------------------------\n");
    std::printf("  SoA SIMD Batch Multiplication — %zu Random Elements per Degree\n", batch_size);
    std::printf("-----------------------------------------------------------------\n");
    run_soa_batch_random_trial<4>(rng, batch_size);
    run_soa_batch_random_trial<8>(rng, batch_size);
    run_soa_batch_random_trial<16>(rng, batch_size);
}

// ============================================================================
// 9. Performance (Speed) Benchmark Assertions
// ============================================================================
static void test_benchmark_performance(std::mt19937_64& rng) {
    double ci_scale = 1.0;
    if (const char* env_scale = std::getenv("NION_CI_SPEED_SCALE")) {
        const double parsed = std::atof(env_scale);
        if (parsed > 0.0) ci_scale = parsed;
    }

    std::printf("\n-----------------------------------------------------------------\n");
    std::printf("  High-Throughput Performance Benchmark Assertions (scale=%.1fx)\n", ci_scale);
    std::printf("-----------------------------------------------------------------\n");

    constexpr std::size_t kItersFast  = 200000;
    constexpr std::size_t kItersLarge = 20000;

    // 1. Complex (N=2) multiplication speed (unit-normalized multiplier prevents subnormal underflow)
    {
        nion<double, 2> a = random_nion<double, 2>(rng, -1.5, 1.5);
        nion<double, 2> b = random_nion<double, 2>(rng, -1.5, 1.5);
        b = b / b.norm();
        const double ns = measure_ns_per_op([&](std::size_t n) {
            nion<double, 2> acc = a;
            for (std::size_t i = 0; i < n; ++i) {
                acc = acc * b;
                do_not_optimize(acc);
            }
        }, kItersFast);
        CI_ASSERT_SPEED(ns, 25.0 * ci_scale, "[N=2 ] Complex multiply latency");
    }

    // 2. Quaternion (N=4) multiplication speed
    {
        nion<double, 4> a = random_nion<double, 4>(rng, -1.5, 1.5);
        nion<double, 4> b = random_nion<double, 4>(rng, -1.5, 1.5);
        b = b / b.norm();
        const double ns = measure_ns_per_op([&](std::size_t n) {
            nion<double, 4> acc = a;
            for (std::size_t i = 0; i < n; ++i) {
                acc = acc * b;
                do_not_optimize(acc);
            }
        }, kItersFast);
        CI_ASSERT_SPEED(ns, 40.0 * ci_scale, "[N=4 ] Quaternion multiply latency");
    }

    // 3. Octonion (N=8) fast-path multiplication speed (and comparison vs Boost if available)
    {
        nion<double, 8> a = random_nion<double, 8>(rng, -1.5, 1.5);
        nion<double, 8> b = random_nion<double, 8>(rng, -1.5, 1.5);
        b = b / b.norm();
        const double nion_ns = measure_ns_per_op([&](std::size_t n) {
            nion<double, 8> acc = a;
            for (std::size_t i = 0; i < n; ++i) {
                acc = acc * b;
                do_not_optimize(acc);
            }
        }, kItersFast);
        CI_ASSERT_SPEED(nion_ns, 80.0 * ci_scale, "[N=8 ] Octonion fast-path multiply latency");

#if NION_HAS_BOOST_MATH
        boost::math::octonion<double> ba(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]);
        boost::math::octonion<double> bb(b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]);
        const double boost_ns = measure_ns_per_op([&](std::size_t n) {
            boost::math::octonion<double> acc = ba;
            for (std::size_t i = 0; i < n; ++i) {
                acc = acc * bb;
                do_not_optimize(acc);
            }
        }, kItersFast);
        CI_ASSERT_SPEED(nion_ns, boost_ns * 1.5 * ci_scale, "[N=8 ] Octonion multiply vs Boost baseline");
#endif
    }

    // 4. Sedenion (N=16) fast-path multiplication speed
    {
        nion<double, 16> a = random_nion<double, 16>(rng, -1.5, 1.5);
        nion<double, 16> b = random_nion<double, 16>(rng, -1.5, 1.5);
        b = b / b.norm();
        const double ns = measure_ns_per_op([&](std::size_t n) {
            nion<double, 16> acc = a;
            for (std::size_t i = 0; i < n; ++i) {
                acc = acc * b;
                do_not_optimize(acc);
            }
        }, kItersLarge);
        CI_ASSERT_SPEED(ns, 250.0 * ci_scale, "[N=16] Sedenion fast-path multiply latency");
    }

    // 5. SoA SIMD Octonion Batch (N=8, M=4096) per-element throughput
    {
        constexpr std::size_t kBatchSize = 4096;
        nion_math::nion_batch_soa<double, 8> ba(kBatchSize);
        nion_math::nion_batch_soa<double, 8> bb(kBatchSize);
        nion_math::nion_batch_soa<double, 8> bc(kBatchSize);

        const double batch_ns = measure_ns_per_op([&](std::size_t n) {
            for (std::size_t i = 0; i < n; ++i) {
                nion_math::nion_batch_soa<double, 8>::multiply_batch(bc, ba, bb);
                do_not_optimize(bc);
            }
        }, 50);
        const double ns_per_octonion = batch_ns / static_cast<double>(kBatchSize);
        CI_ASSERT_SPEED(ns_per_octonion, 20.0 * ci_scale, "[N=8 ] SoA SIMD Octonion batch (per-octonion)");
    }
}

int main(int argc, char** argv) {
    // 1. Determine RNG seed (random by default, or reproducible via NION_TEST_SEED)
    std::uint64_t seed = static_cast<std::uint64_t>(std::random_device{}()) ^
                         static_cast<std::uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if (const char* env_seed = std::getenv("NION_TEST_SEED")) {
        seed = std::strtoull(env_seed, nullptr, 10);
    }

    // 2. Determine number of random trials per function (default 500, configurable via NION_TEST_TRIALS)
    int num_trials = 500;
    if (const char* env_trials = std::getenv("NION_TEST_TRIALS")) {
        const int parsed = std::atoi(env_trials);
        if (parsed > 0) num_trials = parsed;
    }

    // 3. Optional --suite=<n2|n4|n8|n16|n32|simd|speed|all> filter for granular CTest subtests
    std::string suite = "all";
    for (int i = 1; i < argc; ++i) {
        if (std::strncmp(argv[i], "--suite=", 8) == 0) {
            suite = argv[i] + 8;
        }
    }

    std::mt19937_64 rng(seed);

    std::printf("=================================================================\n");
    std::printf("  Nion Automated CI/CD Accuracy & Benchmark Validation Suite\n");
    std::printf("  RNG Seed: %llu (replay with NION_TEST_SEED=%llu) | Trials: %d\n",
                static_cast<unsigned long long>(seed),
                static_cast<unsigned long long>(seed),
                num_trials);
    std::printf("=================================================================\n");

    if (suite == "all" || suite == "n2")    test_degree_2_complex(rng, num_trials);
    if (suite == "all" || suite == "n4")    test_degree_4_quaternion(rng, num_trials);
    if (suite == "all" || suite == "n8")    test_degree_8_octonion(rng, num_trials);
    if (suite == "all" || suite == "n16")   test_degree_16_sedenion(rng, num_trials);
    if (suite == "all" || suite == "n32")   test_degree_32_trigintaduonion(rng, num_trials);
    if (suite == "all" || suite == "simd")  test_soa_simd_batches(rng, num_trials);
    if (suite == "all" || suite == "speed") test_benchmark_performance(rng);

    std::printf("\n=================================================================\n");
    std::printf("  Summary: %d test groups executed, %d failed.\n", g_tests_run, g_tests_failed);
    std::printf("=================================================================\n");

    return (g_tests_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
