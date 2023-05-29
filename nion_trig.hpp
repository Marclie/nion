/*
 Copyright 2022 Marcus Dante Liebenthal

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#ifndef NION_TRIG_H
#define NION_TRIG_H

#include <cmath>
#include <limits>
#include <type_traits>
#include <iostream>
#include <stdexcept>
#include "nion.hpp"

/***************************
    *  NION HYPERBOLIC TRIGONOMETRIC FUNCTIONS *
    ***************************/

/**
 * @brief compute the hyperbolic sine of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic sine of.
 * @return The hyperbolic sine of the nion.
 * @details The hyperbolic sine of a nion is defined as sinh(z) = (e^z - e^-z) / 2.
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> sinh(const nion<T,N> &z);

/**
 * @brief compute the hyperbolic cosine of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic cosine of.
 * @return The hyperbolic cosine of the nion.
 * @details The hyperbolic cosine of a nion is defined as cosh(z) = (e^z + e^-z) / 2.
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> cosh(const nion<T,N> &z);

/**
 * @brief compute the hyperbolic tangent of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic tangent of.
 * @return The hyperbolic tangent of the nion.
 * @details The hyperbolic tangent of a nion is defined as tanh(z) = sinh(z) / cosh(z).
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> tanh(const nion<T,N> &z);

/**
 * @brief compute the hyperbolic cotangent of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic cotangent of.
 * @return The hyperbolic cotangent of the nion.
 * @details The hyperbolic cotangent of a nion is defined as coth(z) = 1 / tanh(z).
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> coth(const nion<T,N> &z);

/**
 * @brief compute the hyperbolic secant of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic secant of.
 * @return The hyperbolic secant of the nion.
 * @details The hyperbolic secant of a nion is defined as sech(z) = 1 / cosh(z).
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> sech(const nion<T,N> &z);

/**
 * @brief compute the hyperbolic cosecant of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic cosecant of.
 * @return The hyperbolic cosecant of the nion.
 * @details The hyperbolic cosecant of a nion is defined as csch(z) = 1 / sinh(z).
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> csch(const nion<T,N> &z);

/***************************
    *  NION TRIGONOMETRIC FUNCTIONS *
    ***************************/

/**
 * @brief compute the sine of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the sine of.
 * @return The sine of the nion.
 * @details The sine of the nion is defined as sin(z) = sin(r + v) = sin(r) * cosh(|v|) + cos(r) * sinh(|v|) * v/|v|.
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://en.wikipedia.org/wiki/Sine_and_cosine#Relationship_to_complex_numbers
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> sin(const nion<T,N> &z);

/**
 * @brief compute the cosine of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the cosine of.
 * @return The cosine of the nion.
 * @details The cosine of the nion is defined as cos(z) = cos(r + v) = cos(r) * cosh(|v|) - sin(r) * sinh(|v|) * v/|v|.
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://en.wikipedia.org/wiki/Sine_and_cosine#Relationship_to_complex_numbers
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> cos(const nion<T,N> &z);

/**
 * @brief compute the tangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the tangent of.
 * @return The tangent of the nion.
 * @details The tangent of the nion is defined as
 * tan(z) = tan(a + bi) = (tan(a) + tanh(b)i) / (1 - tan(a) * tanh(b) i).
 * @see https://en.wikipedia.org/wiki/Proofs_of_trigonometric_identities#Angle_sum_identities
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> tan(const nion<T,N> &z);

/**
 * @brief compute the cotangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the cotangent of.
 * @return The cotangent of the nion.
 * @details The cotangent of the nion is defined as cot(z) = 1 / tan(z).
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> cot(const nion<T,N> &z);

/**
 * @brief compute the secant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the secant of.
 * @return The secant of the nion.
 * @details The secant of the nion is defined as sec(z) = 1 / cos(z).
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> sec(const nion<T,N> &z);

/**
 * @brief compute the cosecant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the cosecant of.
 * @return The cosecant of the nion.
 * @details The cosecant of the nion is defined as csc(z) = 1 / sin(z).
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> csc(const nion<T,N> &z);

/***************************
    *  NION INVERSE HYPERBOLIC TRIGONOMETRIC FUNCTIONS *
    ***************************/

/**
 * @brief compute the inverse hyperbolic sine of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse hyperbolic sine of.
 * @return The inverse hyperbolic sine of the nion.
 * @details The inverse hyperbolic sine of the nion is defined as asinh(z) = ln(z + sqrt(z^2 + 1)).
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions#Principal_values_in_the_complex_plane
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> asinh(const nion<T,N> &z);

/**
 * @brief compute the inverse hyperbolic cosine of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse hyperbolic cosine of.
 * @return The inverse hyperbolic cosine of the nion.
 * @details The inverse hyperbolic cosine of the nion is defined as acosh(z) = sqrt(z-1) / sqrt(1-z) * acos(z).
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://mathworld.wolfram.com/InverseHyperbolicCosine.html
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> acosh(const nion<T,N> &z);

/**
 * @brief compute the inverse hyperbolic tangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse hyperbolic tangent of.
 * @return The inverse hyperbolic tangent of the nion.
 * @details The inverse hyperbolic tangent of the nion is defined as atanh(z) = atanh(r + v) = 1/2 * (ln(1 + z) - ln(1 - z))).
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://mathworld.wolfram.com/InverseHyperbolicTangent.html
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> atanh(const nion<T,N> &z);

/**
 * @brief compute the inverse hyperbolic cotangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse hyperbolic cotangent of.
 * @return The inverse hyperbolic cotangent of the nion.
 * @details The inverse hyperbolic cotangent of the nion is defined as acoth(z) = acoth(r + v) = 1/2 * (ln(1 + 1/z) - ln(1 - 1/z)).
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://mathworld.wolfram.com/InverseHyperbolicCotangent.html
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> acoth(const nion<T,N> &z);

/**
 * @brief compute the inverse hyperbolic secant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse hyperbolic secant of.
 * @return The inverse hyperbolic secant of the nion.
 * @details The inverse hyperbolic secant of the nion is defined as asech(z) = asech(r + v) = ln(sqrt(1/z - 1)*sqrt(1/z + 1) + 1/z).
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://mathworld.wolfram.com/InverseHyperbolicSecant.html
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> asech(const nion<T,N> &z);

/**
 * @brief compute the inverse hyperbolic cosecant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse hyperbolic cosecant of.
 * @return The inverse hyperbolic cosecant of the nion.
 * @details The inverse hyperbolic cosecant of the nion is defined as acsch(z) = acsch(r + v) = ln(sqrt(1+1/z^2) + 1/z).
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://mathworld.wolfram.com/InverseHyperbolicCosecant.html
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> acsch(const nion<T,N> &z);

/***************************
    *  NION INVERSE TRIGONOMETRIC FUNCTIONS *
    ***************************/

/**
 * @brief compute the inverse sine of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse sine of.
 * @return The inverse sine of the nion.
 * @details The inverse sine of the nion is defined as asin(z) = asin(r + v) = v/|v| * ln(sqrt(1 - z^2) - v/|v| * z).
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Logarithmic_forms
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> asin(const nion<T,N> &z);

/**
 * @brief compute the inverse cosine of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse cosine of.
 * @return The inverse cosine of the nion.
 * @details The inverse cosine of the nion is defined as acos(z) = pi/2 - asin(z).
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Extension_to_complex_plane
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> acos(const nion<T,N> &z);

/**
 * @brief compute the inverse tangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse tangent of.
 * @return The inverse tangent of the nion.
 * @details The inverse tangent of the nion is defined as atan(z) = v * atanh(-v*z) = -v/2 * log((v-z)/(v+z))
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Logarithmic_forms
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> atan(const nion<T,N> &z);

/**
 * @brief compute the inverse cotangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse cotangent of.
 * @return The inverse cotangent of the nion.
 * @details The inverse cotangent of the nion is defined as acot(z) = pi/2 - atan(z).
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Extension_to_complex_plane
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> acot(const nion<T,N> &z);

/**
 * @brief compute the inverse secant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse secant of.
 * @return The inverse secant of the nion.
 * @details The inverse secant of the nion is defined as asec(z) = acos(1/z).
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Extension_to_complex_plane
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> asec(const nion<T,N> &z);

/**
 * @brief compute the inverse cosecant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse cosecant of.
 * @return The inverse cosecant of the nion.
 * @details The inverse cosecant of the nion is defined as acsc(z) = asin(1/z).
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Extension_to_complex_plane
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> acsc(const nion<T,N> &z);

/**
 * @brief return atan2 between real and imaginary axis.
 * @tparam T type of the nion.
 * @param z The nion to compute the atan2 of.
 * @return The atan2 of the nion.
 */
template<arith_ops T, std::size_t N = 128, std::size_t M = N>
extern constexpr inline nion<T,N> atan2(const nion<T,N> &y, const nion<T,M> &x);

/***************************
    *  NION GAMMA FUNCTION *
    ***************************/

/**
 * @brief compute the gamma function of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the gamma function of.
 * @return The gamma function of the nion.
 * @details This uses an approximate formula for the gamma function of the nion:
 *     gamma(z) ≈ sqrt(2 π) e^(-z) sqrt(1/(z)) (1/(12 (z) - 1/(10 (z))) + z)^(z)
 * @see https://www.wolframalpha.com/input?i=gamma%28a+%2B+b+i%29
 */
template<arith_ops T, std::size_t N = 128>
extern constexpr inline nion<T,N> gamma(const nion<T,N> &z);

#include "nion_trig.cpp"
#endif //NION_TRIG_H

