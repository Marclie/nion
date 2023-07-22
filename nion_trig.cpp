/*
 Copyright 2023 Marcus Dante Liebenthal

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

#ifndef NION_TRIG_CPP
#define NION_TRIG_CPP

#include "nion_trig.hpp"


/*******************************************
*  NION HYPERBOLIC TRIGONOMETRIC FUNCTIONS *
********************************************/

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> sinh(const nion<T,N> &z) {
    // get polar form of nion
    nion<T,N> i = z.imag();

    // make unit vector
    T i_abs = i.abs();
    T i_norm = sqrt(i_abs);

    // calculate scalars
    T e_z = exp(z.elem_[0]) / 2.0l;
    T e_mz = exp(-z.elem_[0]) / 2.0l;

    // compute denorm_min
    T denorm_min = T(0); // default value
    if constexpr (std::is_arithmetic_v<T>) {
        denorm_min = std::numeric_limits<T>::denorm_min(); // if T is arithmetic, use its denorm_min
    }

    // compute exponential of nion
    if (i_abs <= denorm_min) {
        nion<T,N> sin_nion = i * ((e_z + e_mz) * sin(i_norm));
        sin_nion += cos(i_norm) * (e_z - e_mz);
        return sin_nion;
    } else {
        nion<T,N> sin_nion = i * ((e_z + e_mz) * sin(i_norm) / i_norm);
        sin_nion += cos(i_norm) * (e_z - e_mz);
        return sin_nion;
    }

}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> cosh(const nion<T,N> &z) {
    // get polar form of nion
    nion<T,N> i = z.imag();

    // make unit vector
    T i_abs = i.abs();
    T i_norm = sqrt(i_abs);

    // calculate scalars
    T e_z = exp(z.elem_[0]) / 2.0l;
    T e_mz = exp(-z.elem_[0]) / 2.0l;

    // compute denorm_min
    T denorm_min = T(0); // default value
    if constexpr (std::is_arithmetic_v<T>) {
        denorm_min = std::numeric_limits<T>::denorm_min(); // if T is arithmetic, use its denorm_min
    }

    // compute exponential of nion
    if (i_abs <= denorm_min) {
        nion<T,N> cos_nion = i * ((e_z - e_mz) * sin(i_norm));
        cos_nion += cos(i_norm) * (e_z + e_mz);
        return cos_nion;
    } else {
        nion<T,N> cos_nion = i * ((e_z - e_mz) * sin(i_norm) / i_norm);
        cos_nion += cos(i_norm) * (e_z + e_mz);
        return cos_nion;
    }
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> tanh(const nion<T,N> &z) {
    return (exp(z*2.0l) - 1.0l) / (exp(z*2.0l) + 1.0l);
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> coth(const nion<T,N> &z) {
    return tanh(z).inv();
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> sech(const nion<T,N> &z) {
    return cosh(z).inv();
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> csch(const nion<T,N> &z) {
    return sinh(z).inv();
}

/********************************
*  NION TRIGONOMETRIC FUNCTIONS *
*********************************/


template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> sin(const nion<T,N> &z) {
    // get the polar form of the nion
    T r = real(z);
    nion<T,N> i = imag(z);

    // make unit vector
    T i_abs = i.abs();
    T i_norm = sqrt(i_abs);

    // compute denorm_min
    T denorm_min = T(0); // default value
    if constexpr (std::is_arithmetic_v<T>) {
        denorm_min = std::numeric_limits<T>::denorm_min(); // if T is arithmetic, use its denorm_min
    }

    // compute the sine of the nion
    if (i_abs <= denorm_min)
        return i * (sinh(i_norm) * cos(r)) + sin(r) * cosh(i_norm);
    else
        return i * (sinh(i_norm) * cos(r) / i_norm) + sin(r) * cosh(i_norm);
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> cos(const nion<T,N> &z) {
    // get the polar form of the nion
    T r = real(z);
    nion<T,N> i = imag(z);

    // make unit vector
    T i_abs = i.abs();
    T i_norm = sqrt(i_abs);

    // compute denorm_min
    T denorm_min = T(0); // default value
    if constexpr (std::is_arithmetic_v<T>) {
        denorm_min = std::numeric_limits<T>::denorm_min(); // if T is arithmetic, use its denorm_min
    }

    // compute the cosine of the nion
    if (i_abs <= denorm_min)
        return -i * (sinh(i_norm) * sin(r)) + cos(r) * cosh(i_norm);
    else
        return -i * (sinh(i_norm) * sin(r) / i_norm) + cos(r) * cosh(i_norm);
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> tan(const nion<T,N> &z) {
    // get the polar form of the nion
    T r = real(z);
    nion<T,N> i = imag(z);

    // make unit vector
    T i_abs = i.abs();
    T i_norm = sqrt(i_abs);

    // compute denorm_min
    T denorm_min = T(0); // default value
    if constexpr (std::is_arithmetic_v<T>) {
        denorm_min = std::numeric_limits<T>::denorm_min(); // if T is arithmetic, use its denorm_min
    }

    // compute the tangent of the nion
    if (i_norm <= denorm_min)
        return (tan(r) + i * tanh(i_norm)) / (1 - i * (tan(r) * tanh(i_norm)));
    else
        return (tan(r) + i * (tanh(i_norm) / i_norm)) / (1 - i * (tan(r) / i_norm * tanh(i_norm)));
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> cot(const nion<T,N> &z) {
    return tan(z).inv();
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> sec(const nion<T,N> &z) {
    return cos(z).inv();
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> csc(const nion<T,N> &z) {
    return sin(z).inv();
}

/***************************************************
*  NION INVERSE HYPERBOLIC TRIGONOMETRIC FUNCTIONS *
****************************************************/

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> asinh(const nion<T,N> &z) {
    return log(z + sqrt(sqr(z) + 1.0l));
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> acosh(const nion<T,N> &z) {
    // compute the inverse hyperbolic cosine of the nion
    return 2 * log(sqrt((z-1) * 0.5l) + sqrt((z+1.0l) * 0.5l));

}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> atanh(const nion<T,N> &z) {
    return (log(1.0l + z) - log(1.0l - z)) * 0.5l;
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> acoth(const nion<T,N> &z) {
    return (log(1.0l + inv(z)) - log(1.0l - inv(z))) * 0.5l;
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> asech(const nion<T,N> &z) {
    return log(sqrt(1.0l/sqr(z) - 1.0l) + inv(z));
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> acsch(const nion<T,N> &z) {
    return log(sqrt(1.0l + 1.0l/sqr(z)) + inv(z));
}

/****************************************
*  NION INVERSE TRIGONOMETRIC FUNCTIONS *
*****************************************/

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> atri(const nion<T,N> &z, const nion<T,N> &a, const nion<T,N> &b, const nion<T,N> &c){
    // get the polar form of the nion
    auto [z_abs, i_abs, i] = z.polar();

    return -i * log((a + b * i)/c);
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> asin(const nion<T,N> &z) {


    // get the polar form of the nion
    T r = real(z);
    auto [mag2, phase2, i] = z.polar(); // decompose the nion into its polar form
    T phase = sqrt(phase2); // compute the norm of the imaginary part

//    return i * asinh(-i*r + phase);

    // compute the inv sine of the nion
    return -i * log(sqrt(1.0l - sqr(z)) + (i * r) - phase);
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> acos(const nion<T,N> &z) {
    return M_PI_2 - asin(z);
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> atan(const nion<T,N> &z) {
    // get the polar form of the nion
    T r = real(z);
    nion<T,N> i = z.imag();

    // make unit vector
    T i_abs = i.abs();
    T i_norm = sqrt(i_abs);

    // compute denorm_min
    T denorm_min = T(0); // default value
    if constexpr (std::is_arithmetic_v<T>) {
        denorm_min = std::numeric_limits<T>::denorm_min(); // if T is arithmetic, use its denorm_min
    }

    if (i_abs > denorm_min)
        i /= i_norm;

    // compute the inv tangent of the nion:
    T one = 1;
    return 0.5l * i * (-log((one - i_norm) + (i * r) ) + log((one + i_norm) + (i * -r) ));
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> acot(const nion<T,N> &z) {
    return M_PI_2 - atan(z);
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> asec(const nion<T,N> &z) {
    return acos(inv(z));
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> acsc(const nion<T,N> &z) {
    return asin(inv(z));
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> atan2(const nion<T,N> &y, const nion<T,N> &x) {
    return atan(y / x);
}

/***************************
 *   NION GAMMA FUNCTION   *
 **************************/

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> gamma(const nion<T,N> &z) {
    // compute the gamma function of the nion
    return exp(-z) * sqrt(inv(z))
           * pow(1.0l / (12.0l * z - inv(10.0l * z)) + z, z)
           * sqrt(2.0l * M_PI);
}

#endif // NION_TRIG_CPP