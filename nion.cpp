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

#ifndef NION_CPP
#define NION_CPP

#include <cmath>
#include <limits>
#include <type_traits>
#include <iostream>
#include <stdexcept>
#include "nion.hpp"
#include <tuple>

#define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << message << ": file=" << __FILE__ \
                      << ", line=" << __LINE__ << std::endl; \
            std::terminate(); \
        } \
    } while (false)



template<arith_ops T, std::size_t N>
template< arith_ops S, arith_ops U, std::size_t M, std::size_t P> requires (std::is_convertible_v<S,T> && std::is_convertible_v<U,T>)
constexpr inline nion<T,N> nion<T,N>::make_pair(const nion<S,M> &a, const nion<U,P> &b) {

    D a_size = a.size(), b_size = b.size();
    ASSERT(a_size > 0 && b_size > 0, "The sizes of the nion pair (a, b) must both be greater than zero.");

    // initialize the nion pair
    nion<T,N> pair(a_size + b_size);

    // copy the values into the nion
    pair.elem_.copy(a.elem_, a_size);
    pair.elem_.copy(b.elem_, b_size, a_size);
    return pair;
}

template<arith_ops T, std::size_t N>
template<integral_type Z>
constexpr inline void nion<T,N>::resize(Z size) {
    ASSERT(size > 0, "new nion size must be greater than zero");
    D old_size = elem_.size_;
    elem_.expand(size);

    // zero out the new elements
    elem_.fill(T(), old_size);

    size_ = this->size();

}

/************************************
*  ASSIGNMENT AND ADDITION OPERATORS
*************************************/

template<arith_ops T, std::size_t N>
template<std::size_t M>
constexpr inline void nion<T,N>::operator+=(const nion<T,M> &other) {

    // find the common integer type
    using E = typename nion<T,M>::D;
    using uint_com = std::common_type_t<D, E>;

    // find the min max sizes
    uint_com small_size = std::min(elem_.size_, other.elem_.size_);
    uint_com big_size = std::max(elem_.size_,  other.elem_.size_);
    resize(big_size);

    // for the first smaller_size elements, add the other nion's elements to this nion's elements.
    for (uint_com i = 0; i < small_size; i++)
        elem_[i] += other.elem_[i];

    // copy the remaining values (if any) from the other nion to this nion.
    if (small_size != big_size) {
        if (other.elem_.size_ == big_size) {
            elem_.copy(other.elem_, big_size - small_size, small_size, small_size);
        }
    }

}

template<arith_ops T, std::size_t N>
template<std::size_t M>
constexpr inline void nion<T,N>::operator-=(const nion<T,M> &other) {
    // find the common integer type
    using E = typename nion<T,M>::D;
    using uint_com = std::common_type_t<D, E>;

    // find the min max sizes
    uint_com small_size = std::min(elem_.size_, other.elem_.size_);
    uint_com big_size = std::max(elem_.size_,  other.elem_.size_);
    resize(big_size);

    // for the first smaller_size elements, add the other nion's elements to this nion's elements.
    for (uint_com i = 0; i < small_size; i++)
        elem_[i] -= other.elem_[i];

    // copy and negate the remaining values (if any) from the other nion to this nion.
    if (small_size != big_size) {
        if (other.elem_.size_ == big_size) {
            for (uint_com i = small_size; i < big_size; i++)
                elem_[i] = -other.elem_[i];
        }
    }
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline void nion<T,N>::operator+=(S scalar) {
    elem_[0] += scalar;
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline void nion<T,N>::operator-=(S scalar) {
    elem_[0] -= scalar;
}

/************************************
*        ADDITION OPERATORS
*************************************/

// TODO: make this work for different types than T

template<arith_ops T, std::size_t N>
template<std::size_t M>
constexpr inline nion<T,N> nion<T,N>::operator+(const nion<T,M> &other) const {
    nion<T,N> sum(*this);
    sum += other;
    return sum;
}

template<arith_ops T, std::size_t N>
template<std::size_t M>
constexpr inline nion<T,N> nion<T,N>::operator-(const nion<T,M> &other) const {
    nion<T,N> diff(*this);
    diff -= other;
    return diff;
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> nion<T,N>::operator+(S scalar) const {

    // create a nion to store the sum
    nion<T,N> sum(*this);

    // add the scalar to the real component of the sum nion.
    sum.elem_[0] += scalar;
    return sum;
}

template<arith_ops T, std::size_t N, not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> operator+(S scalar, const nion<T,N> &z) {
    return z + scalar;
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> &nion<T,N>::operator++() {
    elem_[0]++;
    return *this;
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> nion<T,N>::operator-(S scalar) const {
    return *this + -scalar;
}

template<arith_ops T, std::size_t N, not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> operator-(S scalar, const nion<T,N> &z) {
    return scalar + -z;
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> &nion<T,N>::operator--() {
    elem_[0]--;
    return *this;
}

/******************************************
*  ASSIGNMENT AND MULTIPLICATION OPERATORS
*******************************************/

template<arith_ops T, std::size_t N>
template<std::size_t M>
constexpr inline void nion<T,N>::operator*=(const nion<T,M> &other) {
    *this = *this * other;
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline void nion<T,N>::operator*=(S scalar) {
    D size = this->size();
    for (D i = 0; i < size; i++)
        elem_[i] *= scalar;
}

template<arith_ops T, std::size_t N>
template<std::size_t M>
constexpr inline void nion<T,N>::operator/=(const nion<T,M> &other) {
    *this = *this / other;
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline void nion<T,N>::operator/=(S scalar) {
    D size = this->size();
    for (D i = 0; i < size; i++)
        elem_[i] /= scalar;
}

/******************************************
*        MULTIPLICATION OPERATORS
*******************************************/

// create concept to check if T has a `.conj()` method
template<typename T>
concept has_conj = requires(T a) {
    a.conj();
};

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> nion<T,N>::conj() const {
    D size = this->size();
    nion<T,N> conjugate(size); // create a nion to store the conjugate

    // conjugate the first element if T has a `.conj()` method
    if constexpr (has_conj<T>) conjugate.elem_[0] = elem_[0].conj();
    else conjugate.elem_[0] = elem_[0]; // else copy the first element

    // negate all components except the first
    for (D i = 1; i < size; i++)
        conjugate.elem_[i] = -elem_[i]; // negate the component

    return conjugate;
}

template<arith_ops T, std::size_t N>
constexpr inline void nion<T,N>::conj_inplace() {
    // conjugate the first element if T has a `.conj()` method
    if constexpr (has_conj<T>)
        elem_[0] = elem_[0].conj();

    // negate all components except the first
    D size = this->size();
    for (D i = 1; i < size; i++)
        elem_[i] = -elem_[i]; // negate the component
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> nion<T,N>::operator-() const {
    D size = this->size();
    nion<T,N> negated(size); // copy this nion

    // negate all components
    for (D i = 0; i < size; i++)
        negated.elem_[i] = -elem_[i]; // negate the component

    return negated;
}

template<arith_ops T, std::size_t N>
template<arith_ops S, std::size_t M> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> nion<T,N>::operator*(const nion<S,M> &other) const {

    // find common integral type
    using E = typename nion<S,M>::D;
    using uint_com = std::common_type_t<D, E>;

    uint_com size = this->size(), other_size = other.size();
    if (size == other_size) {

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds" // Disable the array out-of-bounds warning
#elif defined(__INTEL_COMPILER)
        #pragma warning(push)
            #pragma warning(disable: 2259) // Disable the array out of bounds warning
#endif
        const T* Lhs = elem_.vals_;
        const S* Rhs = other.elem_.vals_;
        switch (size) {
            case 1: // if this size is 1, then the product is just the scalar product.

                return {Lhs[0] * Rhs[0]};

// if FULL_RECURSION is defined, then the compiler will use recursion to calculate the product of nions for any size.
#ifndef FULL_RECURSION
            case 2: // hard-coded complex product
                return {
                    {Lhs[0] * Rhs[0] - Lhs[1] * Rhs[1]},
                    {Lhs[1] * Rhs[0] + Lhs[0] * Rhs[1]}
                };


            case 4: // hard-coded quaternion product
                return {
                    {Lhs[0] * Rhs[0] - Lhs[1] * Rhs[1] - Lhs[2] * Rhs[2] - Lhs[3] * Rhs[3]},
                    {Lhs[1] * Rhs[0] + Lhs[0] * Rhs[1] - Lhs[3] * Rhs[2] + Lhs[2] * Rhs[3]},
                    {Lhs[2] * Rhs[0] + Lhs[3] * Rhs[1] + Lhs[0] * Rhs[2] - Lhs[1] * Rhs[3]},
                    {Lhs[3] * Rhs[0] - Lhs[2] * Rhs[1] + Lhs[1] * Rhs[2] + Lhs[0] * Rhs[3]}
                };

            case 8: // hard-coded octonion product ( I know, it's a bit much )
                return {
                    {Lhs[0] * Rhs[0] - Lhs[1] * Rhs[1] - Lhs[2] * Rhs[2] - Lhs[3] * Rhs[3] - Lhs[4] * Rhs[4] - Lhs[5] * Rhs[5] - Lhs[6] * Rhs[6] - Lhs[7] * Rhs[7]},
                    {Lhs[1] * Rhs[0] + Lhs[0] * Rhs[1] - Lhs[3] * Rhs[2] + Lhs[2] * Rhs[3] - Lhs[5] * Rhs[4] + Lhs[4] * Rhs[5] + Lhs[7] * Rhs[6] - Lhs[6] * Rhs[7]},
                    {Lhs[2] * Rhs[0] + Lhs[3] * Rhs[1] + Lhs[0] * Rhs[2] - Lhs[1] * Rhs[3] - Lhs[6] * Rhs[4] - Lhs[7] * Rhs[5] + Lhs[4] * Rhs[6] + Lhs[5] * Rhs[7]},
                    {Lhs[3] * Rhs[0] - Lhs[2] * Rhs[1] + Lhs[1] * Rhs[2] + Lhs[0] * Rhs[3] - Lhs[7] * Rhs[4] + Lhs[6] * Rhs[5] - Lhs[5] * Rhs[6] + Lhs[4] * Rhs[7]},
                    {Lhs[4] * Rhs[0] + Lhs[5] * Rhs[1] + Lhs[6] * Rhs[2] + Lhs[7] * Rhs[3] + Lhs[0] * Rhs[4] - Lhs[1] * Rhs[5] - Lhs[2] * Rhs[6] - Lhs[3] * Rhs[7]},
                    {Lhs[5] * Rhs[0] - Lhs[4] * Rhs[1] + Lhs[7] * Rhs[2] - Lhs[6] * Rhs[3] + Lhs[1] * Rhs[4] + Lhs[0] * Rhs[5] + Lhs[3] * Rhs[6] - Lhs[2] * Rhs[7]},
                    {Lhs[6] * Rhs[0] - Lhs[7] * Rhs[1] - Lhs[4] * Rhs[2] + Lhs[5] * Rhs[3] + Lhs[2] * Rhs[4] - Lhs[3] * Rhs[5] + Lhs[0] * Rhs[6] + Lhs[1] * Rhs[7]},
                    {Lhs[7] * Rhs[0] + Lhs[6] * Rhs[1] - Lhs[5] * Rhs[2] - Lhs[4] * Rhs[3] + Lhs[3] * Rhs[4] + Lhs[2] * Rhs[5] - Lhs[1] * Rhs[6] + Lhs[0] * Rhs[7]}
                };

            case 16: // hard-coded sedenion product ( I strongly recommend using nowrap if you're looking at this... )
                return {
                    {Lhs[ 0] * Rhs[ 0] - Lhs[ 1] * Rhs[ 1] - Lhs[ 2] * Rhs[ 2] - Lhs[ 3] * Rhs[ 3] - Lhs[ 4] * Rhs[ 4] - Lhs[ 5] * Rhs[ 5] - Lhs[ 6] * Rhs[ 6] - Lhs[ 7] * Rhs[ 7] - Lhs[ 8] * Rhs[ 8] - Lhs[ 9] * Rhs[ 9] - Lhs[10] * Rhs[10] - Lhs[11] * Rhs[11] - Lhs[12] * Rhs[12] - Lhs[13] * Rhs[13] - Lhs[14] * Rhs[14] - Lhs[15] * Rhs[15]},
                    {Lhs[ 0] * Rhs[ 1] + Lhs[ 1] * Rhs[ 0] + Lhs[ 2] * Rhs[ 3] - Lhs[ 3] * Rhs[ 2] + Lhs[ 4] * Rhs[ 5] - Lhs[ 5] * Rhs[ 4] - Lhs[ 6] * Rhs[ 7] + Lhs[ 7] * Rhs[ 6] + Lhs[ 8] * Rhs[ 9] - Lhs[ 9] * Rhs[ 8] - Lhs[10] * Rhs[11] + Lhs[11] * Rhs[10] - Lhs[12] * Rhs[13] + Lhs[13] * Rhs[12] + Lhs[14] * Rhs[15] - Lhs[15] * Rhs[14]},
                    {Lhs[ 0] * Rhs[ 2] - Lhs[ 1] * Rhs[ 3] + Lhs[ 2] * Rhs[ 0] + Lhs[ 3] * Rhs[ 1] + Lhs[ 4] * Rhs[ 6] + Lhs[ 5] * Rhs[ 7] - Lhs[ 6] * Rhs[ 4] - Lhs[ 7] * Rhs[ 5] + Lhs[ 8] * Rhs[10] + Lhs[ 9] * Rhs[11] - Lhs[10] * Rhs[ 8] - Lhs[11] * Rhs[ 9] - Lhs[12] * Rhs[14] - Lhs[13] * Rhs[15] + Lhs[14] * Rhs[12] + Lhs[15] * Rhs[13]},
                    {Lhs[ 0] * Rhs[ 3] + Lhs[ 1] * Rhs[ 2] - Lhs[ 2] * Rhs[ 1] + Lhs[ 3] * Rhs[ 0] + Lhs[ 4] * Rhs[ 7] - Lhs[ 5] * Rhs[ 6] + Lhs[ 6] * Rhs[ 5] - Lhs[ 7] * Rhs[ 4] + Lhs[ 8] * Rhs[11] - Lhs[ 9] * Rhs[10] + Lhs[10] * Rhs[ 9] - Lhs[11] * Rhs[ 8] - Lhs[12] * Rhs[15] + Lhs[13] * Rhs[14] - Lhs[14] * Rhs[13] + Lhs[15] * Rhs[12]},
                    {Lhs[ 0] * Rhs[ 4] - Lhs[ 1] * Rhs[ 5] - Lhs[ 2] * Rhs[ 6] - Lhs[ 3] * Rhs[ 7] + Lhs[ 4] * Rhs[ 0] + Lhs[ 5] * Rhs[ 1] + Lhs[ 6] * Rhs[ 2] + Lhs[ 7] * Rhs[ 3] + Lhs[ 8] * Rhs[12] + Lhs[ 9] * Rhs[13] + Lhs[10] * Rhs[14] + Lhs[11] * Rhs[15] - Lhs[12] * Rhs[ 8] - Lhs[13] * Rhs[ 9] - Lhs[14] * Rhs[10] - Lhs[15] * Rhs[11]},
                    {Lhs[ 0] * Rhs[ 5] + Lhs[ 1] * Rhs[ 4] - Lhs[ 2] * Rhs[ 7] + Lhs[ 3] * Rhs[ 6] - Lhs[ 4] * Rhs[ 1] + Lhs[ 5] * Rhs[ 0] - Lhs[ 6] * Rhs[ 3] + Lhs[ 7] * Rhs[ 2] + Lhs[ 8] * Rhs[13] - Lhs[ 9] * Rhs[12] + Lhs[10] * Rhs[15] - Lhs[11] * Rhs[14] + Lhs[12] * Rhs[ 9] - Lhs[13] * Rhs[ 8] + Lhs[14] * Rhs[11] - Lhs[15] * Rhs[10]},
                    {Lhs[ 0] * Rhs[ 6] + Lhs[ 1] * Rhs[ 7] + Lhs[ 2] * Rhs[ 4] - Lhs[ 3] * Rhs[ 5] - Lhs[ 4] * Rhs[ 2] + Lhs[ 5] * Rhs[ 3] + Lhs[ 6] * Rhs[ 0] - Lhs[ 7] * Rhs[ 1] + Lhs[ 8] * Rhs[14] - Lhs[ 9] * Rhs[15] - Lhs[10] * Rhs[12] + Lhs[11] * Rhs[13] + Lhs[12] * Rhs[10] - Lhs[13] * Rhs[11] - Lhs[14] * Rhs[ 8] + Lhs[15] * Rhs[ 9]},
                    {Lhs[ 0] * Rhs[ 7] - Lhs[ 1] * Rhs[ 6] + Lhs[ 2] * Rhs[ 5] + Lhs[ 3] * Rhs[ 4] - Lhs[ 4] * Rhs[ 3] - Lhs[ 5] * Rhs[ 2] + Lhs[ 6] * Rhs[ 1] + Lhs[ 7] * Rhs[ 0] + Lhs[ 8] * Rhs[15] + Lhs[ 9] * Rhs[14] - Lhs[10] * Rhs[13] - Lhs[11] * Rhs[12] + Lhs[12] * Rhs[11] + Lhs[13] * Rhs[10] - Lhs[14] * Rhs[ 9] - Lhs[15] * Rhs[ 8]},
                    {Lhs[ 0] * Rhs[ 8] - Lhs[ 1] * Rhs[ 9] - Lhs[ 2] * Rhs[10] - Lhs[ 3] * Rhs[11] - Lhs[ 4] * Rhs[12] - Lhs[ 5] * Rhs[13] - Lhs[ 6] * Rhs[14] - Lhs[ 7] * Rhs[15] + Lhs[ 8] * Rhs[ 0] + Lhs[ 9] * Rhs[ 1] + Lhs[10] * Rhs[ 2] + Lhs[11] * Rhs[ 3] + Lhs[12] * Rhs[ 4] + Lhs[13] * Rhs[ 5] + Lhs[14] * Rhs[ 6] + Lhs[15] * Rhs[ 7]},
                    {Lhs[ 0] * Rhs[ 9] + Lhs[ 1] * Rhs[ 8] - Lhs[ 2] * Rhs[11] + Lhs[ 3] * Rhs[10] - Lhs[ 4] * Rhs[13] + Lhs[ 5] * Rhs[12] + Lhs[ 6] * Rhs[15] - Lhs[ 7] * Rhs[14] - Lhs[ 8] * Rhs[ 1] + Lhs[ 9] * Rhs[ 0] - Lhs[10] * Rhs[ 3] + Lhs[11] * Rhs[ 2] - Lhs[12] * Rhs[ 5] + Lhs[13] * Rhs[ 4] + Lhs[14] * Rhs[ 7] - Lhs[15] * Rhs[ 6]},
                    {Lhs[ 0] * Rhs[10] + Lhs[ 1] * Rhs[11] + Lhs[ 2] * Rhs[ 8] - Lhs[ 3] * Rhs[ 9] - Lhs[ 4] * Rhs[14] - Lhs[ 5] * Rhs[15] + Lhs[ 6] * Rhs[12] + Lhs[ 7] * Rhs[13] - Lhs[ 8] * Rhs[ 2] + Lhs[ 9] * Rhs[ 3] + Lhs[10] * Rhs[ 0] - Lhs[11] * Rhs[ 1] - Lhs[12] * Rhs[ 6] - Lhs[13] * Rhs[ 7] + Lhs[14] * Rhs[ 4] + Lhs[15] * Rhs[ 5]},
                    {Lhs[ 0] * Rhs[11] - Lhs[ 1] * Rhs[10] + Lhs[ 2] * Rhs[ 9] + Lhs[ 3] * Rhs[ 8] - Lhs[ 4] * Rhs[15] + Lhs[ 5] * Rhs[14] - Lhs[ 6] * Rhs[13] + Lhs[ 7] * Rhs[12] - Lhs[ 8] * Rhs[ 3] - Lhs[ 9] * Rhs[ 2] + Lhs[10] * Rhs[ 1] + Lhs[11] * Rhs[ 0] - Lhs[12] * Rhs[ 7] + Lhs[13] * Rhs[ 6] - Lhs[14] * Rhs[ 5] + Lhs[15] * Rhs[ 4]},
                    {Lhs[ 0] * Rhs[12] + Lhs[ 1] * Rhs[13] + Lhs[ 2] * Rhs[14] + Lhs[ 3] * Rhs[15] + Lhs[ 4] * Rhs[ 8] - Lhs[ 5] * Rhs[ 9] - Lhs[ 6] * Rhs[10] - Lhs[ 7] * Rhs[11] - Lhs[ 8] * Rhs[ 4] + Lhs[ 9] * Rhs[ 5] + Lhs[10] * Rhs[ 6] + Lhs[11] * Rhs[ 7] + Lhs[12] * Rhs[ 0] - Lhs[13] * Rhs[ 1] - Lhs[14] * Rhs[ 2] - Lhs[15] * Rhs[ 3]},
                    {Lhs[ 0] * Rhs[13] - Lhs[ 1] * Rhs[12] + Lhs[ 2] * Rhs[15] - Lhs[ 3] * Rhs[14] + Lhs[ 4] * Rhs[ 9] + Lhs[ 5] * Rhs[ 8] + Lhs[ 6] * Rhs[11] - Lhs[ 7] * Rhs[10] - Lhs[ 8] * Rhs[ 5] - Lhs[ 9] * Rhs[ 4] + Lhs[10] * Rhs[ 7] - Lhs[11] * Rhs[ 6] + Lhs[12] * Rhs[ 1] + Lhs[13] * Rhs[ 0] + Lhs[14] * Rhs[ 3] - Lhs[15] * Rhs[ 2]},
                    {Lhs[ 0] * Rhs[14] - Lhs[ 1] * Rhs[15] - Lhs[ 2] * Rhs[12] + Lhs[ 3] * Rhs[13] + Lhs[ 4] * Rhs[10] - Lhs[ 5] * Rhs[11] + Lhs[ 6] * Rhs[ 8] + Lhs[ 7] * Rhs[ 9] - Lhs[ 8] * Rhs[ 6] - Lhs[ 9] * Rhs[ 7] - Lhs[10] * Rhs[ 4] + Lhs[11] * Rhs[ 5] + Lhs[12] * Rhs[ 2] - Lhs[13] * Rhs[ 3] + Lhs[14] * Rhs[ 0] + Lhs[15] * Rhs[ 1]},
                    {Lhs[ 0] * Rhs[15] + Lhs[ 1] * Rhs[14] - Lhs[ 2] * Rhs[13] - Lhs[ 3] * Rhs[12] + Lhs[ 4] * Rhs[11] + Lhs[ 5] * Rhs[10] - Lhs[ 6] * Rhs[ 9] + Lhs[ 7] * Rhs[ 8] - Lhs[ 8] * Rhs[ 7] + Lhs[ 9] * Rhs[ 6] - Lhs[10] * Rhs[ 5] - Lhs[11] * Rhs[ 4] + Lhs[12] * Rhs[ 3] + Lhs[13] * Rhs[ 2] - Lhs[14] * Rhs[ 1] + Lhs[15] * Rhs[ 0]}
                };
#endif
        } // if size is not 1, 2, 8, or 16, we need to do the recursive product.

#ifdef __GNUC__
#pragma GCC diagnostic pop
#elif defined(__INTEL_COMPILER)
#pragma warning(pop)
#endif

    } else { // If the matrix is not square, we need to do the recursive product.

        // if any of the sizes are 1, then the product is just the scalar product.
        if (size == 1) return other * elem_[0];
        if (other_size == 1) return *this * other.elem_[0];

    }

    /// ********** Recursive Product ********** ///

    // get half-sizes for this and other
    uint_com shalf = size - (size >> 1);
    uint_com other_shalf = other_size - (other_size >> 1);

    nion<T,N> a(shalf),       b(size - shalf),
              c(other_shalf), d(other_size - other_shalf);

    // copy the elements of this and other into the new split nions
    a.elem_.copy(elem_, shalf);             b.elem_.copy(elem_, size - shalf, 0, shalf);
    c.elem_.copy(other.elem_, other_shalf); d.elem_.copy(other.elem_, other_size - other_shalf, 0, other_shalf);

    auto a_c = a * c; // compute the product of the top-left and bottom-left matrices
    auto d_a = d * a; // compute the product of the bottom-right and top-left matrices

    // conjugate the c and d matrices in-place (the fastest way to do it)
    c.conj_inplace(); d.conj_inplace();

    auto dH_b = d * b; // compute the product of the bottom-right and top-right matrices
    auto b_cH = b * c; // compute the product of the top-right and bottom-left matrices

    // compute and return the product
    return make_pair(
        a_c - dH_b,
        d_a + b_cH
    );
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> nion<T,N>::operator*(S scalar) const {
    nion<T,N> product(*this);

    // compute the product of each element of the nion with the scalar
    D size = this->size();
    for (D i = 0; i < size; i++)
        product.elem_[i] *= scalar;

    return product;
}

template<arith_ops T, std::size_t N, not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> operator*(S scalar, const nion<T,N> &z) {
    nion<T,N> product(z);

    // compute the product of each element of the nion with the scalar
    using D = typename nion<T,N>::D;
    D size = z.size();
    for (D i = 0; i < size; i++)
        product.elem_[i] = scalar * product.elem_[i];

    return product;
}

template<arith_ops T, std::size_t N>
constexpr inline T nion<T,N>::abs() const {

    T absVal = elem_[0] * elem_[0];

    D size = this->size();
    for (D i = 1; i < size; i++)
        absVal += elem_[i] * elem_[i];

    return absVal;
}

template<arith_ops T, std::size_t N>
constexpr inline T nion<T,N>::norm() const {
    return sqrt(this->abs());
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> nion<T,N>::inv() const {

    nion<T,N> inverse(*this);

    T absolute = this->abs();
    inverse.elem_[0] /= absolute;
    D size = this->size();
    for (D i = 1; i < size; i++)
        inverse.elem_[i] /= -absolute;

    return inverse;
}

template<arith_ops T, std::size_t N>
template<std::size_t M>
constexpr inline nion<T,N> nion<T,N>::operator/(const nion<T,M> &other) const {
    return *this * other.inv();
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> nion<T,N>::operator/(S scalar) const {
    nion<T,N> quotient(*this);

    // compute the product of each element of the nion with the scalar
    D size = this->size();
    for (D i = 0; i < size; i++)
        quotient.elem_[i] /= scalar;

    return quotient;
}

template<arith_ops T, std::size_t N, not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> operator/(S scalar, const nion<T,N> &z) {
    nion<T,N> quotient = z.inv();

    // compute the product of each element of the nion with the scalar
    using D = typename nion<T,N>::D;
    D size = z.size();
    for (D i = 0; i < size; i++)
        quotient.elem_[i] = scalar * quotient.elem_[i];

    return quotient;
}


/******************************************
*            ACCESSOR FUNCTIONS
*******************************************/

template<arith_ops T, std::size_t N>
constexpr inline T &nion<T,N>::operator[](D index) {
    return elem_[index];
}
template<arith_ops T, std::size_t N>
constexpr inline const T &nion<T,N>::operator[](D index) const {
    return elem_[index];
}

template<arith_ops T, std::size_t N>
constexpr inline T nion<T,N>::real() const { return elem_[0]; }

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> nion<T,N>::imag() const {
    nion<T,N> imag(*this);
    imag.elem_[0] = 0;
    return imag;
}

template<arith_ops T, std::size_t N>
constexpr inline T nion<T,N>::arg() const {
    // get the polar form of nion
    T r = elem_[0];

    // compute norms
    T z_abs = abs();
    T i_norm = sqrt(z_abs - r * r);

    // compute argument
    T theta = atan2(i_norm, r);
    return theta;
}

template<arith_ops T, std::size_t N>
constexpr inline std::tuple<T, T, nion<T,N>> nion<T,N>::polar() const {
    // get imaginary part of nion
    nion<T,N> unit = imag();

    // compute norm
    T mag2 = abs();
    T phase2 = mag2 - elem_.vals_[0] * elem_.vals_[0];

    // compute denorm_min
    T denorm_min = T(0); // default value
    if constexpr (std::is_arithmetic_v<T>) {
        denorm_min = std::numeric_limits<T>::denorm_min(); // if T is arithmetic, use its denorm_min
    }

    // compute unit nion
    if (phase2 <= denorm_min) return {mag2, phase2, unit}; // if unit is zero, return unit

    // otherwise, compute unit nion
    return {mag2, phase2, unit / sqrt(phase2)};
}

/******************************************
*            COMPARISON OPERATORS
*******************************************/

/**
* @brief evaluates if two values are similar to machine precision
* @tparam T type of first value
* @tparam S  type of second value
* @param a  first value
* @param b  second value
* @param epsilon tolerance
* @return true if similar, false otherwise
*/
template<arith_ops T, std::size_t N = NION_USE_HEAP, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
constexpr inline bool value_is_similar(T a, S b){
    constexpr T epsilon = std::numeric_limits<T>::epsilon();

    if constexpr (!std::is_arithmetic_v<T> || !std::is_arithmetic_v<S>)
        return abs(a - b) <= epsilon;
    return std::fabs(a - b) <= epsilon;
}

template<arith_ops T, std::size_t N>
template<std::size_t M>
constexpr inline bool nion<T,N>::operator==(const nion<T,M> &other) const {
    if (this == &other) return true;
    if (this->size() != other.size()) return false;

    D size = this->size();
    for (D i = 0; i < size; i++)
        if (!value_is_similar(elem_[i], other.elem_[i])) return false;
    return true;
}

template<arith_ops T, std::size_t N>
template<std::size_t M>
constexpr inline bool nion<T,N>::operator!=(const nion<T,M> &other) const {
    if (this == &other) return false;
    if (this->size() != other.size()) return true;

    D size = this->size();
    for (D i = 0; i < size; i++)
        if (!value_is_similar(elem_[i], other.elem_[i])) return true;
    return false;
}

template<arith_ops T, std::size_t N>
constexpr inline T nion<T,N>::rotate_real() const {
    // this yields the shortest rotation of the nion onto the real axis while preserving the norm and the sign of the real component
    // this is useful for comparing nions with different degrees
    return copysign(norm(), real());
}

template<arith_ops T, std::size_t N>
template<std::size_t M>
constexpr inline bool nion<T,N>::operator>(const nion<T,M> &other) const {
    // nions with degree > 1 are not ordered, but we can arbitrarily order them by their magnitude
    return abs() > other.abs();
}

template<arith_ops T, std::size_t N>
template<std::size_t M>
constexpr inline bool nion<T,N>::operator<(const nion<T,M> &other) const {
    // nions with degree > 1 are not ordered, but we can arbitrarily order them by their magnitude
    return abs() < other.abs();
}

template<arith_ops T, std::size_t N>
template<std::size_t M>
constexpr inline bool nion<T,N>::operator>=(const nion<T,M> &other) const {
    // nions with degree > 1 are not ordered, but we can arbitrarily order them by their magnitude
    if (abs() > other.abs())
        return true;
    return *this == *other;
}

template<arith_ops T, std::size_t N>
template<std::size_t M>
constexpr inline bool nion<T,N>::operator<=(const nion<T,M> &other) const {
    // nions with degree > 1 are not ordered, but we can arbitrarily order them by their magnitude
    if (abs() < other.abs())
        return true;
    return *this == *other;
}

template<arith_ops T, std::size_t N>
constexpr inline bool nion<T,N>::is_real() const {
    D size = this->size();
    for (D i = 1; i < size; i++)
        if (!value_is_similar(elem_[i], 0)) return false;
    return true;
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline bool nion<T,N>::operator==(S scalar) const {
    if (!value_is_similar(real(), scalar)) return false;
    return is_real();
}

template<arith_ops T, std::size_t N, not_nion<T,N> S>
constexpr inline bool operator==(S scalar, const nion<T,N> &z) {
    return z == scalar;
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline bool nion<T,N>::operator!=(S scalar) const {
    if (!value_is_similar(real(), scalar)) {
        return true;
    }
    return !is_real();
}

template<arith_ops T, std::size_t N, not_nion<T,N> S>
constexpr inline bool operator!=(S scalar, const nion<T,N> &z) {
    return z != scalar;
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline bool nion<T,N>::operator>(S scalar) const{
    return abs() > scalar;
}

template<arith_ops T, std::size_t N, not_nion<T,N> S>
constexpr inline bool operator>(S scalar, const nion<T,N> &z) {
    return z < scalar;
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline bool nion<T,N>::operator<(S scalar) const{
    return abs() < scalar;
}

template<arith_ops T, std::size_t N, not_nion<T,N> S>
constexpr inline bool operator<(S scalar, const nion<T,N> &z) {
    return z > scalar;
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline bool nion<T,N>::operator>=(S scalar) const{
    if (*this == nion<T,N>(scalar, this->size()))
        return true;
    return abs() > scalar;
}

template<arith_ops T, std::size_t N, not_nion<T,N> S>
constexpr inline bool operator>=(S scalar, const nion<T,N> &z) {
    return z <= scalar;
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline bool nion<T,N>::operator<=(S scalar) const{
    if (*this == nion<T,N>(scalar, this->size()))
        return true;
    return abs() < scalar;
}

template<arith_ops T, std::size_t N, not_nion<T,N> S>
constexpr inline bool operator<=(S scalar, const nion<T,N> &z) {
    return z >= scalar;
}

/******************************************
*            STREAMING OPERATORS
*******************************************/

template<arith_ops T, std::size_t N>
std::ostream &operator<<(std::ostream &os, const nion<T,N> &z) {
    T component = z.elem_[0];
    os << "(" << component;
    using D = typename nion<T,N>::D;
    D size = z.size();
    for (D i = 1; i < size; i++) {
        component = z.elem_[i];
        os << "," << component;
    }
    os << ")";
    return os;
}

template<arith_ops T, std::size_t N>
std::istream &operator>>(std::istream &is, nion<T,N> &z) {
    using D = typename nion<T,N>::D;
    D size = z.size();
    for (D i = 0; i < size; i++) {
        is >> z.elem_[i];
    }
    return is;
}

/**
* @brief Converts a nion to a string.
* @return The string representation of the nion.
*/
template<arith_ops T, std::size_t N>
inline std::string nion<T,N>::to_string() const {
    std::string nion_string = "(" + std::to_string(elem_[0]);
    D size = this->size();
    for (D i = 1; i < size; i++) {
        nion_string += "," + std::to_string(elem_[i]);
    }
    nion_string += ")";
    return nion_string;
}

/*********************************
*  NION FUNCTION IMPLEMENTATIONS *
**********************************/

template<arith_ops T, std::size_t N>
constexpr inline T real(const nion<T,N> &z) {
    return z.real();
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> imag(const nion<T,N> &z) {
    return z.imag();
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> conj(const nion<T,N> &z) {
    return z.conj();
}

template<arith_ops T, std::size_t N>
constexpr inline T abs(const nion<T,N> &z) {
    return z.abs();
}

template<arith_ops T, std::size_t N>
constexpr inline T norm(const nion<T,N> &z) {
    return z.norm();
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> inv(const nion<T,N> &z) {
    return z.inv();
}

template<arith_ops T, std::size_t N, std::size_t M>
constexpr inline T dot(const nion<T,N> &lhs, const nion<T,M> &rhs) {
    T dotProduct = 0;

    // find common type
    using D = typename nion<T,N>::D;
    using E = typename nion<T,M>::D;
    using uint_com = std::common_type_t<D, E>;

    // find the smaller degree
    uint_com minDegree = std::min(lhs.size(), rhs.size());

    // compute the dot product
    for (uint_com i = 0; i < minDegree; i++)
        dotProduct += lhs.elem_[i] * rhs.elem_[i];
    // the rest of the elements are zero for one of the nions if any

    return dotProduct;
}

/****************************
*  NION ALGEBRAIC FUNCTIONS *
*****************************/

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> exp(const nion<T,N> &z) {

    // get the polar form of nion
    T r = z.real();
    nion<T,N> i = z.imag();

    // make unit vector
    T i_abs = i.abs();
    T i_norm = sqrt(i_abs);

    // compute exponential of nion
    T cos_theta;
    T sin_theta;
    T exp_r = exp(r);

    // compute denorm_min
    T denorm_min = T(0); // default value
    if constexpr (std::is_arithmetic_v<T>) {
        denorm_min = std::numeric_limits<T>::denorm_min(); // if T is arithmetic, use its denorm_min
    }

    if (i_abs < denorm_min)
        return nion<T,N>(exp_r, z.size());

    // compute exponential of nion
    cos_theta = cos(i_norm);
    sin_theta = sin(i_norm);

    return i*(exp_r * sin_theta / i_norm) + exp_r * cos_theta;

}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> log(const nion<T,N> &z) {

    // get the polar form of nion
    T r = z.real();
    nion<T,N> i = z.imag();

    // compute norms
    T z_abs = z.abs();
    T z_norm = sqrt(z_abs);
    T i_abs = z_abs - r * r;
    T i_norm = sqrt(i_abs);
    T theta = atan2(i_norm, r);

    // compute denorm_min
    T denorm_min = T(0); // default value
    if constexpr (std::is_arithmetic_v<T>) {
        denorm_min = std::numeric_limits<T>::denorm_min(); // if T is arithmetic, use its denorm_min
    }

    // compute natural logarithm of nion
    if (i_abs <= denorm_min)
        return log(z_norm) + i * theta;
    else
        return i * (theta / i_norm) + log(z_norm);
}

template<arith_ops T, std::size_t N, not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> pow(const nion<T,N> &base, S power) {
    // get the polar form of nion
    T r = base.real();
    nion<T,N> i = base.imag();

    // compute norms
    T z_abs = base.abs();
    T z_norm = sqrt(z_abs);
    T i_abs = z_abs - r * r;
    T i_norm = sqrt(i_abs);

    T power_t;
    T theta;

    T cos_ptheta;
    T sin_ptheta;

    // compute denorm_min
    T denorm_min = T(0); // default value
    if constexpr (std::is_arithmetic_v<T>) {
        denorm_min = std::numeric_limits<T>::denorm_min(); // if T is arithmetic, use its denorm_min
    }

    if (i_abs >= denorm_min) i /= i_norm;
    else i = nion<T,N>(T(), base.size());

    if constexpr (std::is_integral_v<S>) { // if power is integer, use faster algorithm
        nion<T,N> z = base;
        if (std::signbit(power)) {
            z = inv(z);
            power = -power;
        }

        switch (power) {
            case 0:
                return nion<T,N>(1, z.size());
            case 1:
                return z;
            case 2:
                return sqr(z);
            default:
                power_t = static_cast<T>(power);
                theta = power_t * atan2(i_norm, r);
                cos_ptheta = cos(theta);
                sin_ptheta = sin(theta);
                break;
        }
    } else {
        power_t = static_cast<T>(power);
        theta = power_t * atan2(i_norm, r);
        cos_ptheta = cos(theta); //cos(atan(y/x)) = 1/sqrt(1+y^2/x^2) --> cos(p*atan(y/x)) = ???
        sin_ptheta = sin(theta); //sin(atan(y/x)) = y/(x*sqrt(1 + y^2/x^2)) --> sin(p*atan(y/x)) = ???
    }

    // compute power of nion
    auto result = pow(z_norm, power_t) * (cos_ptheta + i * (sin_ptheta));
    return result;
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> pow(const nion<T,N> &base, const nion<T,N> &power) {
    return exp(log(base) * power);
}

template<arith_ops T, std::size_t N, not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> pow(S base, const nion<T,N> &power) {
    return exp(log(base) * power);
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> sqr(const nion<T,N> &base) {
    // get the polar form of nion
    T r = base.real();
    nion<T,N> i = base.imag();

    // compute norms
    T base_abs = base.abs();
    T base_norm = sqrt(base_abs);

    // make unit vector
    T i_abs = base_abs - r * r;
    T i_norm = sqrt(i_abs);

    T power_t = 2.0l;

    // compute denorm_min
    T denorm_min = T(0); // default value
    if constexpr (std::is_arithmetic_v<T>) {
        denorm_min = std::numeric_limits<T>::denorm_min(); // if T is arithmetic, use its denorm_min
    }

    T x2 = r * r;
    T y2 = i_abs;
    if (x2 + y2 <= denorm_min)
        return nion<T,N>(base_norm * base_norm, base.size()); // if base is zero return zero (0^2 = 0)

    T denom = 1.0l / (x2 + y2);
    T cos_2theta = (x2 - y2) * denom;
    T sin_2theta = 2.0l * r * i_norm * denom;

    // compute power of nion
    if (i_abs <= denorm_min)
        return pow(base_norm, power_t) * (cos_2theta + i * sin_2theta);
    else
        return pow(base_norm, power_t) * (cos_2theta + i * (sin_2theta / i_norm));
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> sqrt(const nion<T,N> &z) {
    return pow(z, 0.5l);
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> cbrt(const nion<T,N> &z) {
    return pow(z, 1.0l / 3.0l);
}

// trigonometric functions
#include "nion_trig.hpp"

#endif //NION_CPP