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

#include <cmath>
#include <limits>
#include <type_traits>
#include <iostream>
#include <stdexcept>
#include "nion.hpp"


#define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << message << ": file=" << __FILE__ \
                      << ", line=" << __LINE__ << std::endl; \
            std::terminate(); \
        } \
    } while (false)

    /***************************
    *  NION CONSTRUCTORS
    ***************************/

    template<arith_ops T, std::size_t N>
    template <arith_ops S> requires (std::is_convertible_v<S, T>)
    constexpr inline nion<T,N>::nion(const S* vals, size_t size) : size_(size) {

        /// check if the degree is greater than zero and less than the maximum size
        ASSERT(size_ > 0, "The size of the nion must be greater than zero.");
        ASSERT(size_ <= N, "The size of the nion is too large. "
                                          "consider increasing template parameter, N.");

        /// copy the values into the nion
        if constexpr (std::same_as<T, S>) {
            memcpy(elem_, vals, size_ * sizeof(T));
        } else {
            for (size_t i = 0; i < size_; ++i) {
                elem_[i] = static_cast<T>(vals[i]);
            }
        }
    }

    template<arith_ops T, std::size_t N>
    constexpr inline nion<T,N>::nion(const std::initializer_list<T> &vals) : size_(vals.size()) {
        /// check if the degree is greater than zero
        ASSERT(size_ > 0, "The size of the nion must be greater than zero.");
        ASSERT(size_ <= N, "The size of the nion is too large. "
                                          "consider increasing template parameter, N.");

        /// copy the values into the nion
        memcpy(elem_, vals.begin(), size_ * sizeof(T));
    }

    template<arith_ops T, std::size_t N>
    constexpr inline nion<T,N>::nion(D size) : size_(size) {
        /// check if the degree is greater than zero
        ASSERT(size_ > 0, "The size of the nion must be greater than zero.");
        ASSERT(size_ <= N, "The size of the nion is too large. "
                                          "consider increasing template parameter, N.");

        /// initialize the values to zero
        zero();
    }

    template<arith_ops T, std::size_t N>
    template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S>)
    constexpr inline nion<T,N>::nion(const nion<S,M> &other) : size_(other.size_) {
        ASSERT(other.size_ <= N, "The size of the nion to copy is too large. "
                                 "consider increasing template parameter, N.");
        /// copy the values into the nion
        if constexpr (std::is_same_v<T,S>) {
            memcpy(elem_, other.elem_, size_ * sizeof(T));
        }
        else {
            for (int i = 0; i < size_; ++i)
                elem_[i] = other.elem_[i];
        }
    }

    template<arith_ops T, std::size_t N>
    template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S>)
    constexpr inline nion<T,N>::nion(nion<S,M> &&other) noexcept: size_(other.size_) {
        ASSERT(other.size_ <= N, "The size of the nion to copy is too large. "
                                 "consider increasing template parameter, N.");
        /// copy the values into the nion
        if constexpr (std::is_same_v<T,S>) {
            memcpy(elem_, other.elem_, size_ * sizeof(T));
        }
        else{
            for (int i = 0; i < size_; ++i) {
                elem_[i] = other.elem_[i];
            }
        }
    }

    template<arith_ops T, std::size_t N>
    template<not_nion<T,N> S> requires (std::is_convertible_v<S,T> && !std::is_pointer_v<S>)
    constexpr inline nion<T,N>::nion(S realVal, D size) : size_(size) {
        // check if the degree is greater than zero
        ASSERT(size_ > 0, "The degree of the nion must be greater than zero.");
        ASSERT(size_ <= N, "The size of the nion is too large. "
                                          "consider increasing template parameter, N.");

        // initialize the values to zero
        zero();

        // set the real part
        elem_[0] = realVal;
    }


    template<arith_ops T, std::size_t N>
    template<std::size_t M, std::size_t P>
    constexpr inline nion<T,N> nion<T,N>::make_pair(const nion<T,M> &a, const nion<T,P> &b) {

        /// initialize the nion pair
        nion<T,N> pair;

        /// set the size of the nion
        pair.size_ = a.size_ + b.size_;

        ASSERT(a.size_ > 0 && b.size_ > 0, "The sizes of the nion pair (a, b) must both be greater than zero.");
        ASSERT(pair.size_ <= N, "The size of the nion is too large. "
                                          "consider increasing template parameter, N.");

        /// copy the values into the nion
        memcpy(pair.elem_, a.elem_, a.size_ * sizeof(T));
        memcpy(pair.elem_ + a.size_, b.elem_, b.size_ * sizeof(T));

        return pair;
    }

    template<arith_ops T, std::size_t N>
    constexpr inline void nion<T,N>::resize(int size) {
        ASSERT(size > 0, "new nion size must be greater than zero");
        size_ = size;

        // set the new elements to zero
        memset(elem_ + size_, 0, (size - size_) * sizeof(T));
    }

    /************************************
    *         ASSIGNMENT OPERATORS
    *************************************/

    template<arith_ops T, std::size_t N>
    template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S>)
    constexpr inline nion<T,N> &nion<T,N>::operator=(const nion<S,M> &other) {
        /// check if the nions are the same
        if (reinterpret_cast<void*>(&other) == reinterpret_cast<void*>(this)) {
            return *this; // return the nion if they are the same
        }

        ASSERT(other.size_ <= N, "The size of the nion to assign is too large. "
                                 "consider increasing template parameter, N.");

        /// set the size
        size_ = other.size_;

        /// copy the values into the nion
        if constexpr (std::is_same_v<T,S>) {
            memcpy(elem_, other.elem_, size_ * sizeof(T));
        }
        else{
            for (int i = 0; i < size_; ++i) {
                elem_[i] = other.elem_[i];
            }
        }

        return *this; // return the nion
    }

    template<arith_ops T, std::size_t N>
    constexpr nion<T,N> &nion<T,N>::operator=(const std::initializer_list<T> &vals) {

        /// set the size of the nion
        size_ = vals.size();

        /// copy the values into the nion
        memcpy(elem_, vals.begin(), size_ * sizeof(T));
        return *this; // return the nion
    }

    template<arith_ops T, std::size_t N>
    template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S>)
    constexpr inline nion<T,N> &nion<T,N>::operator=(nion<S,M> &&other)  noexcept {
        /// check if the nions are the same
        if (reinterpret_cast<void*>(&other) == reinterpret_cast<void*>(this)) {
            return *this; // return the nion if they are the same
        }

        ASSERT(other.size_ <= N, "The size of the nion to assign is too large. "
                                 "consider increasing template parameter, N.");

        /// set the size
        size_ = other.size_;

        /// copy the values into the nion
        if constexpr (std::is_same_v<T,S>) {
            memcpy(elem_, other.elem_, size_ * sizeof(T));
        }
        else{
            for (int i = 0; i < size_; ++i) {
                elem_[i] = other.elem_[i];
            }
        }

        return *this; // return the nion
    }

    template<arith_ops T, std::size_t N>
    template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
    constexpr inline nion<T,N> &nion<T,N>::operator=(S scalar) {
        /// check if the nion is initialized
        if (size_ <= 0) size_ = 1; // set the degree

        zero(); // set the nion to zero
        elem_[0] = scalar; // set the real component
        return *this; // return the nion
    }

    /************************************
    *  ASSIGNMENT AND ADDITION OPERATORS
    *************************************/

    template<arith_ops T, std::size_t N>
    template<std::size_t M>
    constexpr inline void nion<T,N>::operator+=(const nion<T,M> &other) {

        // find common type
        using E = typename nion<T,M>::D;
        using uint_com = std::common_type_t<D, E>;

        // add the components of the other nion to this nion.
        if (size_ >= other.size_) {
            for (uint_com i = 0; i < other.size_; i++)
                elem_[i] += other.elem_[i];
        } else {
            for (uint_com i = 0; i < size_; i++)
                elem_[i] += other.elem_[i];

            // copy the remaining values
            memcpy(elem_ + size_, other.elem_ + size_, (other.size_ - size_) * sizeof(T));

            // set the new size of the nion
            size_ = other.size_;
        }

    }

    template<arith_ops T, std::size_t N>
    template<std::size_t M>
    constexpr inline void nion<T,N>::operator-=(const nion<T,M> &other) {
        *this += -other;
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

    template<arith_ops T, std::size_t N>
    template<std::size_t M>
    constexpr inline nion<T,N> nion<T,N>::operator+(const nion<T,M> &other) const {

        // find common int type
        using E = typename nion<T,M>::D;
        using uint_com = std::common_type_t<D, E>;

        if (size_ >= other.size_) {
            // create a nion to store the difference
            nion<T,N> sum(*this);

            // add the components of the other nion and this nion.
            for (uint_com i = 0; i < other.size_; i++)
                sum.elem_[i] += other.elem_[i];
            return sum;
        } else {
            // create a nion to store the difference
            nion<T,N> sum(other);

            // add the components of the other nion and this nion.
            for (uint_com i = 0; i < size_; i++)
                sum.elem_[i] += elem_[i];

            return sum;
        }
    }

    template<arith_ops T, std::size_t N>
    template<std::size_t M>
    constexpr inline nion<T,N> nion<T,N>::operator-(const nion<T,M> &other) const {
        return *this + -other;
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
        return -z + scalar;
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
        for (D i = 0; i < size_; i++)
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
        for (D i = 0; i < size_; i++)
            elem_[i] /= scalar;
    }

    /******************************************
    *        MULTIPLICATION OPERATORS
    *******************************************/

    template<arith_ops T, std::size_t N>
    constexpr inline nion<T,N> nion<T,N>::conj() const {
        nion<T,N> conjugate(*this); // copy this nion

        // negate all components except the first
        if constexpr (std::is_arithmetic_v<T>) {
            for (D i = 1; i < size_; i++)
                // negate the component by multiplying by -1 (faster than negation, but only works for arithmetic types)
                conjugate.elem_[i] *= -1;
        } else {
            for (D i = 1; i < size_; i++)
                conjugate.elem_[i] = -conjugate.elem_[i]; // negate the component
        }

        return conjugate;
    }

    template<arith_ops T, std::size_t N>
    constexpr inline nion<T,N> nion<T,N>::operator-() const {
        nion<T,N> negated(*this); // copy this nion

        // negate all components
        // negate all components except the first
        if constexpr (std::is_arithmetic_v<T>) {
            for (D i = 0; i < size_; i++)
                // negate the component by multiplying by -1 (faster than negation, but only works for arithmetic types)
                negated.elem_[i] *= -1;
        } else {
            for (D i = 0; i < size_; i++)
                negated.elem_[i] = -negated.elem_[i]; // negate the component
        }

        return negated;
    }

    template<arith_ops T, std::size_t N>
    template<std::size_t M>
    constexpr inline nion<T,N> nion<T,N>::operator*(const nion<T,M> &other) const {

        if (size_ == other.size_) {
            nion<T,N> product(*this); // create a nion to store the product

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds" // Disable the array out of bounds warning
#elif defined(__INTEL_COMPILER)
            #pragma warning(push)
                #pragma warning(disable: 2259) // Disable the array out of bounds warning
#endif
            switch (size_) {
                case 1: // if this size is 1, then the product is just the scalar product.
                    product.elem_[0] = elem_[0] * other.elem_[0];
                    return product;

// if FULL_RECURSION is defined, then the compiler will use recursion to calculate the product of nions for any size.
#ifndef FULL_RECURSION
                case 2: // hard-coded complex product
                    product.elem_[0] = elem_[0] * other.elem_[0] - elem_[1] * other.elem_[1];
                    product.elem_[1] = elem_[1] * other.elem_[0] + elem_[0] * other.elem_[1];
                    return product;

                case 4: // hard-coded quaternion product
                    product.elem_[0] = elem_[0] * other.elem_[0] - elem_[1] * other.elem_[1] - elem_[2] * other.elem_[2] - elem_[3] * other.elem_[3];
                    product.elem_[1] = elem_[1] * other.elem_[0] + elem_[0] * other.elem_[1] - elem_[3] * other.elem_[2] + elem_[2] * other.elem_[3];
                    product.elem_[2] = elem_[2] * other.elem_[0] + elem_[3] * other.elem_[1] + elem_[0] * other.elem_[2] - elem_[1] * other.elem_[3];
                    product.elem_[3] = elem_[3] * other.elem_[0] - elem_[2] * other.elem_[1] + elem_[1] * other.elem_[2] + elem_[0] * other.elem_[3];
                    return product;

                case 8: // hard-coded octonion product ( I know, it's a bit much )
                    product.elem_[0] = elem_[0] * other.elem_[0] - elem_[1] * other.elem_[1] - elem_[2] * other.elem_[2] - elem_[3] * other.elem_[3] - elem_[4] * other.elem_[4] - elem_[5] * other.elem_[5] - elem_[6] * other.elem_[6] - elem_[7] * other.elem_[7];
                    product.elem_[1] = elem_[1] * other.elem_[0] + elem_[0] * other.elem_[1] - elem_[3] * other.elem_[2] + elem_[2] * other.elem_[3] - elem_[5] * other.elem_[4] + elem_[4] * other.elem_[5] + elem_[7] * other.elem_[6] - elem_[6] * other.elem_[7];
                    product.elem_[2] = elem_[2] * other.elem_[0] + elem_[3] * other.elem_[1] + elem_[0] * other.elem_[2] - elem_[1] * other.elem_[3] - elem_[6] * other.elem_[4] - elem_[7] * other.elem_[5] + elem_[4] * other.elem_[6] + elem_[5] * other.elem_[7];
                    product.elem_[3] = elem_[3] * other.elem_[0] - elem_[2] * other.elem_[1] + elem_[1] * other.elem_[2] + elem_[0] * other.elem_[3] - elem_[7] * other.elem_[4] + elem_[6] * other.elem_[5] - elem_[5] * other.elem_[6] + elem_[4] * other.elem_[7];
                    product.elem_[4] = elem_[4] * other.elem_[0] + elem_[5] * other.elem_[1] + elem_[6] * other.elem_[2] + elem_[7] * other.elem_[3] + elem_[0] * other.elem_[4] - elem_[1] * other.elem_[5] - elem_[2] * other.elem_[6] - elem_[3] * other.elem_[7];
                    product.elem_[5] = elem_[5] * other.elem_[0] - elem_[4] * other.elem_[1] + elem_[7] * other.elem_[2] - elem_[6] * other.elem_[3] + elem_[1] * other.elem_[4] + elem_[0] * other.elem_[5] + elem_[3] * other.elem_[6] - elem_[2] * other.elem_[7];
                    product.elem_[6] = elem_[6] * other.elem_[0] - elem_[7] * other.elem_[1] - elem_[4] * other.elem_[2] + elem_[5] * other.elem_[3] + elem_[2] * other.elem_[4] - elem_[3] * other.elem_[5] + elem_[0] * other.elem_[6] + elem_[1] * other.elem_[7];
                    product.elem_[7] = elem_[7] * other.elem_[0] + elem_[6] * other.elem_[1] - elem_[5] * other.elem_[2] - elem_[4] * other.elem_[3] + elem_[3] * other.elem_[4] + elem_[2] * other.elem_[5] - elem_[1] * other.elem_[6] + elem_[0] * other.elem_[7];
                    return product;

                case 16: // hard-coded sedenion product ( I strongly recommend using nowrap if you're looking at this... )
                    product.elem_[ 0] = elem_[ 0] * other.elem_[ 0] - elem_[ 1] * other.elem_[ 1] - elem_[ 2] * other.elem_[ 2] - elem_[ 3] * other.elem_[ 3] - elem_[ 4] * other.elem_[ 4] - elem_[ 5] * other.elem_[ 5] - elem_[ 6] * other.elem_[ 6] - elem_[ 7] * other.elem_[ 7] - elem_[ 8] * other.elem_[ 8] - elem_[ 9] * other.elem_[ 9] - elem_[10] * other.elem_[10] - elem_[11] * other.elem_[11] - elem_[12] * other.elem_[12] - elem_[13] * other.elem_[13] - elem_[14] * other.elem_[14] - elem_[15] * other.elem_[15];
                    product.elem_[ 1] = elem_[ 0] * other.elem_[ 1] + elem_[ 1] * other.elem_[ 0] + elem_[ 2] * other.elem_[ 3] - elem_[ 3] * other.elem_[ 2] + elem_[ 4] * other.elem_[ 5] - elem_[ 5] * other.elem_[ 4] - elem_[ 6] * other.elem_[ 7] + elem_[ 7] * other.elem_[ 6] + elem_[ 8] * other.elem_[ 9] - elem_[ 9] * other.elem_[ 8] - elem_[10] * other.elem_[11] + elem_[11] * other.elem_[10] - elem_[12] * other.elem_[13] + elem_[13] * other.elem_[12] + elem_[14] * other.elem_[15] - elem_[15] * other.elem_[14];
                    product.elem_[ 2] = elem_[ 0] * other.elem_[ 2] - elem_[ 1] * other.elem_[ 3] + elem_[ 2] * other.elem_[ 0] + elem_[ 3] * other.elem_[ 1] + elem_[ 4] * other.elem_[ 6] + elem_[ 5] * other.elem_[ 7] - elem_[ 6] * other.elem_[ 4] - elem_[ 7] * other.elem_[ 5] + elem_[ 8] * other.elem_[10] + elem_[ 9] * other.elem_[11] - elem_[10] * other.elem_[ 8] - elem_[11] * other.elem_[ 9] - elem_[12] * other.elem_[14] - elem_[13] * other.elem_[15] + elem_[14] * other.elem_[12] + elem_[15] * other.elem_[13];
                    product.elem_[ 3] = elem_[ 0] * other.elem_[ 3] + elem_[ 1] * other.elem_[ 2] - elem_[ 2] * other.elem_[ 1] + elem_[ 3] * other.elem_[ 0] + elem_[ 4] * other.elem_[ 7] - elem_[ 5] * other.elem_[ 6] + elem_[ 6] * other.elem_[ 5] - elem_[ 7] * other.elem_[ 4] + elem_[ 8] * other.elem_[11] - elem_[ 9] * other.elem_[10] + elem_[10] * other.elem_[ 9] - elem_[11] * other.elem_[ 8] - elem_[12] * other.elem_[15] + elem_[13] * other.elem_[14] - elem_[14] * other.elem_[13] + elem_[15] * other.elem_[12];
                    product.elem_[ 4] = elem_[ 0] * other.elem_[ 4] - elem_[ 1] * other.elem_[ 5] - elem_[ 2] * other.elem_[ 6] - elem_[ 3] * other.elem_[ 7] + elem_[ 4] * other.elem_[ 0] + elem_[ 5] * other.elem_[ 1] + elem_[ 6] * other.elem_[ 2] + elem_[ 7] * other.elem_[ 3] + elem_[ 8] * other.elem_[12] + elem_[ 9] * other.elem_[13] + elem_[10] * other.elem_[14] + elem_[11] * other.elem_[15] - elem_[12] * other.elem_[ 8] - elem_[13] * other.elem_[ 9] - elem_[14] * other.elem_[10] - elem_[15] * other.elem_[11];
                    product.elem_[ 5] = elem_[ 0] * other.elem_[ 5] + elem_[ 1] * other.elem_[ 4] - elem_[ 2] * other.elem_[ 7] + elem_[ 3] * other.elem_[ 6] - elem_[ 4] * other.elem_[ 1] + elem_[ 5] * other.elem_[ 0] - elem_[ 6] * other.elem_[ 3] + elem_[ 7] * other.elem_[ 2] + elem_[ 8] * other.elem_[13] - elem_[ 9] * other.elem_[12] + elem_[10] * other.elem_[15] - elem_[11] * other.elem_[14] + elem_[12] * other.elem_[ 9] - elem_[13] * other.elem_[ 8] + elem_[14] * other.elem_[11] - elem_[15] * other.elem_[10];
                    product.elem_[ 6] = elem_[ 0] * other.elem_[ 6] + elem_[ 1] * other.elem_[ 7] + elem_[ 2] * other.elem_[ 4] - elem_[ 3] * other.elem_[ 5] - elem_[ 4] * other.elem_[ 2] + elem_[ 5] * other.elem_[ 3] + elem_[ 6] * other.elem_[ 0] - elem_[ 7] * other.elem_[ 1] + elem_[ 8] * other.elem_[14] - elem_[ 9] * other.elem_[15] - elem_[10] * other.elem_[12] + elem_[11] * other.elem_[13] + elem_[12] * other.elem_[10] - elem_[13] * other.elem_[11] - elem_[14] * other.elem_[ 8] + elem_[15] * other.elem_[ 9];
                    product.elem_[ 7] = elem_[ 0] * other.elem_[ 7] - elem_[ 1] * other.elem_[ 6] + elem_[ 2] * other.elem_[ 5] + elem_[ 3] * other.elem_[ 4] - elem_[ 4] * other.elem_[ 3] - elem_[ 5] * other.elem_[ 2] + elem_[ 6] * other.elem_[ 1] + elem_[ 7] * other.elem_[ 0] + elem_[ 8] * other.elem_[15] + elem_[ 9] * other.elem_[14] - elem_[10] * other.elem_[13] - elem_[11] * other.elem_[12] + elem_[12] * other.elem_[11] + elem_[13] * other.elem_[10] - elem_[14] * other.elem_[ 9] - elem_[15] * other.elem_[ 8];
                    product.elem_[ 8] = elem_[ 0] * other.elem_[ 8] - elem_[ 1] * other.elem_[ 9] - elem_[ 2] * other.elem_[10] - elem_[ 3] * other.elem_[11] - elem_[ 4] * other.elem_[12] - elem_[ 5] * other.elem_[13] - elem_[ 6] * other.elem_[14] - elem_[ 7] * other.elem_[15] + elem_[ 8] * other.elem_[ 0] + elem_[ 9] * other.elem_[ 1] + elem_[10] * other.elem_[ 2] + elem_[11] * other.elem_[ 3] + elem_[12] * other.elem_[ 4] + elem_[13] * other.elem_[ 5] + elem_[14] * other.elem_[ 6] + elem_[15] * other.elem_[ 7];
                    product.elem_[ 9] = elem_[ 0] * other.elem_[ 9] + elem_[ 1] * other.elem_[ 8] - elem_[ 2] * other.elem_[11] + elem_[ 3] * other.elem_[10] - elem_[ 4] * other.elem_[13] + elem_[ 5] * other.elem_[12] + elem_[ 6] * other.elem_[15] - elem_[ 7] * other.elem_[14] - elem_[ 8] * other.elem_[ 1] + elem_[ 9] * other.elem_[ 0] - elem_[10] * other.elem_[ 3] + elem_[11] * other.elem_[ 2] - elem_[12] * other.elem_[ 5] + elem_[13] * other.elem_[ 4] + elem_[14] * other.elem_[ 7] - elem_[15] * other.elem_[ 6];
                    product.elem_[10] = elem_[ 0] * other.elem_[10] + elem_[ 1] * other.elem_[11] + elem_[ 2] * other.elem_[ 8] - elem_[ 3] * other.elem_[ 9] - elem_[ 4] * other.elem_[14] - elem_[ 5] * other.elem_[15] + elem_[ 6] * other.elem_[12] + elem_[ 7] * other.elem_[13] - elem_[ 8] * other.elem_[ 2] + elem_[ 9] * other.elem_[ 3] + elem_[10] * other.elem_[ 0] - elem_[11] * other.elem_[ 1] - elem_[12] * other.elem_[ 6] - elem_[13] * other.elem_[ 7] + elem_[14] * other.elem_[ 4] + elem_[15] * other.elem_[ 5];
                    product.elem_[11] = elem_[ 0] * other.elem_[11] - elem_[ 1] * other.elem_[10] + elem_[ 2] * other.elem_[ 9] + elem_[ 3] * other.elem_[ 8] - elem_[ 4] * other.elem_[15] + elem_[ 5] * other.elem_[14] - elem_[ 6] * other.elem_[13] + elem_[ 7] * other.elem_[12] - elem_[ 8] * other.elem_[ 3] - elem_[ 9] * other.elem_[ 2] + elem_[10] * other.elem_[ 1] + elem_[11] * other.elem_[ 0] - elem_[12] * other.elem_[ 7] + elem_[13] * other.elem_[ 6] - elem_[14] * other.elem_[ 5] + elem_[15] * other.elem_[ 4];
                    product.elem_[12] = elem_[ 0] * other.elem_[12] + elem_[ 1] * other.elem_[13] + elem_[ 2] * other.elem_[14] + elem_[ 3] * other.elem_[15] + elem_[ 4] * other.elem_[ 8] - elem_[ 5] * other.elem_[ 9] - elem_[ 6] * other.elem_[10] - elem_[ 7] * other.elem_[11] - elem_[ 8] * other.elem_[ 4] + elem_[ 9] * other.elem_[ 5] + elem_[10] * other.elem_[ 6] + elem_[11] * other.elem_[ 7] + elem_[12] * other.elem_[ 0] - elem_[13] * other.elem_[ 1] - elem_[14] * other.elem_[ 2] - elem_[15] * other.elem_[ 3];
                    product.elem_[13] = elem_[ 0] * other.elem_[13] - elem_[ 1] * other.elem_[12] + elem_[ 2] * other.elem_[15] - elem_[ 3] * other.elem_[14] + elem_[ 4] * other.elem_[ 9] + elem_[ 5] * other.elem_[ 8] + elem_[ 6] * other.elem_[11] - elem_[ 7] * other.elem_[10] - elem_[ 8] * other.elem_[ 5] - elem_[ 9] * other.elem_[ 4] + elem_[10] * other.elem_[ 7] - elem_[11] * other.elem_[ 6] + elem_[12] * other.elem_[ 1] + elem_[13] * other.elem_[ 0] + elem_[14] * other.elem_[ 3] - elem_[15] * other.elem_[ 2];
                    product.elem_[14] = elem_[ 0] * other.elem_[14] - elem_[ 1] * other.elem_[15] - elem_[ 2] * other.elem_[12] + elem_[ 3] * other.elem_[13] + elem_[ 4] * other.elem_[10] - elem_[ 5] * other.elem_[11] + elem_[ 6] * other.elem_[ 8] + elem_[ 7] * other.elem_[ 9] - elem_[ 8] * other.elem_[ 6] - elem_[ 9] * other.elem_[ 7] - elem_[10] * other.elem_[ 4] + elem_[11] * other.elem_[ 5] + elem_[12] * other.elem_[ 2] - elem_[13] * other.elem_[ 3] + elem_[14] * other.elem_[ 0] + elem_[15] * other.elem_[ 1];
                    product.elem_[15] = elem_[ 0] * other.elem_[15] + elem_[ 1] * other.elem_[14] - elem_[ 2] * other.elem_[13] - elem_[ 3] * other.elem_[12] + elem_[ 4] * other.elem_[11] + elem_[ 5] * other.elem_[10] - elem_[ 6] * other.elem_[ 9] + elem_[ 7] * other.elem_[ 8] - elem_[ 8] * other.elem_[ 7] + elem_[ 9] * other.elem_[ 6] - elem_[10] * other.elem_[ 5] - elem_[11] * other.elem_[ 4] + elem_[12] * other.elem_[ 3] + elem_[13] * other.elem_[ 2] - elem_[14] * other.elem_[ 1] + elem_[15] * other.elem_[ 0];
                    return product;
                //case 32: TODO: No way. You can't make me. 32x32 products are just too big. I could do it, but it wouldn't be pretty. I'll leave it to you.
#endif
            } // if size_ is not 1, 2, 8, or 16, we need to do the recursive product.

#ifdef __GNUC__
#pragma GCC diagnostic pop
#elif defined(__INTEL_COMPILER)
#pragma warning(pop)
#endif

        } else { // If the matrix is not square, we need to do the recursive product.

            // if any of the sizes are 1, then the product is just the scalar product.
            if (size_ == 1) {
                return other * elem_[0];
            }
            if (other.size_ == 1) {
                return *this * other.elem_[0];
            }
        }

        /// ********** Recursive Product ********** ///

        // find common int type
        using E = typename nion<T,M>::D;
        using uint_com = std::common_type_t<D, E>;

        // the elements of the first half of the nion
        uint_com this_half_size = size_ - (size_ >> 1);
        uint_com other_half_size = other.size_ - (other.size_ >> 1);

        // pointers to the halves of the nions
        const T *a_ptr = elem_,       *b_ptr = &elem_[this_half_size],
                *c_ptr = other.elem_, *d_ptr = &other.elem_[other_half_size];

        // construct halves of the nions
        nion<T,N-(N>>1)> a(a_ptr, this_half_size),  b(b_ptr, size_ - this_half_size);
        nion<T,M-(M>>1)> c(c_ptr, other_half_size), d(d_ptr, other.size_ - other_half_size);

        /// calculate the cayley-dickson product
        return make_pair(
                (a * c) - (d.conj() * b), // add involution parameter for sign with macro? (split hypercomplex numbers)
                (d * a) + (b * c.conj())
        );
    }

    template<arith_ops T, std::size_t N>
    template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
    constexpr inline nion<T,N> nion<T,N>::operator*(S scalar) const {
        nion<T,N> product;
        product.size_ = size_;

        // compute the product of each element of the nion with the scalar
        for (D i = 0; i < size_; i++)
            product.elem_[i] = elem_[i] * scalar;

        return product;
    }

    template<arith_ops T, std::size_t N, not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
    constexpr inline nion<T,N> operator*(S scalar, const nion<T,N> &z) {
        nion<T,N> product;
        product.size_ = z.size_;

        // compute the product of each element of the nion with the scalar
        using D = typename nion<T,N>::D;
        for (D i = 0; i < z.size_; i++)
            product.elem_[i] = scalar * z.elem_[i];

        return product;
        return z * scalar;
    }

    template<arith_ops T, std::size_t N>
    constexpr inline T nion<T,N>::abs() const {

        T absVal = elem_[0] * elem_[0];

        for (D i = 1; i < size_; i++)
            absVal += elem_[i] * elem_[i];

        return absVal;
    }

    template<arith_ops T, std::size_t N>
    constexpr inline T nion<T,N>::norm() const {
        return sqrt(abs());
    }

    template<arith_ops T, std::size_t N>
    constexpr inline nion<T,N> nion<T,N>::inv() const {

        T absolute = abs();
        nion<T,N> inverse;
        inverse.size_ = size_;

        inverse.elem_[0] = elem_[0] / absolute;
        for (D i = 1; i < size_; i++)
            inverse.elem_[i] = -elem_[i] / absolute;

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
        nion<T,N> quotient;
        quotient.size_ = size_;

        // compute the product of each element of the nion with the scalar
        for (D i = 0; i < size_; i++)
            quotient.elem_[i] = elem_[i] / scalar;

        return quotient;
    }

    template<arith_ops T, std::size_t N, not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
    constexpr inline nion<T,N> operator/(S scalar, const nion<T,N> &z) {
        nion<T,N> quotient;
        quotient.size_ = z.size_;

        // compute the product of each element of the nion with the scalar
        using D = typename nion<T,N>::D;
        for (D i = 0; i < z.size_; i++)
            quotient.elem_[i] = (1.0l / scalar) * z.elem_[i] ;

        return quotient;
    }


    /******************************************
    *            ACCESSOR FUNCTIONS
    *******************************************/

    template<arith_ops T, std::size_t N>
    constexpr inline T nion<T,N>::operator[](D index) const {
        return this->elem_[index];
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
        // get polar form of nion
        T r = elem_[0];

        // compute norms
        T z_abs = abs();
        T i_norm = sqrt(z_abs - r * r);

        // compute argument
        T theta = atan2(i_norm, r);
        return theta;
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
    template<arith_ops T, std::size_t N = 128, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline bool value_is_similar(T a, S b){
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        return std::fabs(a - b) <= epsilon;
    }

    template<arith_ops T, std::size_t N>
    template<std::size_t M>
    constexpr inline bool nion<T,N>::operator==(const nion<T,M> &other) const {
        if (this == &other) return true;
        if (size_ != other.size_) return false;

        for (D i = 0; i < size_; i++)
            if (!value_is_similar(elem_[i], other.elem_[i])) return false;
        return true;
    }

    template<arith_ops T, std::size_t N>
    template<std::size_t M>
    constexpr inline bool nion<T,N>::operator!=(const nion<T,M> &other) const {
        if (this == &other) return false;
        if (size_ != other.size_) return true;

        for (D i = 0; i < size_; i++)
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
        // nions with degree > 1 are not ordered, but we can arbitrarily order them by their rotation to the real line
        // this is not a good idea, but it's better than nothing (and it's what I'm doing for now)
        return rotate_real() > other.rotate_real();
    }

    template<arith_ops T, std::size_t N>
    template<std::size_t M>
    constexpr inline bool nion<T,N>::operator<(const nion<T,M> &other) const {
        // nions with degree > 1 are not ordered, but we can arbitrarily order them by their rotation to the real line
        // this is not a good idea, but it's better than nothing (and it's what I'm doing for now)
        return rotate_real() < other.rotate_real();
    }

    template<arith_ops T, std::size_t N>
    template<std::size_t M>
    constexpr inline bool nion<T,N>::operator>=(const nion<T,M> &other) const {
        if (rotate_real() > other.rotate_real())
            return true;
        return *this == *other;
    }

    template<arith_ops T, std::size_t N>
    template<std::size_t M>
    constexpr inline bool nion<T,N>::operator<=(const nion<T,M> &other) const {
        if (rotate_real() < other.rotate_real())
            return true;
        return *this == *other;
    }

    template<arith_ops T, std::size_t N>
    constexpr inline bool nion<T,N>::is_real() const {
        for (D i = 1; i < size_; i++)
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
        return rotate_real() > scalar;
    }

    template<arith_ops T, std::size_t N, not_nion<T,N> S>
    constexpr inline bool operator>(S scalar, const nion<T,N> &z) {
        return z < scalar;
    }

    template<arith_ops T, std::size_t N>
    template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
    constexpr inline bool nion<T,N>::operator<(S scalar) const{
        return rotate_real() < scalar;
    }

    template<arith_ops T, std::size_t N, not_nion<T,N> S>
    constexpr inline bool operator<(S scalar, const nion<T,N> &z) {
        return z > scalar;
    }

    template<arith_ops T, std::size_t N>
    template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
    constexpr inline bool nion<T,N>::operator>=(S scalar) const{
        if (*this == nion<T,N>(scalar, this->size_))
            return true;
        return rotate_real() > scalar;
    }

    template<arith_ops T, std::size_t N, not_nion<T,N> S>
    constexpr inline bool operator>=(S scalar, const nion<T,N> &z) {
        return z <= scalar;
    }

    template<arith_ops T, std::size_t N>
    template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
    constexpr inline bool nion<T,N>::operator<=(S scalar) const{
        if (*this == nion<T,N>(scalar, this->size_))
            return true;
        return rotate_real() < scalar;
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
        for (D i = 1; i < z.size_; i++) {
            component = z.elem_[i];
            os << "," << component;
        }
        os << ")";
        return os;
    }

    template<arith_ops T, std::size_t N>
    std::istream &operator>>(std::istream &is, nion<T,N> &z) {
        using D = typename nion<T,N>::D;
        for (D i = 0; i < z.size_; i++) {
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
        for (D i = 1; i < size_; i++) {
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
        uint_com minDegree = std::min(lhs.size_, rhs.size_);

        // compute the dot product
        for (uint_com i = 0; i < minDegree; i++)
            dotProduct += lhs.elem_[i] * rhs.elem_[i];

        return dotProduct;
    }

    /*
    template<arith_ops T, std::size_t N>
    constexpr inline nion<T,N>
    cross(const nion<T,N> &lhs, const nion<T,N> &rhs){} //TODO: implement cross product

    template<arith_ops T, std::size_t N>
    constexpr inline nion<T,N>
    wedge(const nion<T,N> &lhs, const nion<T,N> &rhs){} //TODO: implement wedge product

    template<arith_ops T, std::size_t N>
    constexpr inline nion<T,N>
    outer(const nion<T,N> &lhs, const nion<T,N> &rhs){} //TODO: implement outer product
     */


    /****************************
    *  NION ALGEBRAIC FUNCTIONS *
    *****************************/

    template<arith_ops T, std::size_t N>
    constexpr inline nion<T,N> exp(const nion<T,N> &z) {

        // get polar form of nion
        T r = z.real();
        nion<T,N> i = z.imag();

        // make unit vector
        T i_abs = i.abs();
        T i_norm = sqrt(i_abs);

        // compute exponential of nion
        T cos_theta;
        T sin_theta;
        T exp_r = exp(r);

        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
        if (i_abs < denorm_min)
            return nion<T,N>(exp_r, z.size_);

        // compute exponential of nion
        cos_theta = cos(i_norm);
        sin_theta = sin(i_norm);

        return i*(exp_r * sin_theta / i_norm) + exp_r * cos_theta;

    }

    template<arith_ops T, std::size_t N>
    constexpr inline nion<T,N> log(const nion<T,N> &z) {

        // get polar form of nion
        T r = z.real();
        nion<T,N> i = z.imag();

        // compute norms
        T z_abs = z.abs();
        T z_norm = sqrt(z_abs);
        T i_abs = z_abs - r * r;
        T i_norm = sqrt(z_abs - r * r);
        T theta = atan2(i_norm, r);

        // compute natural logarithm of nion
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
        if (i_abs <= denorm_min)
            return log(z_norm) + i * theta;
        else
            return i * (theta / i_norm) + log(z_norm);
    }

    template<arith_ops T, std::size_t N, not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
    constexpr inline nion<T,N> pow(const nion<T,N> &base, S power) {
        return exp(log(base) * power);
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
        // get polar form of nion
        T r = base.real();
        nion<T,N> i = base.imag();

        // compute norms
        T base_abs = base.abs();
        T base_norm = sqrt(base_abs);

        // make unit vector
        T i_abs = base_abs - r * r;
        T i_norm = sqrt(i_abs);

        T power_t = 2.0l;
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();

        T x2 = r * r;
        T y2 = i_abs;
        if (x2 + y2 <= denorm_min)
            return nion<T,N>(base_norm * base_norm, base.size_); // if base is zero return zero (0^2 = 0)

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

        // compute exponential of nion
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
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

        // compute exponential of nion
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
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

        // compute the sine of the nion
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
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

        // compute the cosine of the nion
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
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

        // compute the tangent of the nion
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
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
    constexpr inline nion<T,N> asin(const nion<T,N> &z) {

        // get the polar form of the nion
        T r = real(z);
        nion<T,N> i = z.imag();

        // make unit vector
        T i_abs = i.abs();
        T i_norm = sqrt(i_abs);
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
        if (i_abs > denorm_min)
            i /= i_norm;

        // compute the inv sine of the nion
        return -i * log(sqrt(1.0l - sqr(z)) + (i * r) - i_norm);
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
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
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

