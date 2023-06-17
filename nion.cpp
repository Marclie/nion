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

#ifndef NION_CPP
#define NION_CPP

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
constexpr inline nion<T,N>::nion(const S* vals, std::size_t size) : size_(size) {

    /// check if the degree is greater than zero and less than the maximum size
    ASSERT(size_ > 0, "The size of the nion must be greater than zero.");
    ASSERT(size_ <= N || N == 0,
           "The size of the nion is too large. "
           "consider increasing the template parameter, N, or using the default value for heap allocation."
           );

    if constexpr (on_heap) // if the user wants to use the heap
        elem_ = new T[size_]; // allocate memory on the heap
    // else, the user wants to use the stack

    // copy the values into the nion
    for (std::size_t i = 0; i < size_; ++i) {
        elem_[i] = static_cast<T>(vals[i]);
    }

}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N>::nion(const std::initializer_list<T> &vals) : size_(vals.size()) {
    /// check if the degree is greater than zero
    ASSERT(size_ > 0, "The size of the nion must be greater than zero.");
    ASSERT(size_ <= N || N == 0,
           "The size of the nion is too large. "
           "consider increasing the template parameter, N, or using the default value for heap allocation."
    );

    if constexpr (on_heap) // if the user wants to use the heap
        elem_ = new T[size_]; // allocate memory on the heap
    // else, the user wants to use the stack

    // copy the values into the nion
    for (std::size_t i = 0; i < size_; ++i) {
        elem_[i] = *(vals.begin() + i);
    }
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N>::nion(D size) : size_(size) {
    /// check if the degree is greater than zero
    ASSERT(size_ > 0, "The size of the nion must be greater than zero.");
    ASSERT(size_ <= N || N == 0,
           "The size of the nion is too large. "
           "consider increasing the template parameter, N, or using the default value for heap allocation."
    );

    if constexpr (on_heap) // if the nion is on the heap
        elem_ = new T[size_]; // allocate memory for the nion

    /// initialize the values to zero
    zero();
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N>::nion(const nion<T,N> &other) : size_(other.size_) {
    ASSERT(size_ <= N || N == 0,
           "The size of the nion is too large. "
           "consider increasing the template parameter, N, or using the default value for heap allocation."
    );

    if constexpr (on_heap) // if the user wants to use the heap
        elem_ = new T[size_]; // allocate memory on the heap
    // else, the user wants to use the stack

    // copy the values into the nion
    for (std::size_t i = 0; i < size_; ++i) {
        elem_[i] = other.elem_[i];
    }
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> &nion<T,N>::operator=(const nion<T,N> &other) {
    // check if the nions are the same
    if (&other == this) return *this; // return the nion if they are the same
    
    // set the size
    resize(other.size_); // resize the nion

    // copy the values into the nion
    for (std::size_t i = 0; i < size_; ++i) {
        elem_[i] = other.elem_[i];
    }
    return *this; // return the nion
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N>::nion(nion<T,N> &&other) noexcept : size_(other.size_) {
    ASSERT(size_ <= N || N == 0,
           "The size of the nion is too large. "
           "consider increasing the template parameter, N, or using the default value for heap allocation."
    );

    if constexpr (on_heap) { // if the user wants to use the heap
        elem_ = other.elem_; // allocate memory on the heap
        other.elem_ = nullptr;
        return;
    } else { // the user wants to use the stack, which requires a deep copy
        for (std::size_t i = 0; i < size_; ++i) {
            elem_[i] = other.elem_[i];
        }
    }
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> &nion<T,N>::operator=(nion<T,N> &&other)  noexcept {
    /// check if the nions are the same
    if (&other == this) return *this; // return the nion if they are the same

    /// set the size
    size_ = other.size_;

    if constexpr (on_heap) {// if the user wants to use the heap
        // just copy the pointer and delete this one
        if (elem_) delete[] elem_; // free the memory if it is not null
        elem_ = other.elem_;
        other.elem_ = nullptr;
    } else { // the user wants to use the stack, and we need to deep copy the values

        ASSERT(size_ <= N || N == 0,
               "The size of the nion is too large. "
               "consider increasing the template parameter, N, or using the default value for heap allocation."
        );

        // copy the values into the nion
        for (std::size_t i = 0; i < size_; ++i) {
            elem_[i] = other.elem_[i];
        }
    }
    return *this; // return the nion
}

template<arith_ops T, std::size_t N>
template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S> && (M != N || !std::is_same_v<T, S>))
constexpr inline nion<T,N>::nion(const nion<S,M> &other) : size_(other.size_) {
    ASSERT(other.size_ <= N || N == 0, "The size of the nion to copy is too large. "
                             "consider increasing the template parameter, N, or using the default value for heap allocation.");

    if constexpr (on_heap) // if the user wants to use the heap
        elem_ = new T[size_]; // allocate memory on the heap
    // else, the user wants to use the stack

    // copy the values into the nion
    for (D i = 0; i < size_; ++i)
        elem_[i] = static_cast<T>(other.elem_[i]);

}

template<arith_ops T, std::size_t N>
template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S> && (M != N || !std::is_same_v<T, S>))
constexpr inline nion<T,N> &nion<T,N>::operator=(const nion<S,M> &other) {

    // resize the nion (if needed)
    resize(other.size_);

    // copy the values into the nion
    for (D i = 0; i < size_; ++i)
        elem_[i] = static_cast<T>(other.elem_[i]);

    return *this; // return the nion
}

template<arith_ops T, std::size_t N>
template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S> && (M != N || !std::is_same_v<T, S>))
constexpr inline nion<T,N>::nion(nion<S,M> &&other) noexcept: size_(other.size_) {
    ASSERT(other.size_ <= N || N == 0, "The size of the nion to copy is too large. "
                                       "consider increasing the template parameter, N, or using the default value for heap allocation.");

    if constexpr (on_heap && nion<S,M>::on_heap) { // if the user wants to use the heap
        if (elem_) delete[] elem_; // free the memory if it is not null
        elem_ = other.elem_; // allocate memory on the heap
        other.elem_ = nullptr;
    } else { //else the user wants to use the stack, which requires a copy
        *this = other; // copy the other nion
    }
}

template<arith_ops T, std::size_t N>
template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S> && (M != N || !std::is_same_v<T, S>))
constexpr inline nion<T,N> &nion<T,N>::operator=(nion<S,M> &&other)  noexcept {

    ASSERT(other.size_ <= N || N == 0, "The size of the nion to assign is too large. "
                                       "consider increasing the template parameter, N, or using the default value for heap allocation.");

    /// set the size
    size_ = other.size_;

    if constexpr (on_heap && nion<S,M>::on_heap) {// if the user wants to use the heap
        if (elem_) delete[] elem_; // free the memory if it is not null
        elem_ = other.elem_; // allocate memory on the heap
        other.elem_ = nullptr;
        return *this;
    } //else the user wants to use the stack, which requires a copy
    *this = other; // copy the other nion

    return *this; // return the nion
}



/************************************
*         OTHER CONSTRUCTORS        *
*************************************/

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T> && !std::is_pointer_v<S>)
constexpr inline nion<T,N>::nion(S realVal, D size) : size_(size) {
    // check if the degree is greater than zero
    ASSERT(size_ > 0, "The degree of the nion must be greater than zero.");
    ASSERT(size_ <= N || N == 0, "The size of the nion is too large. "
                                      "consider increasing the template parameter, N, or using the default value for heap allocation.");

    if constexpr (on_heap) // if the nion is on the heap
        elem_ = new T[size_]; // allocate memory for the nion

    // initialize the values to zero
    zero();

    // set the real part
    elem_[0] = static_cast<T>(realVal);
}

template<arith_ops T, std::size_t N>
template<arith_ops S, std::size_t M>
requires(std::is_convertible_v<T, S>)
constexpr nion<T,N>::operator nion<S,M>() {
    if constexpr (std::is_same_v<T, S>) // if types are the same, just copy the array
        return nion<T,M>(elem_, size_);
    else if constexpr (N == M) // if max sizes are the same, just copy the array
        return nion<S,N>(elem_, size_);
    else { // cast elem_ to S type and return new nion
        elem_type new_elem_;
        if constexpr (nion<S,M>::on_heap)  // if M is NION_USE_HEAP, then use heap for memory, else use stack
            new_elem_ = new S[size_];

        for (D i = 0; i < size_; i++)
            new_elem_[i] = static_cast<S>(elem_[i]);

        // if using heap, free new_elem_ to avoid memory leak
        if constexpr (nion<S,M>::on_heap){
            nion<S,M> cast = nion<S,M>(new_elem_, size_);
            if (elem_) delete[] new_elem_; // free new_elem_ to avoid memory leak
            return cast;
        }

        // if using stack, just return new nion
        return nion<S,M>(new_elem_, size_);
    }
}


template<arith_ops T, std::size_t N>
template<std::size_t M, std::size_t P>
constexpr inline nion<T,N> nion<T,N>::make_pair(const nion<T,M> &a, const nion<T,P> &b) {

    /// initialize the nion pair
    nion<T,N> pair;

    /// set the size of the nion
    pair.size_ = a.size_ + b.size_;

    ASSERT(a.size_ > 0 && b.size_ > 0, "The sizes of the nion pair (a, b) must both be greater than zero.");
    ASSERT(pair.size_ <= N || N == 0, "The size of the nion is too large. "
                                      "consider increasing the template parameter, N, or using the default value for heap allocation.\n"
                                      "The value of N is" + std::to_string(N) + ".\n"
                                      "The value of a.size() is " + std::to_string(a.size_) + ".\n"
                                      "The value of b.size() is " + std::to_string(b.size_) + ".");

    if constexpr (on_heap) // if the user wants to use the heap
        pair.elem_ = new T[pair.size_]; // allocate memory on the heap
    // else, the user wants to use the stack

    /// copy the values into the nion
    for (D i = 0; i < a.size_; i++) pair.elem_[i] = a.elem_[i];
    for (D i = 0; i < b.size_; i++) pair.elem_[i + a.size_] = b.elem_[i];


    return pair;
}

template<arith_ops T, std::size_t N>
constexpr inline void nion<T,N>::resize(int size) {
    ASSERT(size > 0, "new nion size must be greater than zero");
    // TODO: make this more efficient by allocating more memory than needed and only resizing when needed.
    if (size_ < size) {
        if constexpr (on_heap) { // if the user wants to use the heap

            elem_type new_elem_ = new T[size]; // allocate memory on the heap (`new` calls the constructor)

            // copy the values into the nion
            for (D i = 0; i < size_; ++i)
                new_elem_[i] = elem_[i];

            if (elem_) delete[] elem_; // free the old memory
            elem_ = new_elem_; // set the new memory
        } else {
            ASSERT(size <= N || N == 0, // else, the user wants to use the stack
                   "The size of the nion is too large. "
                   "consider increasing the template parameter, N, or using the default value for heap allocation."
            );

            // set the new values to zero (default value)
            for (D i = size_; i < size; ++i)
                elem_[i] = T();
        }
    }
    
    size_ = size;
}

template<arith_ops T, std::size_t N>
constexpr nion<T,N> &nion<T,N>::operator=(const std::initializer_list<T> &vals) {
    // set the size of the nion
    resize(vals.size());

    // copy the values into the nion
    for (std::size_t i = 0; i < size_; ++i)
        elem_[i] = *(vals.begin() + i);

    return *this; // return the nion
}


template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> &nion<T,N>::operator=(S scalar) {
    // check if the nion is initialized
    if (size_ <= 0) {
        size_ = 1; // set the degree
        if constexpr (on_heap) { // if the user wants to use the heap
            if (!elem_) elem_ = new T[size_]; // allocate memory on the heap
        }
    }

    zero(); // reset the nion
    elem_[0] = static_cast<T>(scalar); // set the real component
    return *this; // return the nion
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

    // find the smaller size
    uint_com smaller_size;

    if (size_ < other.size_){
        smaller_size = size_;

        // resize this nion to fit the other nion
        resize(other.size_);

    } else smaller_size = other.size_;

    // for the first smaller_size elements, add the other nion's elements to this nion's elements.
    for (uint_com i = 0; i < smaller_size; i++)
        elem_[i] += other.elem_[i];

    // copy the remaining values (if any) from the other nion to this nion.
    for (uint_com i = smaller_size; i < other.size_; i++)
        elem_[i] = other.elem_[i];

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
    nion<T,N> sum(*this);
    sum += other;
    return sum;
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
    return -scalar + z;
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

// create concept to check if T has a `.conj()` method
template<typename T>
concept has_conj = requires(T a) {
    a.conj();
};

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> nion<T,N>::conj() const {
    nion<T,N> conjugate(*this); // copy this nion

    // conjugate the first element if T has a `.conj()` method
    if constexpr (has_conj<T>)
        conjugate.elem_[0] = conjugate.elem_[0].conj();
    // else do nothing

    // negate all components except the first
    for (D i = 1; i < size_; i++)
        conjugate.elem_[i] = -conjugate.elem_[i]; // negate the component

    return conjugate;
}

template<arith_ops T, std::size_t N>
constexpr inline void nion<T,N>::conj_inplace() {
    // conjugate the first element if T has a `.conj()` method
    if constexpr (has_conj<T>)
        elem_[0] = elem_[0].conj();
    // else do nothing

    // negate all components except the first
    for (D i = 1; i < size_; i++)
        elem_[i] = -elem_[i]; // negate the component
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> nion<T,N>::operator-() const {
    nion<T,N> negated(*this); // copy this nion

    // negate all components
    for (D i = 0; i < size_; i++)
            negated.elem_[i] = -negated.elem_[i]; // negate the component

    return negated;
}

template<arith_ops T, std::size_t N>
template<std::size_t M> // set S
constexpr inline nion<T,N> nion<T,N>::operator*(const nion<T,M> &other) const {

    if (size_ == other.size_) {

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds" // Disable the array out-of-bounds warning
#elif defined(__INTEL_COMPILER)
        #pragma warning(push)
            #pragma warning(disable: 2259) // Disable the array out of bounds warning
#endif
        switch (size_) {
            case 1: // if this size is 1, then the product is just the scalar product.

                return {elem_[0] * other.elem_[0]};

// if FULL_RECURSION is defined, then the compiler will use recursion to calculate the product of nions for any size.
#ifndef FULL_RECURSION
            case 2: // hard-coded complex product
                return {
                    {elem_[0] * other.elem_[0] - elem_[1] * other.elem_[1]},
                    {elem_[1] * other.elem_[0] + elem_[0] * other.elem_[1]}
                };


            case 4: // hard-coded quaternion product
                return {
                    {elem_[0] * other.elem_[0] - elem_[1] * other.elem_[1] - elem_[2] * other.elem_[2] - elem_[3] * other.elem_[3]},
                    {elem_[1] * other.elem_[0] + elem_[0] * other.elem_[1] - elem_[3] * other.elem_[2] + elem_[2] * other.elem_[3]},
                    {elem_[2] * other.elem_[0] + elem_[3] * other.elem_[1] + elem_[0] * other.elem_[2] - elem_[1] * other.elem_[3]},
                    {elem_[3] * other.elem_[0] - elem_[2] * other.elem_[1] + elem_[1] * other.elem_[2] + elem_[0] * other.elem_[3]}
                };

            case 8: // hard-coded octonion product ( I know, it's a bit much )
                return {
                    {elem_[0] * other.elem_[0] - elem_[1] * other.elem_[1] - elem_[2] * other.elem_[2] - elem_[3] * other.elem_[3] - elem_[4] * other.elem_[4] - elem_[5] * other.elem_[5] - elem_[6] * other.elem_[6] - elem_[7] * other.elem_[7]},
                    {elem_[1] * other.elem_[0] + elem_[0] * other.elem_[1] - elem_[3] * other.elem_[2] + elem_[2] * other.elem_[3] - elem_[5] * other.elem_[4] + elem_[4] * other.elem_[5] + elem_[7] * other.elem_[6] - elem_[6] * other.elem_[7]},
                    {elem_[2] * other.elem_[0] + elem_[3] * other.elem_[1] + elem_[0] * other.elem_[2] - elem_[1] * other.elem_[3] - elem_[6] * other.elem_[4] - elem_[7] * other.elem_[5] + elem_[4] * other.elem_[6] + elem_[5] * other.elem_[7]},
                    {elem_[3] * other.elem_[0] - elem_[2] * other.elem_[1] + elem_[1] * other.elem_[2] + elem_[0] * other.elem_[3] - elem_[7] * other.elem_[4] + elem_[6] * other.elem_[5] - elem_[5] * other.elem_[6] + elem_[4] * other.elem_[7]},
                    {elem_[4] * other.elem_[0] + elem_[5] * other.elem_[1] + elem_[6] * other.elem_[2] + elem_[7] * other.elem_[3] + elem_[0] * other.elem_[4] - elem_[1] * other.elem_[5] - elem_[2] * other.elem_[6] - elem_[3] * other.elem_[7]},
                    {elem_[5] * other.elem_[0] - elem_[4] * other.elem_[1] + elem_[7] * other.elem_[2] - elem_[6] * other.elem_[3] + elem_[1] * other.elem_[4] + elem_[0] * other.elem_[5] + elem_[3] * other.elem_[6] - elem_[2] * other.elem_[7]},
                    {elem_[6] * other.elem_[0] - elem_[7] * other.elem_[1] - elem_[4] * other.elem_[2] + elem_[5] * other.elem_[3] + elem_[2] * other.elem_[4] - elem_[3] * other.elem_[5] + elem_[0] * other.elem_[6] + elem_[1] * other.elem_[7]},
                    {elem_[7] * other.elem_[0] + elem_[6] * other.elem_[1] - elem_[5] * other.elem_[2] - elem_[4] * other.elem_[3] + elem_[3] * other.elem_[4] + elem_[2] * other.elem_[5] - elem_[1] * other.elem_[6] + elem_[0] * other.elem_[7]}
                };

            case 16: // hard-coded sedenion product ( I strongly recommend using nowrap if you're looking at this... )
                return {
                    {elem_[ 0] * other.elem_[ 0] - elem_[ 1] * other.elem_[ 1] - elem_[ 2] * other.elem_[ 2] - elem_[ 3] * other.elem_[ 3] - elem_[ 4] * other.elem_[ 4] - elem_[ 5] * other.elem_[ 5] - elem_[ 6] * other.elem_[ 6] - elem_[ 7] * other.elem_[ 7] - elem_[ 8] * other.elem_[ 8] - elem_[ 9] * other.elem_[ 9] - elem_[10] * other.elem_[10] - elem_[11] * other.elem_[11] - elem_[12] * other.elem_[12] - elem_[13] * other.elem_[13] - elem_[14] * other.elem_[14] - elem_[15] * other.elem_[15]},
                    {elem_[ 0] * other.elem_[ 1] + elem_[ 1] * other.elem_[ 0] + elem_[ 2] * other.elem_[ 3] - elem_[ 3] * other.elem_[ 2] + elem_[ 4] * other.elem_[ 5] - elem_[ 5] * other.elem_[ 4] - elem_[ 6] * other.elem_[ 7] + elem_[ 7] * other.elem_[ 6] + elem_[ 8] * other.elem_[ 9] - elem_[ 9] * other.elem_[ 8] - elem_[10] * other.elem_[11] + elem_[11] * other.elem_[10] - elem_[12] * other.elem_[13] + elem_[13] * other.elem_[12] + elem_[14] * other.elem_[15] - elem_[15] * other.elem_[14]},
                    {elem_[ 0] * other.elem_[ 2] - elem_[ 1] * other.elem_[ 3] + elem_[ 2] * other.elem_[ 0] + elem_[ 3] * other.elem_[ 1] + elem_[ 4] * other.elem_[ 6] + elem_[ 5] * other.elem_[ 7] - elem_[ 6] * other.elem_[ 4] - elem_[ 7] * other.elem_[ 5] + elem_[ 8] * other.elem_[10] + elem_[ 9] * other.elem_[11] - elem_[10] * other.elem_[ 8] - elem_[11] * other.elem_[ 9] - elem_[12] * other.elem_[14] - elem_[13] * other.elem_[15] + elem_[14] * other.elem_[12] + elem_[15] * other.elem_[13]},
                    {elem_[ 0] * other.elem_[ 3] + elem_[ 1] * other.elem_[ 2] - elem_[ 2] * other.elem_[ 1] + elem_[ 3] * other.elem_[ 0] + elem_[ 4] * other.elem_[ 7] - elem_[ 5] * other.elem_[ 6] + elem_[ 6] * other.elem_[ 5] - elem_[ 7] * other.elem_[ 4] + elem_[ 8] * other.elem_[11] - elem_[ 9] * other.elem_[10] + elem_[10] * other.elem_[ 9] - elem_[11] * other.elem_[ 8] - elem_[12] * other.elem_[15] + elem_[13] * other.elem_[14] - elem_[14] * other.elem_[13] + elem_[15] * other.elem_[12]},
                    {elem_[ 0] * other.elem_[ 4] - elem_[ 1] * other.elem_[ 5] - elem_[ 2] * other.elem_[ 6] - elem_[ 3] * other.elem_[ 7] + elem_[ 4] * other.elem_[ 0] + elem_[ 5] * other.elem_[ 1] + elem_[ 6] * other.elem_[ 2] + elem_[ 7] * other.elem_[ 3] + elem_[ 8] * other.elem_[12] + elem_[ 9] * other.elem_[13] + elem_[10] * other.elem_[14] + elem_[11] * other.elem_[15] - elem_[12] * other.elem_[ 8] - elem_[13] * other.elem_[ 9] - elem_[14] * other.elem_[10] - elem_[15] * other.elem_[11]},
                    {elem_[ 0] * other.elem_[ 5] + elem_[ 1] * other.elem_[ 4] - elem_[ 2] * other.elem_[ 7] + elem_[ 3] * other.elem_[ 6] - elem_[ 4] * other.elem_[ 1] + elem_[ 5] * other.elem_[ 0] - elem_[ 6] * other.elem_[ 3] + elem_[ 7] * other.elem_[ 2] + elem_[ 8] * other.elem_[13] - elem_[ 9] * other.elem_[12] + elem_[10] * other.elem_[15] - elem_[11] * other.elem_[14] + elem_[12] * other.elem_[ 9] - elem_[13] * other.elem_[ 8] + elem_[14] * other.elem_[11] - elem_[15] * other.elem_[10]},
                    {elem_[ 0] * other.elem_[ 6] + elem_[ 1] * other.elem_[ 7] + elem_[ 2] * other.elem_[ 4] - elem_[ 3] * other.elem_[ 5] - elem_[ 4] * other.elem_[ 2] + elem_[ 5] * other.elem_[ 3] + elem_[ 6] * other.elem_[ 0] - elem_[ 7] * other.elem_[ 1] + elem_[ 8] * other.elem_[14] - elem_[ 9] * other.elem_[15] - elem_[10] * other.elem_[12] + elem_[11] * other.elem_[13] + elem_[12] * other.elem_[10] - elem_[13] * other.elem_[11] - elem_[14] * other.elem_[ 8] + elem_[15] * other.elem_[ 9]},
                    {elem_[ 0] * other.elem_[ 7] - elem_[ 1] * other.elem_[ 6] + elem_[ 2] * other.elem_[ 5] + elem_[ 3] * other.elem_[ 4] - elem_[ 4] * other.elem_[ 3] - elem_[ 5] * other.elem_[ 2] + elem_[ 6] * other.elem_[ 1] + elem_[ 7] * other.elem_[ 0] + elem_[ 8] * other.elem_[15] + elem_[ 9] * other.elem_[14] - elem_[10] * other.elem_[13] - elem_[11] * other.elem_[12] + elem_[12] * other.elem_[11] + elem_[13] * other.elem_[10] - elem_[14] * other.elem_[ 9] - elem_[15] * other.elem_[ 8]},
                    {elem_[ 0] * other.elem_[ 8] - elem_[ 1] * other.elem_[ 9] - elem_[ 2] * other.elem_[10] - elem_[ 3] * other.elem_[11] - elem_[ 4] * other.elem_[12] - elem_[ 5] * other.elem_[13] - elem_[ 6] * other.elem_[14] - elem_[ 7] * other.elem_[15] + elem_[ 8] * other.elem_[ 0] + elem_[ 9] * other.elem_[ 1] + elem_[10] * other.elem_[ 2] + elem_[11] * other.elem_[ 3] + elem_[12] * other.elem_[ 4] + elem_[13] * other.elem_[ 5] + elem_[14] * other.elem_[ 6] + elem_[15] * other.elem_[ 7]},
                    {elem_[ 0] * other.elem_[ 9] + elem_[ 1] * other.elem_[ 8] - elem_[ 2] * other.elem_[11] + elem_[ 3] * other.elem_[10] - elem_[ 4] * other.elem_[13] + elem_[ 5] * other.elem_[12] + elem_[ 6] * other.elem_[15] - elem_[ 7] * other.elem_[14] - elem_[ 8] * other.elem_[ 1] + elem_[ 9] * other.elem_[ 0] - elem_[10] * other.elem_[ 3] + elem_[11] * other.elem_[ 2] - elem_[12] * other.elem_[ 5] + elem_[13] * other.elem_[ 4] + elem_[14] * other.elem_[ 7] - elem_[15] * other.elem_[ 6]},
                    {elem_[ 0] * other.elem_[10] + elem_[ 1] * other.elem_[11] + elem_[ 2] * other.elem_[ 8] - elem_[ 3] * other.elem_[ 9] - elem_[ 4] * other.elem_[14] - elem_[ 5] * other.elem_[15] + elem_[ 6] * other.elem_[12] + elem_[ 7] * other.elem_[13] - elem_[ 8] * other.elem_[ 2] + elem_[ 9] * other.elem_[ 3] + elem_[10] * other.elem_[ 0] - elem_[11] * other.elem_[ 1] - elem_[12] * other.elem_[ 6] - elem_[13] * other.elem_[ 7] + elem_[14] * other.elem_[ 4] + elem_[15] * other.elem_[ 5]},
                    {elem_[ 0] * other.elem_[11] - elem_[ 1] * other.elem_[10] + elem_[ 2] * other.elem_[ 9] + elem_[ 3] * other.elem_[ 8] - elem_[ 4] * other.elem_[15] + elem_[ 5] * other.elem_[14] - elem_[ 6] * other.elem_[13] + elem_[ 7] * other.elem_[12] - elem_[ 8] * other.elem_[ 3] - elem_[ 9] * other.elem_[ 2] + elem_[10] * other.elem_[ 1] + elem_[11] * other.elem_[ 0] - elem_[12] * other.elem_[ 7] + elem_[13] * other.elem_[ 6] - elem_[14] * other.elem_[ 5] + elem_[15] * other.elem_[ 4]},
                    {elem_[ 0] * other.elem_[12] + elem_[ 1] * other.elem_[13] + elem_[ 2] * other.elem_[14] + elem_[ 3] * other.elem_[15] + elem_[ 4] * other.elem_[ 8] - elem_[ 5] * other.elem_[ 9] - elem_[ 6] * other.elem_[10] - elem_[ 7] * other.elem_[11] - elem_[ 8] * other.elem_[ 4] + elem_[ 9] * other.elem_[ 5] + elem_[10] * other.elem_[ 6] + elem_[11] * other.elem_[ 7] + elem_[12] * other.elem_[ 0] - elem_[13] * other.elem_[ 1] - elem_[14] * other.elem_[ 2] - elem_[15] * other.elem_[ 3]},
                    {elem_[ 0] * other.elem_[13] - elem_[ 1] * other.elem_[12] + elem_[ 2] * other.elem_[15] - elem_[ 3] * other.elem_[14] + elem_[ 4] * other.elem_[ 9] + elem_[ 5] * other.elem_[ 8] + elem_[ 6] * other.elem_[11] - elem_[ 7] * other.elem_[10] - elem_[ 8] * other.elem_[ 5] - elem_[ 9] * other.elem_[ 4] + elem_[10] * other.elem_[ 7] - elem_[11] * other.elem_[ 6] + elem_[12] * other.elem_[ 1] + elem_[13] * other.elem_[ 0] + elem_[14] * other.elem_[ 3] - elem_[15] * other.elem_[ 2]},
                    {elem_[ 0] * other.elem_[14] - elem_[ 1] * other.elem_[15] - elem_[ 2] * other.elem_[12] + elem_[ 3] * other.elem_[13] + elem_[ 4] * other.elem_[10] - elem_[ 5] * other.elem_[11] + elem_[ 6] * other.elem_[ 8] + elem_[ 7] * other.elem_[ 9] - elem_[ 8] * other.elem_[ 6] - elem_[ 9] * other.elem_[ 7] - elem_[10] * other.elem_[ 4] + elem_[11] * other.elem_[ 5] + elem_[12] * other.elem_[ 2] - elem_[13] * other.elem_[ 3] + elem_[14] * other.elem_[ 0] + elem_[15] * other.elem_[ 1]},
                    {elem_[ 0] * other.elem_[15] + elem_[ 1] * other.elem_[14] - elem_[ 2] * other.elem_[13] - elem_[ 3] * other.elem_[12] + elem_[ 4] * other.elem_[11] + elem_[ 5] * other.elem_[10] - elem_[ 6] * other.elem_[ 9] + elem_[ 7] * other.elem_[ 8] - elem_[ 8] * other.elem_[ 7] + elem_[ 9] * other.elem_[ 6] - elem_[10] * other.elem_[ 5] - elem_[11] * other.elem_[ 4] + elem_[12] * other.elem_[ 3] + elem_[13] * other.elem_[ 2] - elem_[14] * other.elem_[ 1] + elem_[15] * other.elem_[ 0]}
                };
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

    // Helpers for creating nions
    constexpr auto make_pair_heap = [](auto& lhalf, auto& rhalf, T* elem, uint_com size) {
        uint_com half_size = size - (size >> 1);
        lhalf.size_ = half_size; rhalf.size_ = size - half_size;
        lhalf.elem_ = elem; rhalf.elem_ = elem + half_size;
    };

    constexpr auto make_this_pair_stack = [](const T* elem, std::size_t size) {
        uint_com half_size = size - (size >> 1);
        return std::pair<nion<T, N-(N>>1)>, nion<T, N-(N>>1)>>
            {{elem, half_size}, {elem + half_size, size - half_size}};
    };

    constexpr auto make_other_pair_stack = [](const T* elem, std::size_t size) {
        uint_com half_size = size - (size >> 1);
        return std::pair<nion<T, M-(M>>1)>, nion<T, M-(M>>1)>>
                {{elem, half_size}, {elem + half_size, size - half_size}};
    };

    /// compute product depending on where the nions are stored (stack or heap)

    if constexpr (on_heap && other.on_heap) {
        // if both nions are on the heap, only copy the pointers

        nion<T,N> a, b, c, d;
        make_pair_heap(a, b, elem_, size_); // left half
        make_pair_heap(c, d, other.elem_, other.size_); // right half

        // compute the product
        nion<T,N> result = make_pair(
                a*c - d.conj() * b,
                d*a + b*c.conj()
        );

        // set null pointers
        a.elem_ = nullptr; b.elem_ = nullptr; c.elem_ = nullptr; d.elem_ = nullptr;
        return result; // return the product
    } else if constexpr (on_heap) {
        // if only the first nion is on the heap, copy the pointers to the first half of the nion

        nion<T,N> a, b;
        make_pair_heap(a, b, elem_, size_); // left half on heap
        auto [c, d] = make_other_pair_stack(other.elem_, other.size_); // right half on stack

        auto ac = a*c; auto da = d*a; // compute the products
        d.conj_inplace(); c.conj_inplace(); // conjugate the right half in place
        auto db = d*b; auto bc = b*c; // compute the products

        // compute the product
        nion<T,N> result = make_pair(
                ac - db,
                da + bc
        );

        // set null pointers
        a.elem_ = nullptr; b.elem_ = nullptr;
        return result; // return the product
    } else if constexpr (other.on_heap){
        // if only the second nion is on the heap, copy the pointers to the second half of the nion

        auto [a, b] = make_this_pair_stack(elem_, size_); // left half on stack
        nion<T,N> c, d;
        make_pair_heap(c, d, other.elem_, other.size_); // right half on heap

        // compute the product
        nion<T,N> result = make_pair(
                a*c - d.conj() * b,
                d*a + b*c.conj()
        );

        // set null pointers
        c.elem_ = nullptr; d.elem_ = nullptr;
        return result; // return the product
    } else { // if both nions are on the stack, copy the elements to the stack
        auto [a, b] = make_this_pair_stack(elem_, size_); // left half
        auto [c, d] = make_other_pair_stack(other.elem_, other.size_); // right half

        auto ac = a*c; auto da = d*a; // compute the products
        d.conj_inplace(); c.conj_inplace(); // conjugate the right half in place
        auto db = d*b; auto bc = b*c; // compute the products

        // compute the product
        nion<T,N> result = make_pair(
                ac - db,
                da + bc
        );

        return result; // return the product
    }
}

template<arith_ops T, std::size_t N>
template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> nion<T,N>::operator*(S scalar) const {
    nion<T,N> product(*this);

    // compute the product of each element of the nion with the scalar
    for (D i = 0; i < size_; i++)
        product.elem_[i] *= scalar;

    return product;
}

template<arith_ops T, std::size_t N, not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> operator*(S scalar, const nion<T,N> &z) {
    nion<T,N> product(z);

    // compute the product of each element of the nion with the scalar
    using D = typename nion<T,N>::D;
    for (D i = 0; i < z.size_; i++)
        product.elem_[i] = scalar * product.elem_[i];

    return product;
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
    return sqrt(this->abs());
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> nion<T,N>::inv() const {

    nion<T,N> inverse(*this);

    T absolute = this->abs();
    inverse.elem_[0] /= absolute;
    for (D i = 1; i < size_; i++)
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
    for (D i = 0; i < size_; i++)
        quotient.elem_[i] /= scalar;

    return quotient;
}

template<arith_ops T, std::size_t N, not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
constexpr inline nion<T,N> operator/(S scalar, const nion<T,N> &z) {
    nion<T,N> quotient = z.inv();

    // compute the product of each element of the nion with the scalar
    using D = typename nion<T,N>::D;
    for (D i = 0; i < z.size_; i++)
        quotient.elem_[i] = scalar * quotient.elem_[i];

    return quotient;
}


/******************************************
*            ACCESSOR FUNCTIONS
*******************************************/

template<arith_ops T, std::size_t N>
constexpr inline T nion<T,N>::operator[](D index) const {
    if (index >= size_)
        throw std::out_of_range("nion index out of range");
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
    // get polar form of nion
    T r = elem_[0];

    // compute norms
    T z_abs = abs();
    T i_norm = sqrt(z_abs - r * r);

    // compute argument
    T theta = atan2(i_norm, r);
    return theta;
}

template<arith_ops T, std::size_t N>
constexpr inline nion<T,N> nion<T,N>::unit() const {
    // get imaginary part of nion
    nion<T,N> i = imag();

    // compute norm
    T i_abs = i.abs();

    // compute denorm_min
    T denorm_min = T(0); // default value
    if constexpr (std::is_arithmetic_v<T>) {
        denorm_min = std::numeric_limits<T>::denorm_min(); // if T is arithmetic, use its denorm_min
    }

    // compute unit nion
    if (i_abs <= denorm_min) return i; // if i is zero, return i

    // otherwise, compute unit nion
    T i_norm = sqrt(i_abs);
    return i / i_norm;
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
    if (*this == nion<T,N>(scalar, this->size_))
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
    if (*this == nion<T,N>(scalar, this->size_))
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

    // compute denorm_min
    T denorm_min = T(0); // default value
    if constexpr (std::is_arithmetic_v<T>) {
        denorm_min = std::numeric_limits<T>::denorm_min(); // if T is arithmetic, use its denorm_min
    }

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
    // get polar form of nion
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
    else i = nion<T,N>(0, base.size_);

    if constexpr (std::is_integral_v<S>) { // if power is integer, use faster algorithm
        nion<T,N> z = base;
        if (std::signbit(power)) {
            z = inv(z);
            power = -power;
        }

        switch (power) {
            case 0:
                return nion<T,N>(1, z.size_);
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
    return pow(z_norm, power_t) * (cos_ptheta + i * (sin_ptheta));
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

    // compute denorm_min
    T denorm_min = T(0); // default value
    if constexpr (std::is_arithmetic_v<T>) {
        denorm_min = std::numeric_limits<T>::denorm_min(); // if T is arithmetic, use its denorm_min
    }

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

// trigonometric functions
#include "nion_trig.hpp"

#endif //NION_CPP