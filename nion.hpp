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

#ifndef NION_NION_HPP
#define NION_NION_HPP

#include <cmath>
#include <limits>
#include <type_traits>
#include <iostream>
#include <stdexcept>
#include <cstring>
#include <concepts>

// Make a concept that checks if a type has arithmetic operations
template<typename T>
concept arith_ops = requires(T a, T b) {
    { a + b } -> std::same_as<T>;
    { a - b } -> std::same_as<T>;
    { a * b } -> std::same_as<T>;
    { a / b } -> std::same_as<T>;
};

// make a concept that checks if a type is integral
template<typename Z>
concept integral_type = std::is_integral_v<Z>;

/**
 * @file nion.hpp
 * @brief Templated class that implements Cayley-Dickson division algebra.
 * @author Marcus Dante Liebenthal
 * @version 1.0
 * @date 2022-10-08
 */
template<arith_ops T, // type of the coefficients
        std::size_t N = 32> // type of the size
struct nion; // forward declaration

// create a concept that checks if a type is the same as this nion type
template<typename S, typename T, std::size_t N>
concept not_nion = !std::is_same_v<T, nion<T, N>>;

template<arith_ops T, std::size_t N>
struct nion {

    // determine minimum number of bits required to represent N at compile time
    static constexpr std::size_t N_BITS = std::numeric_limits<std::size_t>::digits - __builtin_clzl(N);

    // set the internal integer type to the smallest width that can hold N
    using D = std::conditional_t<N_BITS <= 8, uint8_t, // max size of 256
            std::conditional_t<N_BITS <= 16, uint16_t, // max size of 65536
            std::conditional_t<N_BITS <= 32, uint32_t, // max size of 4294967296
            std::conditional_t<N_BITS <= 64, uint64_t, // max size of 18446744073709551616
            std::conditional_t<N_BITS <= 128, __uint128_t, void>>>>>; // max size of 340282366920938463463374607431768211456


    /// coefficients
    T elem_[N]; // declare array of coefficients on stack (where max size is N)
    D size_;

    /***************************
    *  NION CONSTRUCTORS
    ***************************/

    /**
     * @brief Default constructor.
     * @details Constructs a null nion object.
     */
    constexpr inline nion<T,N>() : size_(1) { zero();}

    /**
     * @brief cast nion if both arith_ops types are different and max size is different
     */
    template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S>)
    constexpr inline explicit operator nion<S,M>(){
        if constexpr (std::is_same_v<T, S>) // if types are the same, just copy the array
            return nion<T,M>(elem_, size_);
        else { // cast elem_ to S type and return new nion
            S new_elem_[M];
            for (D i = 0; i < size_; i++)
                new_elem_[i] = static_cast<S>(elem_[i]);
            return nion<S, N>(new_elem_, size_);
        }
    }

    /**
     * @brief Construct a new nion object from vector
     * @param components The components of the nion.
     * @param size The size of the nion.
     * @details Constructs a nion object from a vector of components.
     * @note The size of the nion is the number of components.
     * @note The size of the nion must be greater than zero.
     */
    template <arith_ops S> requires (std::is_convertible_v<S, T>)
    constexpr inline nion<T,N>(const S *vals, size_t size);

    /**
     * @brief Construct a new nion object from vector
     * @param components The components of the nion as an initializer list.
     * @note The size of the nion is determined by the size of the initializer list.
     */
    constexpr inline nion<T,N>(const std::initializer_list<T> &vals);

    /**
     * @brief Construct an empty nion object
     * @param size The size of the nion.
     */
    constexpr inline explicit nion<T,N>(D size);

    /**
     * @brief zero nion
     */
    constexpr inline void zero(){
        memset(elem_, 0, sizeof(T) * size_);
    };

    /**
     * @brief copy constructor
     * @param other The nion to copy.
     * @return A copy of the nion.
     * @note This is a deep copy.
     */
    template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S>)
    constexpr inline explicit nion<T,N>(const nion<S,M> &other);

    /**
     * @brief move constructor
     * @param other The nion to copy.
     * @return A copy of the nion.
     * @note This is a deep copy.
     */
    template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S>)
    constexpr inline explicit nion<T,N>(nion<S,M> &&other) noexcept;

    /**
     * @brief Construct a new nion object from a scalar with no imaginary components.
     * @param size The size of the nion.
     * @param scalar The scalar value of the nion.
     * @return nion<T,N> The nion object.
     * @note This is a convenience function for creating a nion from a scalar.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T> && !std::is_pointer_v<S>)
    constexpr inline nion<T,N>(S realVal, D size);

    /**
     * @brief Destroy the nion object
     */
    ~nion() = default;

    /**
     * @brief Construct nion from half size nions. q = (a,b)
     * @param a The first half size nion in pair.
     * @param b The second half size nion in pair.
     * @return The nion constructed from the half size nions.
     * @note This is a convenience function for constructing a nion from pairing two half size nions.
     */
    template<std::size_t M, std::size_t P> // M is the size of the first nion, P is the size of the second nion.
    static constexpr inline nion<T,N> make_pair(const nion<T,M> &a, const nion<T,P> &b);

    /**
     * @brief resizes the nion to the given size.
     * @param size The new size of the nion.
     */
    constexpr inline void resize(int size);

    /**
     * @brief copy assignment operator
     * @param other The nion to copy.
     * @return A copy of the nion.
     * @note This is a deep copy.
     */
    template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S>)
    constexpr inline nion<T,N> &operator=(const nion<S,M> &other);

    /**
     * @brief assignment operator from initializer list
     * @param vals The initializer list to copy.
     * @return A nion with the values of the initializer list.
     */
    constexpr inline nion<T,N> &operator=(const std::initializer_list<T> &vals);

    /**
     * @brief override move assignment operator
     * @param other The nion to move.
     * @return A copy of the nion.
     * @note This is a shallow copy.
     */
    template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S>)
    constexpr inline nion<T,N> &operator=(nion<S,M> &&other) noexcept;

    /**
     * @brief convert scalar to nion
     * @param scalar The scalar to convert to a nion.
     * @return The nion constructed from the scalar.
     * @note This is a convenience function for constructing a nion from a scalar.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline nion<T,N> &operator=(S scalar);

    /**
     * @brief overload the - operator for a nion.
     * @return The negation of the nion.
     */
    constexpr inline nion<T,N> operator-() const;

    /**
     * @brief overload the [] operator for a nion.
     * @param index The index of the component to get.
     * @return The component at the index passed by value.
     */
    constexpr inline T operator[](D index) const;

    /**
     * @brief Get the conjugate of the nion.
     * @return The conjugate of the nion.
     * @detail (a,b)* = (a*,-b)
     */
    constexpr inline nion<T,N> conj() const;

    /**
     * @brief overload the += operator for nions.
     * @param other The nion to add to this nion.
     * @return The sum of this nion and the other nion inplace.
     */
    template<std::size_t M> // M is the size of the other nion.
    constexpr inline void operator+=(const nion<T,M> &other);

    /**
     * @breif overload the -= operator for nions.
     * @param other The nion to substract from this nion.
     * @return The subtraction of this nion and the other nion inplace.
     */
    template<std::size_t M> // M is the size of the other nion.
    constexpr inline void operator-=(const nion<T,M> &other);

    /**
     * @breif overload the *= operator for nions.
     * @param other The nion to multiply this nion by.
     * @return The product of this nion and the other nion inplace.
     */
    template<std::size_t M> // M is the size of the other nion.
    constexpr inline void operator*=(const nion<T,M> &other);

    /**
     * @breif overload the /= operator for nions.
     * @param other The nion to divide this nion by.
     * @return The division of this nion and the other nion inplace.
     */
    template<std::size_t M> // M is the size of the other nion.
    constexpr inline void operator/=(const nion<T,M> &other);

    /**
     * @brief overload the + operator for nions.
     * @param other The nion to add to this nion.
     * @return The sum of this nion and the other nion.
     */
    template<std::size_t M> // M is the size of the other nion.
    constexpr inline nion<T,N> operator+(const nion<T,M> &other) const;

    /**
     * @brief overload the - operator for nions.
     * @param other The nion to substract this nion by.
     * @return The subtraction of this nion and the other nion.
     */
    template<std::size_t M> // M is the size of the other nion.
    constexpr inline nion<T,N> operator-(const nion<T,M> &other) const;


    /**
     * @brief overload the * operator for nions.
     * @param other The nion to multiply this nion by.
     * @return The product of this nion and the other nion.
     * @details The product of two nions is defined as (a,b)(c,d) = (a c - d* b, d a + b c*) = (z,w), where * is the conjugate.
     * and a, b, c, d are the nions with half the size of the original nions.
     * @note product has the same size as the larger size of the two nions.
     * @note This is recursive function and will call itself until the size is 1.
     */
    template<std::size_t M> // M is the size of the other nion.
    constexpr inline nion<T,N> operator*(const nion<T,M> &other) const;

    /**
     * @brief compute the inverse of the nion.
     * @return The inverse of the nion.
     */
    constexpr inline nion<T,N> inv() const;

    /**
     * @brief overload the / operator for nions.
     * @param other The nion to divide this nion by.
     * @return The division of this nion and the other nion.
     */
    template<std::size_t M> // M is the size of the other nion.
    constexpr inline nion<T,N> operator/(const nion<T,M> &other) const;

    /**
     * @brief absolute value of the nion.
     * @return The absolute value of the nion.
     * @details The absolute value of the nion is the sum of the squares of its components.
     */
    constexpr inline T abs() const;

    /**
     * @brief norm of the nion.
     * @return The norm of the nion.
     * @details The norm of a nion is the sqrt of the sum of the squares of its components.
     */
    constexpr inline T norm() const;

    /**
     * @brief shortest rotation of the nion to the real axis.
     * @tparam T type of the nion.
     * @param nion The nion to rotate
     * @return The shortest rotation of the nion to the real axis (preserves the norm and sign of the real component).
     * @details computes the length of the nion and does the shortest rotation to the real line.
     */
    constexpr inline T rotate_real() const;

    /**
     * @brief return real part of the nion.
     * @return The real part of the nion.
     */
    constexpr inline T real() const;

    /**
     * @brief Calculate the imaginary components of a nion.
     * @return The imaginary components of the nion.
     */
    constexpr inline nion<T,N> imag() const;

    /**
     * @brief Compute the argument of the nion.
     * @return The argument of the nion.
     * @details The argument of a nion is the angle between the real axis and the nion.
     */
    constexpr inline T arg() const;

    /**
     * @brief overload the == operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if the nions are equal, false otherwise.
     * @details Two nions are equal if they have the same size and the same components.
     */
    template<std::size_t M> // M is the size of the other nion.
    constexpr inline bool operator==(const nion<T,M> &other) const;

    /**
     * @brief overload the != operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if the nions are not equal, false otherwise.
     * @details Two nions are equal if they have the same size and the same components.
     */
    template<std::size_t M> // M is the size of the other nion.
    constexpr inline bool operator!=(const nion<T,M> &other) const;

    /**
     * @brief overload the > operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if this nion is greater than the other nion, false otherwise.
     * @details sorting is undefined for nions with sizes greater than 1. However, we can still compare
     *         nions with sizes greater than 1 by comparing the projections of the nions onto the real line.
     */
    template<std::size_t M> // M is the size of the other nion.
    constexpr inline bool operator>(const nion<T,M> &other) const;

    /**
     * @brief overload the < operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if this nion is less than the other nion, false otherwise.
     * @details sorting is undefined for nions with sizes greater than 1. However, we can still compare
     *         nions with sizes greater than 1 by comparing the rotations of the nions onto the real line.
     */
    template<std::size_t M> // M is the size of the other nion.
    constexpr inline bool operator<(const nion<T,M> &other) const;

    /**
     * @brief overload the >= operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if this nion is greater than or equal to the other nion, false otherwise.
     * @details sorting is undefined for nions with sizes greater than 1. However, we can still compare
     *         nions with sizes greater than 1 by comparing the rotations of the nions onto the real line.
     */
    template<std::size_t M> // M is the size of the other nion.
    constexpr inline bool operator>=(const nion<T,M> &other) const;

    /**
     * @brief overload the <= operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if this nion is less than or equal to the other nion, false otherwise.
     * @details sorting is undefined for nions with sizes greater than 1. However, we can still compare
     *         nions with sizes greater than 1 by comparing the rotations of the nions onto the real line.
     */
    template<std::size_t M> // M is the size of the other nion.
    constexpr inline bool operator<=(const nion<T,M> &other) const;

    /**
     * @brief overload the + operator for nions with scalars.
     * @tparam S The type of the scalar.
     * @param scalar The scalar to add this nion by.
     * @return The sum of this nion and the scalar.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline nion<T,N> operator+(S scalar) const;

    /**
     * @brief overload the ++ operator for nions.
     * @return The nion incremented by 1.
     */
    constexpr inline nion<T,N> &operator++();

    /**
     * @brief overload the - operator for nions with scalars.
     * @tparam S The type of the scalar.
     * @param scalar The scalar to subtract this nion by.
     * @return The subtraction of this nion and the scalar.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline nion<T,N> operator-(S scalar) const;

    /**
     * @brief overload the -- operator for nions.
     * @return The real component of the nion decremented by 1.
     */
    constexpr inline nion<T,N> &operator--();

    /**
     * @brief overload the * operator for nions with scalars.
     * @tparam S The type of the scalar.
     * @param scalar The scalar to multiply this nion by.
     * @return The product of this nion and the scalar.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline nion<T,N> operator*(S scalar) const;

    /**
     * @brief overload the / operator for nions with scalars.
     * @tparam S The type of the scalar.
     * @param scalar The scalar to divide this nion by.
     * @return The division of this nion and the scalar.
     */
    template<not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
    constexpr inline nion<T,N> operator/(S scalar) const;

    /**
     * @brief overload the += operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param other The scalar to add to this nion.
     * @return The sum of this nion and the scalar inplace.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline void operator+=(S scalar);

    /**
     * @brief overload the -= operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param other The scalar to substract from this nion.
     * @return The subtraction of this nion and the scalar inplace.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline void operator-=(S scalar);

    /**
     * @breif overload the *= operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param scalar The scalar to multiply this nion by.
     * @return The product of this nion and the scalar inplace.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline void operator*=(S scalar);

    /**
     * @breif overload the /= operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param scalar The scalar to divide this nion by.
     * @return The division of this nion and the scalar inplace.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline void operator/=(S scalar);

    /**
     * @brief overload the == operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param scalar The scalar to compare this nion to.
     * @return True if the nion is equal to the scalar, false otherwise.
     * @details A nion is equal to a scalar if the scalar is equal to the first component of the nion
     *          and all other components are equal to zero.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline bool operator==(S scalar) const;

    /**
     * @brief overload the != operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param scalar The scalar to compare this nion to.
     * @return False if the nion is equal to the scalar, true otherwise.
     * @details A nion is equal to a scalar if the scalar is equal to the first component of the nion
     *          and all other components are equal to zero.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline bool operator!=(S scalar) const;

    /**
     * @brief overload the > operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param scalar The scalar to compare this nion to.
     * @return True if the nion is greater than the scalar, false otherwise.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline bool operator>(S scalar) const;

    /**
     * @brief overload the < operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param scalar The scalar to compare this nion to.
     * @return True if the nion is less than the scalar, false otherwise.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline bool operator<(S scalar) const;

    /**
     * @brief overload the >= operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param scalar The scalar to compare this nion to.
     * @return True if the nion is greater than or equal to the scalar, false otherwise.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline bool operator>=(S scalar) const;

    /**
     * @brief overload the <= operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param scalar The scalar to compare this nion to.
     * @return True if the nion is less than or equal to the scalar, false otherwise.
     */
    template<not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    constexpr inline bool operator<=(S scalar) const;


    /**
     * @brief comparison if nion is a real number.
     * @return True if the nion is a real number, false otherwise.
     * @details A nion is a real number if all components except the first are equal to zero.
     */
    constexpr inline bool is_real() const;

    /**
    * @brief Converts a nion to a string.
    * @return The string representation of the nion.
    */
    inline std::string to_string() const;
}; // struct nion

/***************************
    *  NION OPERATOR OVERLOADS
    ***************************/

/**
 * @brief overload the * operator for lhs scalars and rhs nions.
 * @tparam T type of the nion.
 * @tparam S type of the scalar.
 * @param scalar type of the scalar.
 * @param z The nion to multiply the scalar by.
 * @return
 */
    template<arith_ops T, std::size_t N = 128, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline nion<T,N> operator*(S scalar, const nion<T,N> &z);

/**
 * @brief overload the / operator for lhs scalars and rhs nions.
 * @tparam T type of the nion.
 * @tparam S type of the scalar.
 * @param scalar type of the scalar.
 * @param z The nion to divide the scalar by.
 * @return
 */
    template<arith_ops T, std::size_t N = 128, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline nion<T,N> operator/(S scalar, const nion<T,N> &z);

/**
* @brief overload the + operator for lhs scalars and rhs nions.
 * @tparam S type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to add the scalar by.
*/
    template<arith_ops T, std::size_t N = 128, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline nion<T,N> operator+(S scalar, const nion<T,N> &z);

/**
 * @brief overload the - operator for lhs scalars and rhs nions.
 * @tparam S The type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to subtract the scalar by.
*/
    template<arith_ops T, std::size_t N = 128, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline nion<T,N> operator-(S scalar, const nion<T,N> &z);

/**
 * @brief overload the == operator for lhs scalars and rhs nions.
 * @tparam S type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to compare the scalar to.
 * @return True if the nion is equal to the scalar, false otherwise.
 * @details A nion is equal to a scalar if the scalar is equal to the first component of the nion
 *         and all other components are equal to zero.
 */
    template<arith_ops T, std::size_t N = 128, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline bool operator==(S scalar, const nion<T,N> &z);

/**
 * @brief overload the != operator for lhs scalars and rhs nions.
 * @tparam S type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to compare the scalar to.
 * @return False if the nion is equal to the scalar, true otherwise.
 * @details A nion is equal to a scalar if the scalar is equal to the first component of the nion
 *         and all other components are equal to zero.
 */
    template<arith_ops T, std::size_t N = 128, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline bool operator!=(S scalar, const nion<T,N> &z);

/**
 * @brief overload the > operator for lhs scalars and rhs nions.
 * @tparam S type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to compare the scalar to.
 * @return True if the scalar is greater than the nion, false otherwise.
 */
    template<arith_ops T, std::size_t N = 128, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline bool operator>(S scalar, const nion<T,N> &z);

/**
 * @brief overload the < operator for lhs scalars and rhs nions.
 * @tparam S type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to compare the scalar to.
 * @return True if the scalar is less than the nion, false otherwise.
 */
    template<arith_ops T, std::size_t N = 128, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline bool operator<(S scalar, const nion<T,N> &z);

/**
 * @brief overload the >= operator for lhs scalars and rhs nions.
 * @tparam S type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to compare the scalar to.
 * @return True if the scalar is greater than or equal to the nion, false otherwise.
 */
    template<arith_ops T, std::size_t N = 128, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline bool operator>=(S scalar, const nion<T,N> &z);

/**
 * @brief overload the <= operator for lhs scalars and rhs nions.
 * @tparam S type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to compare the scalar to.
 * @return True if the scalar is less than or equal to the nion, false otherwise.
 */
    template<arith_ops T, std::size_t N = 128, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline bool operator<=(S scalar, const nion<T,N> &z);

/**
 * @brief overload the << operator for nions.
 * @param os The output stream.
 * @param z The nion to print.
 * @return The output stream.
 */
    template<arith_ops T, std::size_t N = 128>
    extern std::ostream &operator<<(std::ostream &os, const nion<T,N> &z);

/**
 * @brief overload the >> operator for nions.
 * @param is The input stream.
 * @param z The nion to read into.
 * @return The input stream.
 */
    template<arith_ops T, std::size_t N = 128>
    extern std::istream &operator>>(std::istream &is, nion<T,N> &z);

/***************************
    *  NION FUNCTION IMPLEMENTATIONS *
    ***************************/

/**
 * @brief Calculate the real part of a nion.
 * @tparam T type of the nion.
 * @param z The nion to calculate the real part of.
 * @return The real part of the nion.
 */
    template<arith_ops T, std::size_t N = 128>
    extern constexpr inline T real(const nion<T,N> &z);

/**
 * @brief Calculate the imaginary part of a nion.
 * @tparam T type of the nion.
 * @param z The nion to calculate the imaginary part of.
 * @return The imaginary part of the nion.
 */
    template<arith_ops T, std::size_t N = 128>
    extern constexpr inline nion<T,N> imag(const nion<T,N> &z);

/**
 * @brief compute the conjugate of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the conjugate of.
 * @return The conjugate of the nion.
 */
    template<arith_ops T, std::size_t N = 128>
    extern constexpr inline nion<T,N> conj(const nion<T,N> &z);

/**
 * @brief compute the absolute value of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the absolute value of.
 * @return The absolute value of the nion.
 */
    template<arith_ops T, std::size_t N = 128>
    extern constexpr inline T abs(const nion<T,N> &z);

/**
 * @brief compute the norm of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the norm of.
 * @return The norm of the nion.
 */
    template<arith_ops T, std::size_t N = 128>
    extern constexpr inline T norm(const nion<T,N> &z);

/**
 * @brief compute the inverse of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse of.
 * @return The inverse of the nion.
 */
    template<arith_ops T, std::size_t N = 128>
    extern constexpr inline nion<T,N> inv(const nion<T,N> &z);

/**
 * @brief compute the dot product of two nions.
 * @tparam T type of the nions.
 * @param lhs The left hand side nion.
 * @param rhs The right hand side nion.
 * @return The dot product of the nions.
 */
    template<arith_ops T, std::size_t N = 128, std::size_t M = N>
    extern constexpr inline T dot(const nion<T,N> &lhs, const nion<T,M> &rhs);


/***************************
    *  NION ALGEBRAIC FUNCTIONS *
    ***************************/

/**
 * @brief compute exponential of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the exponential of.
 * @return The exponential of the nion.
 * @details The exponential of a nion is defined as e^z = e^r * (cos|v| + v/|v| * sin|v|).
 *          where a is the real component and v is the imaginary components.
 */
    template<arith_ops T, std::size_t N = 128>
    extern constexpr inline nion<T,N> exp(const nion<T,N> &z);

/**
 * @brief compute the principle logarithm of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the natural logarithm of.
 * @return The natural logarithm of the nion.
 * @details The natural logarithm of a nion is defined as ln(z) = ln(|z|) + v/|v| * atan(|v|/r).
 *          where a is the real component and v is the imaginary components.
 */
    template<arith_ops T, std::size_t N = 128>
    extern constexpr inline nion<T,N> log(const nion<T,N> &z);

/**
 * @brief compute the power of a nion with scalar.
 * @tparam T type of the nion.
 * @tparam S type of the power.
 * @param base The nion to compute the power of.
 * @param power The power to raise the nion to.
 * @return The power of the nion.
 * @details The power of a nion is defined as z^p = e^(p * ln(z)).
 */
    template<arith_ops T, std::size_t N = 128, not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
    extern constexpr inline nion<T,N> pow(const nion<T,N> &base, S power);

/**
 * @brief compute the power of a nion with nion.
 * @tparam T type of the nion.
 * @param base The nion to compute the power of.
 * @param power The power to raise the nion to.
 * @return The power of the nion.
 * @details The power of a nion is defined as z^p = e^(p * ln(z)).
 */
    template<arith_ops T, std::size_t N = 128, std::size_t M = N>
    extern constexpr inline nion<T,N> pow(const nion<T,N> &base, const nion<T,M> &power);

/**
 * @brief compute the power of a scalar with a nion.
 * @tparam T type of the nion.
 * @tparam S type of the power.
 * @param power The nion to use as a power.
 * @param base The scalar to raise by the nion
 * @return The power of the nion.
 * @details The power of with a nion is defined as x^z = e^(z * ln(x)).
 */
    template<arith_ops T, std::size_t N = 128, not_nion<T,N> S = T>
    extern constexpr inline nion<T,N> pow(S base, const nion<T,N> &power);

/**
 * @brief compute the square of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the square  of.
 * @return The square of the nion.
 */
    template<arith_ops T, std::size_t N = 128>
    extern constexpr inline nion<T,N> sqr(const nion<T,N> &base);

/**
 * @brief compute the square root of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the square root of.
 * @return The square root of the nion.
 * @details The square root of a nion is defined as sqrt(z) = z^(1/2).
 */
    template<arith_ops T, std::size_t N = 128>
    extern constexpr inline nion<T,N> sqrt(const nion<T,N> &z);

/**
 * @brief compute the cube root of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the cube root of.
 * @return The cube root of the nion.
 * @details The cube root of a nion is defined as cbrt(z) = z^(1/3).
 */
    template<arith_ops T, std::size_t N = 128>
    extern constexpr inline nion<T,N> cbrt(const nion<T,N> &z);

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



#include "nion.cpp"

#endif //NION_NION_HPP
