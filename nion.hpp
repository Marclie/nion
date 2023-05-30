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

#ifndef NION_HPP
#define NION_HPP

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

// macro to determine if heap should be used for memory (zero enables heap)
#define NION_USE_HEAP 0

/**
 * @file nion.hpp
 * @brief Templated class that implements Cayley-Dickson division algebra.
 * @author Marcus Dante Liebenthal
 * @version 1.0
 * @date 2022-10-08
 */
template<arith_ops T, // type of the coefficients (required to have arithmetic operations)
        std::size_t N = NION_USE_HEAP> // type of the size
                                       // (default is max size, where heap is used for memory instead of stack)
struct nion; // forward declaration

// create a concept that checks if a type is the same as this nion type
template<typename S, typename T, std::size_t N>
concept not_nion = !std::is_same_v<T, nion<T, N>>;

template<arith_ops T, std::size_t N>
struct nion {

    // determine minimum number of bits required to represent N at compile time
    static constexpr std::size_t N_BITS = std::numeric_limits<std::size_t>::digits - __builtin_clzl(N);
    static constexpr bool on_heap = N == NION_USE_HEAP; // determine if heap should be used for memory

    // set the internal integer type to the smallest width that can hold N
    using D = std::conditional_t<on_heap, std::size_t, // default to std::size_t if N is 0 for heap allocation
              std::conditional_t<N_BITS <= 8, uint8_t,     // max size of 256
              std::conditional_t<N_BITS <= 16, uint16_t,   // max size of 65536
              std::conditional_t<N_BITS <= 32, uint32_t,   // max size of 4294967296
              std::size_t>>>>; // max size of 18446744073709551616
            
    /// coefficients
    // if N is NION_USE_HEAP, then use heap for memory, else use stack
    using elem_type = std::conditional_t<N == NION_USE_HEAP, T*, T[N]>; // container for coefficients

    elem_type elem_; // declare array of coefficients on stack (where max size is N)
    D size_; // number of coefficients

    /***************************
    *  NION CONSTRUCTORS
    ***************************/

    /**
     * @brief Default constructor.
     * @details Constructs a nion object of size 1 with all coefficients set to 0.
     */
    constexpr inline nion<T,N>() : size_(1) {
        if constexpr (on_heap) elem_ = (T*)malloc(sizeof(T)); // allocate memory on heap
        else elem_[0] = T(); // set first coefficient to default value (0)
    }

    /**
     * @brief Destroy the nion object
     */
    ~nion() {
        if constexpr (on_heap) { // if the nion is on the heap
            if (elem_) free(elem_); // free the memory
            elem_ = nullptr; // set the pointer to null
        }
        // else the nion is on the stack and will be freed automatically
    }

    /**
     * @brief cast nion if both arith_ops types are different and max size is different
     */
    template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S>)
    constexpr inline explicit operator nion<S,M>();

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
        memset(elem_, 0, sizeof(T) * size_); // set all coefficients to zero
    };




    /**
     * @brief copy constructor for different nion types
     * @param other The nion to copy.
     * @return A copy of the nion.
     * @note This is a deep copy.
     */
    constexpr inline nion<T,N>(const nion<T,N> &other);

    /**
     * @brief copy assignment operator for different nion types
     * @param other The nion to copy.
     * @return A copy of the nion.
     * @note This is a deep copy.
     */
    constexpr inline nion<T,N> &operator=(const nion<T,N> &other);

    /**
     * @brief move constructor for different nion types
     * @param other The nion to move.
     * @return A moved nion.
     * @note This is a shallow move. The source nion is left in a valid but unspecified state.
     */
    constexpr inline nion<T,N>(nion<T,N> &&other) noexcept;

    /**
     * @brief move assignment operator for different nion types
     * @param other The nion to move.
     * @return The current nion after moving.
     * @note This is a shallow move. The source nion is left in a valid but unspecified state.
     */
    constexpr inline nion<T,N> &operator=(nion<T,N> &&other) noexcept;






    /**
     * @brief copy constructor for different nion types
     * @param other The nion to copy.
     * @return A copy of the nion.
     * @note This is a deep copy.
     */
    template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S> && (M != N || !std::is_same_v<T, S>))
    constexpr inline explicit nion<T,N>(const nion<S,M> &other);

    /**
     * @brief copy assignment operator for different nion types
     * @param other The nion to copy.
     * @return A copy of the nion.
     * @note This is a deep copy.
     */
    template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S> && (M != N || !std::is_same_v<T, S>))
    constexpr inline nion<T,N> &operator=(const nion<S,M> &other);

    /**
     * @brief move constructor for different nion types
     * @param other The nion to move.
     * @return A moved nion.
     * @note This is a shallow move. The source nion is left in a valid but unspecified state.
     */
    template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S> && (M != N || !std::is_same_v<T, S>))
    constexpr inline explicit nion<T,N>(nion<S,M> &&other) noexcept;

    /**
     * @brief move assignment operator for different nion types
     * @param other The nion to move.
     * @return The current nion after moving.
     * @note This is a shallow move. The source nion is left in a valid but unspecified state.
     */
    template<arith_ops S, std::size_t M> requires (std::is_convertible_v<T, S> && (M != N || !std::is_same_v<T, S>))
    constexpr inline nion<T,N> &operator=(nion<S,M> &&other) noexcept;







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
     * @brief assignment operator from initializer list
     * @param vals The initializer list to copy.
     * @return A nion with the values of the initializer list.
     */
    constexpr inline nion<T,N> &operator=(const std::initializer_list<T> &vals);

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
     * @brief make this nion the conjugate of itself.
     * @return reference to this nion.
     * @detail (a,b)* = (a*,-b)
     */
    constexpr inline void conj_inplace();

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
     * @param lhalf The nion to multiply this nion by.
     * @return The product of this nion and the lhalf nion.
     * @details The product of two nions is defined as (a,b)(c,d) = (a c - d* b, d a + b c*) = (z,w), where * is the conjugate.
     * and a, b, c, d are the nions with half the size of the original nions.
     * @note product has the same size as the larger size of the two nions.
     * @note This is recursive function and will call itself until the size is 1.
     */
    template<std::size_t M> // M is the size of the lhalf nion.
    constexpr inline nion<T,N> operator*(const nion<T,M> &lhalf) const;

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
     * @brief Compute the unit nion.
     * @return The unit nion.
     */
    constexpr inline nion<T,N> unit() const;

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
};

// struct nion

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
    template<arith_ops T, std::size_t N = NION_USE_HEAP, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline nion<T,N> operator*(S scalar, const nion<T,N> &z);

/**
 * @brief overload the / operator for lhs scalars and rhs nions.
 * @tparam T type of the nion.
 * @tparam S type of the scalar.
 * @param scalar type of the scalar.
 * @param z The nion to divide the scalar by.
 * @return
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline nion<T,N> operator/(S scalar, const nion<T,N> &z);

/**
* @brief overload the + operator for lhs scalars and rhs nions.
 * @tparam S type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to add the scalar by.
*/
    template<arith_ops T, std::size_t N = NION_USE_HEAP, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline nion<T,N> operator+(S scalar, const nion<T,N> &z);

/**
 * @brief overload the - operator for lhs scalars and rhs nions.
 * @tparam S The type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to subtract the scalar by.
*/
    template<arith_ops T, std::size_t N = NION_USE_HEAP, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
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
    template<arith_ops T, std::size_t N = NION_USE_HEAP, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
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
    template<arith_ops T, std::size_t N = NION_USE_HEAP, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline bool operator!=(S scalar, const nion<T,N> &z);

/**
 * @brief overload the > operator for lhs scalars and rhs nions.
 * @tparam S type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to compare the scalar to.
 * @return True if the scalar is greater than the nion, false otherwise.
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline bool operator>(S scalar, const nion<T,N> &z);

/**
 * @brief overload the < operator for lhs scalars and rhs nions.
 * @tparam S type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to compare the scalar to.
 * @return True if the scalar is less than the nion, false otherwise.
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline bool operator<(S scalar, const nion<T,N> &z);

/**
 * @brief overload the >= operator for lhs scalars and rhs nions.
 * @tparam S type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to compare the scalar to.
 * @return True if the scalar is greater than or equal to the nion, false otherwise.
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline bool operator>=(S scalar, const nion<T,N> &z);

/**
 * @brief overload the <= operator for lhs scalars and rhs nions.
 * @tparam S type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to compare the scalar to.
 * @return True if the scalar is less than or equal to the nion, false otherwise.
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP, not_nion<T,N> S = T> requires (std::is_convertible_v<S,T>)
    extern constexpr inline bool operator<=(S scalar, const nion<T,N> &z);

/**
 * @brief overload the << operator for nions.
 * @param os The output stream.
 * @param z The nion to print.
 * @return The output stream.
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP>
    extern std::ostream &operator<<(std::ostream &os, const nion<T,N> &z);

/**
 * @brief overload the >> operator for nions.
 * @param is The input stream.
 * @param z The nion to read into.
 * @return The input stream.
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP>
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
    template<arith_ops T, std::size_t N = NION_USE_HEAP>
    constexpr inline T real(const nion<T,N> &z);

/**
 * @brief Calculate the imaginary part of a nion.
 * @tparam T type of the nion.
 * @param z The nion to calculate the imaginary part of.
 * @return The imaginary part of the nion.
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP>
    constexpr inline nion<T,N> imag(const nion<T,N> &z);

/**
 * @brief compute the conjugate of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the conjugate of.
 * @return The conjugate of the nion.
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP>
    constexpr inline nion<T,N> conj(const nion<T,N> &z);

/**
 * @brief compute the absolute value of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the absolute value of.
 * @return The absolute value of the nion.
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP>
    constexpr inline T abs(const nion<T,N> &z);

/**
 * @brief compute the norm of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the norm of.
 * @return The norm of the nion.
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP>
    constexpr inline T norm(const nion<T,N> &z);

/**
 * @brief compute the inverse of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse of.
 * @return The inverse of the nion.
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP>
    constexpr inline nion<T,N> inv(const nion<T,N> &z);

/**
 * @brief compute the dot product of two nions.
 * @tparam T type of the nions.
 * @param lhs The left hand side nion.
 * @param rhs The right hand side nion.
 * @return The dot product of the nions.
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP, std::size_t M = N>
    constexpr inline T dot(const nion<T,N> &lhs, const nion<T,M> &rhs);


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
    template<arith_ops T, std::size_t N = NION_USE_HEAP>
    constexpr inline nion<T,N> exp(const nion<T,N> &z);

/**
 * @brief compute the principle logarithm of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the natural logarithm of.
 * @return The natural logarithm of the nion.
 * @details The natural logarithm of a nion is defined as ln(z) = ln(|z|) + v/|v| * atan(|v|/r).
 *          where a is the real component and v is the imaginary components.
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP>
    constexpr inline nion<T,N> log(const nion<T,N> &z);

/**
 * @brief compute the power of a nion with scalar.
 * @tparam T type of the nion.
 * @tparam S type of the power.
 * @param base The nion to compute the power of.
 * @param power The power to raise the nion to.
 * @return The power of the nion.
 * @details The power of a nion is defined as z^p = e^(p * ln(z)).
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP, not_nion<T,N> S> requires (std::is_convertible_v<S,T>)
    constexpr inline nion<T,N> pow(const nion<T,N> &base, S power);

/**
 * @brief compute the power of a nion with nion.
 * @tparam T type of the nion.
 * @param base The nion to compute the power of.
 * @param power The power to raise the nion to.
 * @return The power of the nion.
 * @details The power of a nion is defined as z^p = e^(p * ln(z)).
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP, std::size_t M = N>
    constexpr inline nion<T,N> pow(const nion<T,N> &base, const nion<T,M> &power);

/**
 * @brief compute the power of a scalar with a nion.
 * @tparam T type of the nion.
 * @tparam S type of the power.
 * @param power The nion to use as a power.
 * @param base The scalar to raise by the nion
 * @return The power of the nion.
 * @details The power of with a nion is defined as x^z = e^(z * ln(x)).
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP, not_nion<T,N> S = T>
    constexpr inline nion<T,N> pow(S base, const nion<T,N> &power);

/**
 * @brief compute the square of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the square  of.
 * @return The square of the nion.
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP>
    constexpr inline nion<T,N> sqr(const nion<T,N> &base);

/**
 * @brief compute the square root of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the square root of.
 * @return The square root of the nion.
 * @details The square root of a nion is defined as sqrt(z) = z^(1/2).
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP>
    constexpr inline nion<T,N> sqrt(const nion<T,N> &z);

/**
 * @brief compute the cube root of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the cube root of.
 * @return The cube root of the nion.
 * @details The cube root of a nion is defined as cbrt(z) = z^(1/3).
 */
    template<arith_ops T, std::size_t N = NION_USE_HEAP>
    constexpr inline nion<T,N> cbrt(const nion<T,N> &z);

#include "nion.cpp"

#endif //NION_HPP
