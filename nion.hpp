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
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>
#include <type_traits>
#include <complex>
#include <cstring>

/**
 * @file nion.hpp
 * @brief Templated class that implements Cayley-Dickson division algebra.
 * @author Marcus Dante Liebenthal
 * @version 1.0
 * @date 2022-10-08
 */

// concept such that nion templates are only instantiated for arithmetic types
template <typename T>
concept arithmetic = std::is_arithmetic<T>::value;
template<arithmetic T>
struct nion {
    T* components;
    uint_fast16_t order;
    static constexpr T epsilon = std::numeric_limits<T>::epsilon();

    /**
     * @brief Default constructor.
     * @details Constructs a nion object with only real component with value 0.
     */
    constexpr nion<T>() : components((T*)malloc(sizeof(T))), order(1) {
        components[0] = 0;
    }

    /**
     * @brief Construct a new nion object from vector
     * @param components The components of the nion.
     * @param order The order of the nion.
     */
    constexpr explicit nion<T>(const T* components, uint_fast16_t order) : order(order) {

        // check if the order is greater than zero
        if (order <= 0) {
            throw std::invalid_argument("The order of the nions must be greater than zero.");
        }

        this->components = (T*) malloc(order * sizeof(T));
        std::memcpy(this->components, components, order * sizeof(T));
    };

    /**
     * @brief Construct a new nion object from brace initializer list
     * @param order The order of the nion.
     */
    constexpr nion<T>(const std::initializer_list<T> &components) : order(components.size()) {
        if (order <= 0) {
            throw std::invalid_argument("The order of the nions must be greater than zero.");
        }
        this->components = (T*) malloc(order * sizeof(T));
        std::memcpy(this->components, components.begin(), order * sizeof(T));
    };

    /**
     * @brief Construct an empty nion object
     * @param order The order of the nion.
     */
    constexpr explicit nion<T>(uint_fast16_t order) : order(order) {
        // check if the order is greater than zero
        if (order <= 0) {
            throw std::invalid_argument("The order of the nions must be greater than zero.");
        }
        this->components = (T*) malloc(this->order * sizeof(T));
        std::memset(this->components, 0, order * sizeof(T));
    };

    /**
     * @brief copy constructor
     * @param other The nion to copy.
     * @return A copy of the nion.
     * @note This is a deep copy.
     */
    constexpr nion<T>(const nion<T> &other) : order(other.order) {
        this->components = (T*) malloc(other.order * sizeof(T));
        std::memcpy(this->components, other.components, other.order * sizeof(T));
    };

    /**
     * @brief Construct a new nion object from a scalar with no imaginary components.
     * @param order The order of the nion.
     * @param scalar The scalar value of the nion.
     * @return nion<T> The nion object.
     * @note This is a convenience function for creating a nion from a scalar.
     */
    template<arithmetic S>
    constexpr explicit nion<T>(S realVal, uint_fast16_t order) : order(order) {
        // check if the order is greater than zero
        if (order <= 0) {
            throw std::invalid_argument("The order of the nions must be greater than zero.");
        }
        this->components = (T*) malloc(order * sizeof(T));
        std::memset(this->components, 0, order * sizeof(T));
        components[0] = realVal;
    };

    /**
     * @brief Construct a new nion object from std::complex
     * @param order The order of the nion.
     * @param complex The std::complex object.
     * @return nion<T> The nion object.
     */
    constexpr explicit nion<T>(std::complex<T> complex) : order(2) {
        this->components = (T*) malloc(2 * sizeof(T));
        memset(this->components, 0, 2 * sizeof(T));

        components[0] = complex.real();
        components[1] = complex.imag();
    };

    /**
     * @brief Destroy the nion object
     */
    ~nion<T>() { free(components); components = nullptr; };

    /**
     * @brief get the left pair of the nion. a in (a,b)
     * @tparam T type of the nion.
     * @param nion The nion to get the left pair of.
     * @return The left pair of the nion.
     */
    constexpr nion<T> left() const {
        uint_fast16_t halfOrder = order / 2;
        nion<T> a(order - halfOrder);
        std::memcpy(a.components, components, (order - halfOrder) * sizeof(T));
        return a;
    };

    /**
     * @brief get the right pair of the nion. b in (a,b)
     * @tparam T type of the nion.
     * @param nion The nion to get the right pair of.
     * @return The right pair of the nion.
     */
    constexpr nion<T> right() const {
        uint_fast16_t halfOrder = order / 2;
        nion<T> b(halfOrder);
        std::memcpy(b.components, components + (order - halfOrder), halfOrder * sizeof(T));
        return b;
    };


    /**
     * @brief Construct nion from half order nions. q = (a,b)
     * @param a The first half order nion in pair.
     * @param b The second half order nion in pair.
     * @return The nion constructed from the half order nions.
     * @note This is a convenience function for constructing a nion from pairing two half order nions.
     */
    constexpr nion<T>(const nion<T> &a, const nion<T> &b) : order(a.order + b.order) {
        this->components = (T*) malloc(order * sizeof(T));
        std::memset(this->components, 0, order * sizeof(T));

        std::memcpy(this->components, a.components, a.order * sizeof(T));
        std::memcpy(this->components + a.order, b.components, b.order * sizeof(T));
    };

    /**
     * @brief copy assignment operator
     * @param other The nion to copy.
     * @return A copy of the nion.
     * @note This is a deep copy.
     */
    constexpr nion<T> &operator=(const nion<T> &other)  noexcept {
        if (&other == this) {
            return *this;
        }
        free(this->components);
        this->order = other.order;
        this->components = (T *) malloc(order * sizeof(T));
        std::memcpy(this->components, other.components, order * sizeof(T));
        return *this;
    };

    /**
     * @brief override move assignment operator
     * @param other The nion to move.
     * @return A copy of the nion.
     * @note This is a shallow copy.
     */
    constexpr nion<T> &operator=(nion<T> &&other)  noexcept {
        if (&other != this) {
            this->order = other.order;
            this->components = other.components;
            other.components = nullptr;
        }
        return *this;
    }

    /**
     * @brief convert scalar to nion
     * @param scalar The scalar to convert to a nion.
     * @return The nion constructed from the scalar.
     * @note This is a convenience function for constructing a nion from a scalar.
     */
    template<arithmetic S>
    constexpr nion<T> &operator=(S scalar) {
        *this = nion<T>(static_cast<T>(scalar), 1);
        return *this;
    };

    /**
     * @brief overload the - operator for a nion.
     * @return The negation of the nion.
     */
    constexpr nion<T> operator-() const {
        nion<T> negated = *this;
        #pragma simd
        for (uint_fast16_t i = 0; i < order; i++) {
            negated[i] *= -1;
        }
        return negated;
    };

    /**
     * @brief overload the [] operator for a nion.
     * @param index The index of the component to get.
     * @return The component at the index.
     */
    constexpr T &operator[](uint_fast16_t index) const {
        return components[index];
    }

    /**
     * @brief Get the conjugate of the nion.
     * @return The conjugate of the nion.
     * @detail (a,b)* = (a*,-b)
     * @note This is a recursive function.
     */
    constexpr nion<T> conj() const {
        nion<T> conjugate = *this;
        // negate all components except the first
        #pragma simd
        for (uint_fast16_t i = 1; i < order; i++) {
            conjugate[i] *= -1;
        }
        return conjugate;
    };

    /**
     * @breif overload the += operator for nions.
     * @param other The nion to add to this nion.
     * @return The sum of this nion and the other nion inplace.
     */
    constexpr void operator+=(const nion <T> &other) {
        // if the order is less than the other nion, resize this nion.
        if (order < other.order)
            resize(other.order);

        // add the components of the other nion to this nion.
        #pragma simd
        for (uint_fast16_t i = 0; i < other.order; i++) {
            components[i] += other[i];
        }
    };

    /**
     * @breif overload the -= operator for nions.
     * @param other The nion to substract from this nion.
     * @return The subtraction of this nion and the other nion inplace.
     */
    constexpr void operator-=(const nion <T> &other) {
        // if the order is less than the other nion, resize this nion.
        if (this->order < other.order)
            resize(other.order);

        // substract the components of the other nion from this nion.
        #pragma simd
        for (uint_fast16_t i = 0; i < other.order; i++) {
            components[i] -= other[i];
        }
    };

    /**
     * @breif overload the *= operator for nions.
     * @param other The nion to multiply this nion by.
     * @return The product of this nion and the other nion inplace.
     */
    constexpr void operator*=(const nion <T> &other) {
        nion<T> product = *this * other;
        free(this->components);
        *this = product;
    };

    /**
     * @breif overload the /= operator for nions.
     * @param other The nion to divide this nion by.
     * @return The division of this nion and the other nion inplace.
     */
    constexpr void operator/=(const nion <T> &other) {
        nion<T> quotient = *this / other;
        free(this->components);
        *this = quotient;
    };

    /**
     * @brief overload the + operator for nions.
     * @param other The nion to add to this nion.
     * @return The sum of this nion and the other nion.
     */
    constexpr nion<T> operator+(const nion <T> &other) const {
        nion<T> sum = *this;
        sum += other;
        return sum;
    };

    /**
     * @brief overload the - operator for nions.
     * @param other The nion to substract this nion by.
     * @return The subtraction of this nion and the other nion.
     */
    constexpr nion<T> operator-(const nion <T> &other) const {
        nion<T> difference = *this;
        difference -= other;
        return difference;
    };

    /**
     * @brief overload the * operator for nions.
     * @param other The nion to multiply this nion by.
     * @return The product of this nion and the other nion.
     * @details The product of two nions is defined as (a,b)(c,d) = (a c - d* b, d a + b c*) = (z,w), where * is the conjugate.
     * and a, b, c, d are the nions with half the order of the original nions.
     * @note This is only defined for nions of even order.
     * @note product has the same order as the larger order of the two nions.
     * @note This is recursive function.
     * TODO: Replace recursion with stack based iteration.
     */
    constexpr nion<T> operator*(const nion <T> &other) const {
        // check if the orders are the same
        if (order == other.order) {

            // if the order is one, the product is the scalar product
            if (order == 1) {
                return nion<T>(real() * other.real(), 1);
            }


            // extract the dual elements of the nions with half the order
            nion<T> a = std::move(this->left());
            nion<T> b = std::move(this->right());
            nion<T> c = std::move(other.left());
            nion<T> d = std::move(other.right());

            // calculate the product
            nion<T> z = (a * c) - (d.conj() * b);
            nion<T> w = (d * a) + (b * c.conj());

            return nion<T>(z, w);
        } else { // orders are different

            // resize the nion with the smaller order to the order of the nion with the larger order.
            if (order > other.order) {
                return *this * other.resize(order);
            } else {
                return this->resize(other.order) * other;
            }
        }
    };

    /**
     * @brief compute the inverse of the nion.
     * @return The inverse of the nion.
     */
    constexpr nion<T> inv() const {
        return conj() / abs();
    };

    /**
     * @brief overload the / operator for nions.
     * @param other The nion to divide this nion by.
     * @return The division of this nion and the other nion.
     */
    constexpr nion<T> operator/(const nion <T> &other) const {
        return *this * other.inv();
    };

    /**
     * @brief compute the dot product with another nion.
     * @param other The nion to compute the dot product with.
     * @return The dot product of this nion and the other nion.
     */
    constexpr T dot(const nion<T> &other) const {
        T dotProduct = 0;
        #pragma simd
        for (uint_fast16_t i = 0; i < order; i++) {
            dotProduct += components[i] * other[i];
        }
        return dotProduct;
    };

    /**
     * @brief absolute value of the nion.
     * @return The absolute value of the nion.
     * @details The absolute value of the nion is the sum of the squares of its components.
     */
    constexpr T abs() const {
        T absVal = 0;

        #pragma simd
        for (uint_fast16_t i = 0; i < order; i++) {
            absVal += components[i] * components[i];
        }

        return absVal;
    };

    /**
     * @brief norm of the nion.
     * @return The norm of the nion.
     * @details The norm of a nion is the sqrt of the sum of the squares of its components.
     */
    constexpr T norm() const { return sqrt(abs()); };

    /**
     * @brief project nion to the real line.
     * @tparam T type of the nion.
     * @param nion The nion to project to the real line.
     * @return The projection of the nion to the real line.
     * @details computes the length of the nion and does the shortest rotation to the real line.
     */
    constexpr T proj() const { return copysign(norm(), real()); }

    /**
     * @brief return real part of the nion.
     * @return The real part of the nion.
     */
    constexpr T real() const { return components[0]; }

    /**
     * @brief Calculate the imaginary components of a nion.
     * @return The imaginary components of the nion.
     */
    constexpr nion<T> imag() const {
        nion<T> imaginary = *this;
        imaginary[0] = 0;
        return imaginary;
    }

    /**
     * @brief resize the nion to a new order.
     * @param newOrder The new order of the nion.
     * @return The nion converted to the new order.
     * @note newOrder must be larger than the current order.
     */
    constexpr nion<T> resize(uint_fast16_t newOrder) const {

        // keep the same order if the new order is smaller than the current order
        if (newOrder <= order) {
            return *this;
        }

        nion<T> converted(newOrder);
        std::memcpy(converted.components, components, sizeof(T) * order);
        return converted;
    };


    /**
     * @brief overload the == operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if the nions are equal, false otherwise.
     * @details Two nions are equal if they have the same order and the same components.
     */
    constexpr bool operator==(const nion <T> &other) const {
        if (order != other.order) {
            return false;
        }
        #pragma simd
        for (uint_fast16_t i = 0; i < order; i++) {
            if (std::abs(components[i] - other[i]) >= epsilon) {
                return false;
            }
        }
        return true;
    };

    /**
     * @brief overload the != operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if the nions are not equal, false otherwise.
     * @details Two nions are equal if they have the same order and the same components.
     */
    constexpr bool operator!=(const nion <T> &other) const {
        if (order != other.order) {
            return true;
        }
        #pragma simd
        for (uint_fast16_t i = 0; i < order; i++) {
            if (std::abs(components[i] - other[i]) >= epsilon) {
                return true;
            }
        }
        return false;
    };

    /**
     * @brief overload the > operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if this nion is greater than the other nion, false otherwise.
     * @details sorting is undefined for nions with orders greater than 1. However, we can still compare
     *         nions with orders greater than 1 by comparing the projections of the nions onto the real line.
     */
    constexpr bool operator>(const nion <T> &other) const {
        return proj() > other.proj();
    };

    /**
     * @brief overload the < operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if this nion is less than the other nion, false otherwise.
     * @details sorting is undefined for nions with orders greater than 1. However, we can still compare
     *         nions with orders greater than 1 by comparing the projections of the nions onto the real line.
     */
    constexpr bool operator<(const nion <T> &other) const {
        return proj() < other.proj();
    };

    /**
     * @brief overload the >= operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if this nion is greater than or equal to the other nion, false otherwise.
     * @details sorting is undefined for nions with orders greater than 1. However, we can still compare
     *         nions with orders greater than 1 by comparing the projections of the nions onto the real line.
     */
    constexpr bool operator>=(const nion <T> &other) const {
        return proj() >= other.proj();
    };

    /**
     * @brief overload the <= operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if this nion is less than or equal to the other nion, false otherwise.
     * @details sorting is undefined for nions with orders greater than 1. However, we can still compare
     *         nions with orders greater than 1 by comparing the projections of the nions onto the real line.
     */
    constexpr bool operator<=(const nion <T> &other) const {
        return proj() <= other.proj();
    };

    /**
     * @brief overload the + operator for nions with scalars.
     * @tparam S The type of the scalar.
     * @param scalar The scalar to add this nion by.
     * @return The sum of this nion and the scalar.
     */
    template<arithmetic S>
    constexpr nion<T> operator+(S scalar) const {
        nion<T> sum = *this;
        sum.components[0] += scalar;
        return sum;
    };

    /**
     * @brief overload the - operator for nions with scalars.
     * @tparam S The type of the scalar.
     * @param scalar The scalar to subtract this nion by.
     * @return The subtraction of this nion and the scalar.
     */
    template<arithmetic S>
    constexpr nion<T> operator-(S scalar) const {
        return *this + (-scalar);
    };

    /**
     * @brief overload the * operator for nions with scalars.
     * @tparam S The type of the scalar.
     * @param scalar The scalar to multiply this nion by.
     * @return The product of this nion and the scalar.
     */
    template<arithmetic S>
    constexpr nion<T> operator*(S scalar) const {
        nion<T> product = *this;
        #pragma simd
        for (uint_fast16_t i = 0; i < order; i++) {
            product[i] *= scalar;
        }
        return product;
    };

    /**
     * @brief overload the / operator for nions with scalars.
     * @tparam S The type of the scalar.
     * @param scalar The scalar to divide this nion by.
     * @return The division of this nion and the scalar.
     */
    template<arithmetic S>
    constexpr nion<T> operator/(S scalar) const {
        return *this * (1.0l / scalar);
    };

    /**
     * @brief overload the += operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param other The scalar to add to this nion.
     * @return The sum of this nion and the scalar inplace.
     */
    template<arithmetic S>
    constexpr void operator+=(S scalar) const {
        *this[0] += scalar;
    };

    /**
     * @brief overload the -= operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param other The scalar to substract from this nion.
     * @return The subtraction of this nion and the scalar inplace.
     */
    template<arithmetic S>
    constexpr void operator-=(S scalar) const {
        *this[0] -= scalar;
    };

    /**
     * @breif overload the *= operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param scalar The scalar to multiply this nion by.
     * @return The product of this nion and the scalar inplace.
     */
    template<arithmetic S>
    constexpr void operator*=(S scalar) {
        #pragma simd
        for (uint_fast16_t i = 0; i < order; i++) {
            components[i] *= scalar;
        }
    };

    /**
     * @breif overload the /= operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param scalar The scalar to divide this nion by.
     * @return The division of this nion and the scalar inplace.
     */
    template<arithmetic S>
    constexpr void operator/=(S scalar) {
        # pragma simd
        for (uint_fast16_t i = 0; i < order; i++) {
            components[i] /= scalar;
        }
    };

    /**
     * @brief overload the == operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param scalar The scalar to compare this nion to.
     * @return True if the nion is equal to the scalar, false otherwise.
     * @details A nion is equal to a scalar if the scalar is equal to the first component of the nion
     *          and all other components are equal to zero.
     */
    template<arithmetic S>
    constexpr bool operator==(S scalar) const {
        if (real() != scalar) {
            return false;
        }
        #pragma simd
        for (uint_fast16_t i = 1; i < order; i++) {
            if (components[i] != 0) {
                return false;
            }
        }
        return true;
    };

    /**
     * @brief overload the != operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param scalar The scalar to compare this nion to.
     * @return False if the nion is equal to the scalar, true otherwise.
     * @details A nion is equal to a scalar if the scalar is equal to the first component of the nion
     *          and all other components are equal to zero.
     */
    template<arithmetic S>
    constexpr bool operator!=(S scalar) const {
        if (real() != scalar) {
            return true;
        }
        #pragma simd
        for (uint_fast16_t i = 1; i < order; i++) {
            if (components[i] != 0) {
                return true;
            }
        }
        return false;
    };
};


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
template<arithmetic T, arithmetic S>
constexpr nion<T> operator*(S scalar, const nion<T> &z) {
    return z * scalar;
}

/**
 * @brief overload the / operator for lhs scalars and rhs nions.
 * @tparam T type of the nion.
 * @tparam S type of the scalar.
 * @param scalar type of the scalar.
 * @param z The nion to divide the scalar by.
 * @return
 */
template<arithmetic T, arithmetic S>
constexpr nion<T> operator/(S scalar, const nion<T> &z) {
    return z.inv() * scalar;
}

/**
* @brief overload the + operator for lhs scalars and rhs nions.
 * @tparam S type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to add the scalar by.
*/
template<arithmetic T, arithmetic S>
constexpr nion<T> operator+(S scalar, const nion<T> &z) {
    return z + scalar;
}

/**
 * @brief overload the - operator for lhs scalars and rhs nions.
 * @tparam S The type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to subtract the scalar by.
*/
template<arithmetic T, arithmetic S>
constexpr nion<T> operator-(S scalar, const nion<T> &z) {
    return -z + scalar;
}

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
template<arithmetic T, arithmetic S>
bool operator==(S scalar, const nion<T> &z) {
    return z == scalar;
}

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
template<arithmetic T, arithmetic S>
bool operator!=(S scalar, const nion<T> &z) {
    return z != scalar;
}

/**
     * @brief overload the << operator for nions.
     * @param os The output stream.
     * @param z The nion to print.
     * @return The output stream.
     */
    template<arithmetic T>
std::ostream &operator<<(std::ostream &os, const nion<T> &z) {
    T component = z.components[0];
    os << component;

    for (int i = 1; i < z.order; i++) {
        component = z.components[i];
        os << " + " << component << " e" << i;
    }
    return os;
}

/**
 * @brief overload the >> operator for nions.
 * @param is The input stream.
 * @param z The nion to read into.
 * @return The input stream.
 */
template<arithmetic T>
std::istream &operator>>(std::istream &is, nion<T> &z) {
    for (int i = 0; i < z.order; i++) {
        is >> z.components[i];
    }
    return is;
}

/***************************
    *  NION FUNCTION IMPLEMENTATIONS *
    ***************************/

/**
 * @brief Calculate the real part of a nion.
 * @tparam T type of the nion.
 * @param z The nion to calculate the real part of.
 * @return The real part of the nion.
 */
template<arithmetic T>
T real(const nion<T> &z) {
    return z.real();
}

/**
 * @brief Calculate the imaginary part of a nion.
 * @tparam T type of the nion.
 * @param z The nion to calculate the imaginary part of.
 * @return The imaginary part of the nion.
 */
template<arithmetic T>
constexpr nion<T> imag(const nion<T> &z) {
    return z.imag();
}

/**
 * @brief compute the conjugate of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the conjugate of.
 * @return The conjugate of the nion.
 */
template<arithmetic T>
constexpr nion<T> conj(const nion<T> &z) {
    return z.conj();
}

/**
 * @brief compute the absolute value of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the absolute value of.
 * @return The absolute value of the nion.
 */
template<arithmetic T>
T abs(const nion<T> &z) {
    return z.abs();
}

/**
 * @brief compute the norm of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the norm of.
 * @return The norm of the nion.
 */
template<arithmetic T>
T norm(const nion<T> &z) {
    return z.norm();
}

/**
 * @brief compute the inverse of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse of.
 * @return The inverse of the nion.
 */
template<arithmetic T>
constexpr nion<T> inv(const nion<T> &z) {
    return z.inv();
}

/**
 * @brief compute the dot product of two nions.
 * @tparam T type of the nions.
 * @param lhs The left hand side nion.
 * @param rhs The right hand side nion.
 * @return The dot product of the nions.
 */
template<arithmetic T>
T dot(const nion<T> &lhs, const nion<T> &rhs) {
    return lhs.dot(rhs);
}

/**
 * @brief compute the cross product of two nions.
 * @tparam T type of the nions.
 * @param lhs The left hand side nion.
 * @param rhs The right hand side nion.
 * @return The cross product of the nions.
 */
template<arithmetic T>
nion<T> cross(const nion<T> &lhs, const nion<T> &rhs); //TODO: implement cross product

/**
 * @brief compute the wedge product of two nions.
 * @tparam T type of the nions.
 * @param lhs The left hand side nion.
 * @param rhs The right hand side nion.
 * @return The wedge product of the nions.
 */
template<arithmetic T>
nion<T> wedge(const nion<T> &lhs, const nion<T> &rhs); //TODO: implement wedge product

/**
 * @brief compute the outer product of two nions.
 * @tparam T type of the nions.
 * @param lhs The left hand side nion.
 * @param rhs The right hand side nion.
 * @return The outer product of the nions.
 */
template<arithmetic T>
nion<T> outer(const nion<T> &lhs, const nion<T> &rhs); //TODO: implement outer product


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
template<arithmetic T>
constexpr nion<T> exp(const nion<T> &z) noexcept {

    // get polar form of nion
    T r = z.real();
    nion<T> v = z.imag();

    // make unit vector
    T v_norm = v.norm();
    if (v_norm != 0)
        v /= v_norm;

    // compute exponential of nion
    return exp(r) * (v * sin(v_norm) + cos(v_norm));
}

/**
 * @brief compute the natural logarithm of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the natural logarithm of.
 * @return The natural logarithm of the nion.
 * @details The natural logarithm of a nion is defined as ln(z) = ln(|z|) + v/|v| * arccos(r/|z|).
 *          where a is the real component and v is the imaginary components.
 */
template<arithmetic T>
constexpr nion<T> log(const nion<T> &z) noexcept {

    // get polar form of nion
    T r = z.real();
    nion<T> v = z.imag();

    // compute norms
    T z_norm = z.norm();
    T v_norm = v.norm();

    // make unit vector
    if (v_norm != 0)
        v /= v_norm;

    // compute natural logarithm of nion
    return log(z_norm) + v * atan2(v_norm, r);
}

/**
 * @brief compute the power of a nion with scalar.
 * @tparam T type of the nion.
 * @tparam S type of the power.
 * @param z The nion to compute the power of.
 * @param power The power to raise the nion to.
 * @return The power of the nion.
 * @details The power of a nion is defined as z^p = e^(p * ln(z)).
 */
template<arithmetic T, arithmetic S>
constexpr nion<T> pow(const nion<T> &z, S power) noexcept {
    if (power == 0) {
        return nion<T>(z.order, 1);
    } else if (std::abs(power - 1) <= z.epsilon) {
        return nion<T>(z);
    } else if (std::abs(power + 1) <= z.epsilon) {
        return inv(z);
    } else {
        return exp(power * log(z));
    }
}

/**
 * @brief compute the power of a nion with nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the power of.
 * @param power The power to raise the nion to.
 * @return The power of the nion.
 * @details The power of a nion is defined as z^p = e^(p * ln(z)).
 */
template<arithmetic T>
constexpr nion<T> pow(const nion<T> &z, const nion<T> &power) noexcept {
    return exp(power * log(z));
}

/**
 * @brief compute the power of a scalar with a nion.
 * @tparam T type of the nion.
 * @tparam S type of the power.
 * @param z The nion to use as a power.
 * @param x The scalar to raise by the nion
 * @return The power of the nion.
 * @details The power of with a nion is defined as x^z = e^(z * ln(x)).
 */
template<arithmetic T, arithmetic S>
constexpr nion<T> pow(S x, const nion<T> &z) noexcept {
    if (x == 0) {
        return nion<T>(z.order, 0);
    } else if (std::abs(x - 1) <= z.epsilon) {
        return nion<T>(z.order, 1);
    } else {
        return exp(z * log(x));
    }
}

/**
 * @brief compute the square root of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the square root of.
 * @return The square root of the nion.
 * @details The square root of a nion is defined as sqrt(z) = z^(1/2).
 */
template<arithmetic T>
constexpr nion<T> sqrt(const nion<T> &z) noexcept{
    return pow(z, 0.5f);
}

/**
 * @brief compute the cube root of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the cube root of.
 * @return The cube root of the nion.
 * @details The cube root of a nion is defined as cbrt(z) = z^(1/3).
 */
template<arithmetic T>
constexpr nion<T> cbrt(const nion<T> &z) noexcept{
    return pow(z, 1.0/3.0f);
}

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
template<arithmetic T>
constexpr nion<T> sinh(const nion<T> &z) noexcept{
    return (exp(z) - exp(-z)) * 0.5f;
}

/**
 * @brief compute the hyperbolic cosine of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic cosine of.
 * @return The hyperbolic cosine of the nion.
 * @details The hyperbolic cosine of a nion is defined as cosh(z) = (e^z + e^-z) / 2.
 */
template<arithmetic T>
constexpr nion<T> cosh(const nion<T> &z) noexcept{
    return (exp(z) + exp(-z)) * 0.5f;
}

/**
 * @brief compute the hyperbolic tangent of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic tangent of.
 * @return The hyperbolic tangent of the nion.
 * @details The hyperbolic tangent of a nion is defined as tanh(z) = sinh(z) / cosh(z).
 */
template<arithmetic T>
constexpr nion<T> tanh(const nion<T> &z) noexcept{
    return sinh(z) / cosh(z);
}

/**
 * @brief compute the hyperbolic cotangent of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic cotangent of.
 * @return The hyperbolic cotangent of the nion.
 * @details The hyperbolic cotangent of a nion is defined as coth(z) = cosh(z) / sinh(z).
 */
template<arithmetic T>
constexpr nion<T> coth(const nion<T> &z) noexcept{
    return cosh(z) / sinh(z);
}

/**
 * @brief compute the hyperbolic secant of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic secant of.
 * @return The hyperbolic secant of the nion.
 * @details The hyperbolic secant of a nion is defined as sech(z) = 1 / cosh(z).
 */
template<arithmetic T>
constexpr nion<T> sech(const nion<T> &z) noexcept{
    return cosh(z).inv();
}

/**
 * @brief compute the hyperbolic cosecant of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic cosecant of.
 * @return The hyperbolic cosecant of the nion.
 * @details The hyperbolic cosecant of a nion is defined as csch(z) = 1 / sinh(z).
 */
template<arithmetic T>
constexpr nion<T> csch(const nion<T> &z) noexcept{
    return sinh(z).inv();
}

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
template<arithmetic T>
constexpr nion<T> sin(const nion<T> &z) noexcept{
    // get the polar form of the nion
    T r = real(z);
    nion<T> v = imag(z);

    // make unit vector
    T v_norm = v.norm();
    if (v_norm != 0)
        v /= v_norm;

    // compute the sine of the nion
    return sin(r) * cosh(v_norm) + v * sinh(v_norm) * cos(r);
}

/**
 * @brief compute the cosine of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the cosine of.
 * @return The cosine of the nion.
 * @details The cosine of the nion is defined as cos(z) = cos(r + v) = cos(r) * cosh(|v|) - sin(r) * sinh(|v|) * v/|v|.
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://en.wikipedia.org/wiki/Sine_and_cosine#Relationship_to_complex_numbers
 */
template<arithmetic T>
constexpr nion<T> cos(const nion<T> &z) noexcept{
    // get the polar form of the nion
    T r = real(z);
    nion<T> v = imag(z);

    // make unit vector
    T v_norm = v.norm();
    if (v_norm != 0)
        v /= v_norm;

    // compute the cosine of the nion
    return cos(r) * cosh(v_norm) - v * sinh(v_norm) * sin(r);
}

/**
 * @brief compute the tangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the tangent of.
 * @return The tangent of the nion.
 * @details The tangent of the nion is defined as tan(z) = sin(z) / cos(z).
 */
template<arithmetic T>
constexpr nion<T> tan(const nion<T> &z) noexcept{
    return sin(z) / cos(z);
}

/**
 * @brief compute the cotangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the cotangent of.
 * @return The cotangent of the nion.
 * @details The cotangent of the nion is defined as cot(z) = cos(z) / sin(z).
 */
template<arithmetic T>
constexpr nion<T> cot(const nion<T> &z) noexcept{
    return cos(z) / sin(z);
}

/**
 * @brief compute the secant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the secant of.
 * @return The secant of the nion.
 * @details The secant of the nion is defined as sec(z) = 1 / cos(z).
 */
template<arithmetic T>
constexpr nion<T> sec(const nion<T> &z) noexcept{
    return cos(z).inv();
}

/**
 * @brief compute the cosecant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the cosecant of.
 * @return The cosecant of the nion.
 * @details The cosecant of the nion is defined as csc(z) = 1 / sin(z).
 */
template<arithmetic T>
constexpr nion<T> csc(const nion<T> &z) noexcept{
    return sin(z).inv();
}

/***************************
    *  NION INVERSE HYPERBOLIC TRIGONOMETRIC FUNCTIONS *
    ***************************/

/**
 * @brief compute the inverse hyperbolic sine of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse hyperbolic sine of.
 * @return The inverse hyperbolic sine of the nion.
 * @details The inverse hyperbolic sine of the nion is defined as asinh(z) = asinh(r + v) = ln(z + sqrt(1 + z^2)).
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://mathworld.wolfram.com/InverseHyperbolicSine.html
 */
template<arithmetic T>
constexpr nion<T> asinh(const nion<T> &z) noexcept{
    return log(z + sqrt(1 + z*z));
}

/**
 * @brief compute the inverse hyperbolic cosine of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse hyperbolic cosine of.
 * @return The inverse hyperbolic cosine of the nion.
 * @details The inverse hyperbolic cosine of the nion is defined as acosh(z) = acosh(r + v) = ln(z + sqrt(z + 1)*sqrt(z - 1)).
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://mathworld.wolfram.com/InverseHyperbolicCosine.html
 */
template<arithmetic T>
constexpr nion<T> acosh(const nion<T> &z) noexcept{
    return log(z + sqrt(z + 1) * sqrt(z - 1));
}

/**
 * @brief compute the inverse hyperbolic tangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse hyperbolic tangent of.
 * @return The inverse hyperbolic tangent of the nion.
 * @details The inverse hyperbolic tangent of the nion is defined as atanh(z) = atanh(r + v) = 1/2 * (ln(1 + z) - ln(1 - z))).
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://mathworld.wolfram.com/InverseHyperbolicTangent.html
 */
template<arithmetic T>
constexpr nion<T> atanh(const nion<T> &z) noexcept{
    return 0.5l * (log(1 + z) - log(1 - z));
}

/**
 * @brief compute the inverse hyperbolic cotangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse hyperbolic cotangent of.
 * @return The inverse hyperbolic cotangent of the nion.
 * @details The inverse hyperbolic cotangent of the nion is defined as acoth(z) = acoth(r + v) = 1/2 * (ln(1 + 1/z) - ln(1 - 1/z)).
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://mathworld.wolfram.com/InverseHyperbolicCotangent.html
 */
template<arithmetic T>
constexpr nion<T> acoth(const nion<T> &z) noexcept{
    return 0.5l * (log(1 + 1/z) - log(1 - 1/z));
}

/**
 * @brief compute the inverse hyperbolic secant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse hyperbolic secant of.
 * @return The inverse hyperbolic secant of the nion.
 * @details The inverse hyperbolic secant of the nion is defined as asech(z) = asech(r + v) = ln(sqrt(1/z - 1)*sqrt(1/z + 1) + 1/z).
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://mathworld.wolfram.com/InverseHyperbolicSecant.html
 */
template<arithmetic T>
constexpr nion<T> asech(const nion<T> &z) noexcept{
    return log(sqrt(1/z - 1) * sqrt(1/z + 1) + 1/z);
}

/**
 * @brief compute the inverse hyperbolic cosecant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse hyperbolic cosecant of.
 * @return The inverse hyperbolic cosecant of the nion.
 * @details The inverse hyperbolic cosecant of the nion is defined as acsch(z) = acsch(r + v) = ln(sqrt(1+1/z^2) + 1/z).
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://mathworld.wolfram.com/InverseHyperbolicCosecant.html
 */
template<arithmetic T>
constexpr nion<T> acsch(const nion<T> &z) noexcept{
    return log(sqrt(1 + pow(z, -2)) + 1/z);
}

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
template<arithmetic T>
constexpr nion<T> asin(const nion<T> &z) noexcept {
    // get the polar form of the nion
    nion<T> v = z.imag();

    // make unit vector
    T v_norm = v.norm();
    if (v_norm != 0)
        v /= v_norm;

    // compute the inv sine of the nion
    return v * log(sqrt(1 - z*z) - v * z);
}

/**
 * @brief compute the inverse cosine of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse cosine of.
 * @return The inverse cosine of the nion.
 * @details The inverse cosine of the nion is defined as acos(z) = pi/2 - asin(z).
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Extension_to_complex_plane
 */
template<arithmetic T>
constexpr nion<T> acos(const nion<T> &z) noexcept {
    return M_PI_2l - asin(z);
}

/**
 * @brief compute the inverse tangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse tangent of.
 * @return The inverse tangent of the nion.
 * @details The inverse tangent of the nion is defined as atan(z) = v * atanh(-v*z);
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Logarithmic_forms
 */
template<arithmetic T>
constexpr nion<T> atan(const nion<T> &z) noexcept {
    // get the polar form of the nion
    nion<T> v = z.imag();

    // make unit vector
    T v_norm = v.norm();
    if (v_norm != 0)
        v /= v_norm;

    // compute the inv tangent of the nion
    return v * atanh(-v*z);
}

/**
 * @brief compute the inverse cotangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse cotangent of.
 * @return The inverse cotangent of the nion.
 * @details The inverse cotangent of the nion is defined as acot(z) = pi/2 - atan(z).
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Extension_to_complex_plane
 */
template<arithmetic T>
constexpr nion<T> acot(const nion<T> &z) noexcept {
    return M_PI_2l - atan(z);
}

/**
 * @brief compute the inverse secant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse secant of.
 * @return The inverse secant of the nion.
 * @details The inverse secant of the nion is defined as asec(z) = acos(1/z).
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Extension_to_complex_plane
 */
template<arithmetic T>
constexpr nion<T> asec(const nion<T> &z) noexcept {
    return acos(1/z);
}

/**
 * @brief compute the inverse cosecant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse cosecant of.
 * @return The inverse cosecant of the nion.
 * @details The inverse cosecant of the nion is defined as acsc(z) = asin(1/z).
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Extension_to_complex_plane
 */
template<arithmetic T>
constexpr nion<T> acsc(const nion<T> &z) noexcept {
    return asin(1/z);
}

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
template<arithmetic T>
constexpr nion<T> gamma(const nion<T> &z) noexcept {
    // compute the gamma function of the nion
    return sqrt(2 * M_PI) * exp(-z) * sqrt(1/z) * pow(1/(12 * z - 1/(10 * z)) + z, z);
}

#endif //NION_NION_HPP