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
 * @brief evaluates if two values are similar to machine precision
 * @tparam T type of first value
 * @tparam S  type of second value
 * @param a  first value
 * @param b  second value
 * @param epsilon tolerance
 * @return true if similar, false otherwise
 */
template<typename T, typename D = uint_fast16_t, typename S = T>
constexpr inline bool value_is_similar(const T a, const S b, const T epsilon = std::numeric_limits<T>::epsilon()) {
    return std::fabs(a - static_cast<T>(b)) <= epsilon;
}


/**
 * @file nion.hpp
 * @brief Templated class that implements Cayley-Dickson division algebra.
 * @author Marcus Dante Liebenthal
 * @version 1.0
 * @date 2022-10-08
 */
template<typename T, typename D = uint_fast16_t>
struct nion {
    T* components;
    D degree;
    static_assert(std::is_arithmetic_v<T>, "nion only supports arithmetic types");

    /**
     * @brief Default constructor.
     * @details Constructs a null nion object.
     */
    constexpr inline nion<T, D>() : components(nullptr), degree(0) {};

    /**
     * @brief Construct a new nion object from vector
     * @param components The components of the nion.
     * @param degree The degree of the nion.
     */
    constexpr inline explicit nion<T, D>(T* components, D degree) : degree(degree) {

        // check if the degree is greater than zero
        if (degree <= 0) {
            throw std::invalid_argument("The degree of the nion must be greater than zero.");
        }

        this->components = (T*) malloc(degree * sizeof(T));
        memcpy(this->components, components, degree * sizeof(T));
    };

    /**
     * @brief Construct a new nion object from brace initializer list
     * @param degree The degree of the nion.
     */
    constexpr inline nion<T, D>(const std::initializer_list<T> &components) : degree(components.size()) {
        if (degree <= 0) {
            throw std::invalid_argument("The degree of the nion must be greater than zero.");
        }
        this->components = (T*) malloc(degree * sizeof(T));
        memcpy(this->components, components.begin(), degree * sizeof(T));
    };

    /**
     * @brief Construct an empty nion object
     * @param degree The degree of the nion.
     */
    constexpr inline explicit nion<T, D>(D degree) : degree(degree) {
        // check if the degree is greater than zero
        if (degree <= 0) {
            throw std::invalid_argument("The degree of the nions must be greater than zero.");
        }
        this->components = (T*) calloc(this->degree, sizeof(T));
    };

    /**
     * @brief copy constructor
     * @param other The nion to copy.
     * @return A copy of the nion.
     * @note This is a deep copy.
     */
    constexpr inline nion<T, D>(const nion<T, D> &other) : degree(other.degree) {
        this->components = (T*) malloc(other.degree * sizeof(T));
        memcpy(this->components, other.components, other.degree * sizeof(T));
    };

    /**
     * @brief Construct a new nion object from a scalar with no imaginary components.
     * @param degree The degree of the nion.
     * @param scalar The scalar value of the nion.
     * @return nion<T, D> The nion object.
     * @note This is a convenience function for creating a nion from a scalar.
     */
    template<typename S = T>
    constexpr inline explicit nion<T, D>(S realVal, D degree) : degree(degree) {
        // check if the degree is greater than zero
        if (degree <= 0) {
            throw std::invalid_argument("The degree of the nions must be greater than zero.");
        }
        this->components = (T*) calloc(degree, sizeof(T));
        components[0] = static_cast<T>(realVal);
    };

    /**
     * @brief Construct a new nion object from std::complex
     * @param degree The degree of the nion.
     * @param complex The std::complex object.
     * @return nion<T, D> The nion object.
     */
    constexpr inline explicit nion<T, D>(std::complex<T> complex) : degree(2) {
        this->components = (T*) calloc(2, sizeof(T));

        components[0] = complex.real();
        components[1] = complex.imag();
    };

    /**
     * @brief Destroy the nion object
     */
    ~nion<T, D>() { free(components); components = nullptr; };

    /**
     * @brief Construct nion from half degree nions. q = (a,b)
     * @param a The first half degree nion in pair.
     * @param b The second half degree nion in pair.
     * @return The nion constructed from the half degree nions.
     * @note This is a convenience function for constructing a nion from pairing two half degree nions.
     */
    constexpr inline nion<T, D>(const nion<T, D> &a, const nion<T, D> &b) : degree(a.degree + b.degree) {
        this->components = (T*) malloc(degree * sizeof(T));
        std::memset(this->components, 0, degree * sizeof(T));

        memcpy(this->components, a.components, a.degree * sizeof(T));
        memcpy(this->components + a.degree, b.components, b.degree * sizeof(T));
    };

    /**
     * @brief copy assignment operator
     * @param other The nion to copy.
     * @return A copy of the nion.
     * @note This is a deep copy.
     */
    constexpr inline nion<T, D> &operator=(const nion<T, D> &other)  noexcept {
        if (&other == this) {
            return *this;
        }

        this->degree = other.degree;
        if (this->components == nullptr) {
            this->components = (T *) malloc(other.degree * sizeof(T));
        } else {
            this->components = (T *) realloc(this->components, other.degree * sizeof(T));
        }

        memcpy(this->components, other.components, degree * sizeof(T));
        return *this;
    };

    /**
     * @brief override move assignment operator
     * @param other The nion to move.
     * @return A copy of the nion.
     * @note This is a shallow copy.
     */
    constexpr inline nion<T, D> &operator=(nion<T, D> &&other)  noexcept {
        if (&other != this) {
            free(this->components);
            this->degree = other.degree;
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
    template<typename S>
    constexpr inline nion<T, D> &operator=(S scalar) {
        free(this->components);
        this->components = nullptr;
        *this = nion<T, D>(static_cast<T>(scalar), 1);
        return *this;
    };

    /**
     * @brief overload the - operator for a nion.
     * @return The negation of the nion.
     */
    constexpr inline nion<T, D> operator-() const {
        nion<T, D> negated;
        negated.degree = this->degree;
        negated.components = (T*) malloc(this->degree * sizeof(T));

        #pragma vector aligned
        for (D i = 0; i < degree; i++) {
            negated[i] = -components[i];
        }

        return negated;
    };

    /**
     * @brief overload the [] operator for a nion.
     * @param index The index of the component to get.
     * @return The component at the index.
     */
    constexpr inline T &operator[](D index) const {
        return components[index];
    }

    /**
     * @brief Get the conjugate of the nion.
     * @return The conjugate of the nion.
     * @detail (a,b)* = (a*,-b)
     * @note This is a recursive function.
     */
    constexpr inline nion<T, D> conj() const {
        nion<T, D> conjugate;
        conjugate.degree = this->degree;
        conjugate.components = (T*) malloc(this->degree * sizeof(T));

        // negate all components except the first
        conjugate[0] = this->components[0];
        #pragma vector aligned
        for (D i = 1; i < degree; i++) {
            conjugate[i] = -components[i];
        }

        return conjugate;
    };

    /**
     * @breif overload the += operator for nions.
     * @param other The nion to add to this nion.
     * @return The sum of this nion and the other nion inplace.
     */
    constexpr inline void operator+=(const nion <T> &other) {
        // if the degree is less than the other nion, resize this nion.
        if (degree < other.degree)
            resize(other.degree);

        // add the components of the other nion to this nion.
        #pragma vector aligned
        for (D i = 0; i < other.degree; i++) {
            components[i] += other[i];
        }
    };

    /**
     * @breif overload the -= operator for nions.
     * @param other The nion to substract from this nion.
     * @return The subtraction of this nion and the other nion inplace.
     */
    constexpr inline void operator-=(const nion <T> &other) {
        // if the degree is less than the other nion, resize this nion.
        if (this->degree < other.degree)
            resize(other.degree);

        // substract the components of the other nion from this nion.
        #pragma vector aligned
        for (D i = 0; i < other.degree; i++) {
            components[i] -= other[i];
        }
    };

    /**
     * @breif overload the *= operator for nions.
     * @param other The nion to multiply this nion by.
     * @return The product of this nion and the other nion inplace.
     */
    constexpr inline void operator*=(const nion <T> &other) {
        *this = *this * other;
    };

    /**
     * @breif overload the /= operator for nions.
     * @param other The nion to divide this nion by.
     * @return The division of this nion and the other nion inplace.
     */
    constexpr inline void operator/=(const nion <T> &other) {
        *this = *this / other;
    };

    /**
     * @brief overload the + operator for nions.
     * @param other The nion to add to this nion.
     * @return The sum of this nion and the other nion.
     */
    constexpr inline nion<T, D> operator+(const nion <T> &other) const {
        nion<T, D> sum = *this;
        sum += other;
        return sum;
    };

    /**
     * @brief overload the - operator for nions.
     * @param other The nion to substract this nion by.
     * @return The subtraction of this nion and the other nion.
     */
    constexpr inline nion<T, D> operator-(const nion <T> &other) const {
        nion<T, D> difference = *this;
        difference -= other;
        return difference;
    };


    /**
     * @brief overload the * operator for nions.
     * @param other The nion to multiply this nion by.
     * @return The product of this nion and the other nion.
     * @details The product of two nions is defined as (a,b)(c,d) = (a c - d* b, d a + b c*) = (z,w), where * is the conjugate.
     * and a, b, c, d are the nions with half the degree of the original nions.
     * @note product has the same degree as the larger degree of the two nions.
     * @note This is recursive function.
     */
    constexpr inline nion<T, D> operator*(const nion <T> &other) const {

        switch (degree){
            case 1:
                // if the degree is 1, then the product is just the scalar product.
                return other * components[0];
            case 2:
                // hard-coded traditional complex product
                if (other.degree == 2)
                    return nion<T, D>({components[0] * other[0] - components[1] * other[1],
                                    components[1] * other[0] + components[0] * other[1]});
                break;
            default:
                if (other.degree == 1) // if the degree is 1, then the product is just the scalar product.
                    return *this * other[0];
                break;
        }


        // extract the dual elements of the nions with half the degree
        D half_degree = degree / 2;
        nion<T, D> a(this->components,  half_degree);
        nion<T, D> b(this->components + half_degree, degree - half_degree);
        nion<T, D> c(other.components,  half_degree);
        nion<T, D> d(other.components + half_degree, degree - half_degree);


        // calculate the product
        return nion<T, D>(
                (a * c) - (d.conj() * b), // consider adding involution parameter for sign
                (d * a) + (b * c.conj())
                );
    };

    /**
     * @brief compute the inverse of the nion.
     * @return The inverse of the nion.
     */
    constexpr inline nion<T, D> inv() const {
        return conj() / abs();
    };

    /**
     * @brief overload the / operator for nions.
     * @param other The nion to divide this nion by.
     * @return The division of this nion and the other nion.
     */
    constexpr inline nion<T, D> operator/(const nion <T> &other) const {
        return *this * other.inv();
    };

    /**
     * @brief absolute value of the nion.
     * @return The absolute value of the nion.
     * @details The absolute value of the nion is the sum of the squares of its components.
     */
    constexpr inline T abs() const {
        T absVal = 0;

        #pragma vector aligned
        for (D i = 0; i < degree; i++) {
            absVal += components[i] * components[i];
        }

        return absVal;
    };

    /**
     * @brief norm of the nion.
     * @return The norm of the nion.
     * @details The norm of a nion is the sqrt of the sum of the squares of its components.
     */
    constexpr inline T norm() const { return sqrt(abs()); };

    /**
     * @brief project nion to the real line.
     * @tparam T type of the nion.
     * @param nion The nion to project to the real line.
     * @return The projection of the nion to the real line.
     * @details computes the length of the nion and does the shortest rotation to the real line.
     */
    constexpr inline T proj() const { return copysign(norm(), real()); }

    /**
     * @brief return real part of the nion.
     * @return The real part of the nion.
     */
    constexpr inline T real() const { return components[0]; }

    /**
     * @brief Calculate the imaginary components of a nion.
     * @return The imaginary components of the nion.
     */
    constexpr inline nion<T, D> imag() const {
        nion<T, D> imag;

        T* imag_components = (T*) malloc(degree * sizeof(T));
        imag.components = imag_components;
        imag.degree = degree;

        imag[0] = 0;
        memcpy(imag_components + 1, components + 1, (degree - 1) * sizeof(T));

        return imag;
    }

    /**
     * @brief resize the nion to a new degree.
     * @param newDegree The new degree of the nion.
     * @return The nion converted to the new degree.
     * @note newDegree must be larger than the current degree.
     */
    constexpr inline void resize(D newDegree) {

        degree = newDegree;
        components = (T *) realloc(components, newDegree * sizeof(T));

        // set the new components to zero
        memset(components + degree, 0, (newDegree - degree) * sizeof(T));
    };


    /**
     * @brief overload the == operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if the nions are equal, false otherwise.
     * @details Two nions are equal if they have the same degree and the same components.
     */
    constexpr inline bool operator==(const nion <T> &other) const {
        if (degree != other.degree) {
            return false;
        }
        #pragma vector aligned
        for (D i = 0; i < degree; i++) {
            if (!value_is_similar(components[i], other[i])) {
                return false;
            }
        }
        return true;
    };

    /**
     * @brief overload the != operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if the nions are not equal, false otherwise.
     * @details Two nions are equal if they have the same degree and the same components.
     */
    constexpr inline bool operator!=(const nion <T> &other) const {
        if (degree != other.degree) {
            return true;
        }
        #pragma vector aligned
        for (D i = 0; i < degree; i++) {
            if (!value_is_similar(components[i], other[i])) {
                return true;
            }
        }
        return false;
    };

    /**
     * @brief overload the > operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if this nion is greater than the other nion, false otherwise.
     * @details sorting is undefined for nions with degrees greater than 1. However, we can still compare
     *         nions with degrees greater than 1 by comparing the projections of the nions onto the real line.
     */
    constexpr inline bool operator>(const nion <T> &other) const {
        return proj() > other.proj();
    };

    /**
     * @brief overload the < operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if this nion is less than the other nion, false otherwise.
     * @details sorting is undefined for nions with degrees greater than 1. However, we can still compare
     *         nions with degrees greater than 1 by comparing the projections of the nions onto the real line.
     */
    constexpr inline bool operator<(const nion <T> &other) const {
        return proj() < other.proj();
    };

    /**
     * @brief overload the >= operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if this nion is greater than or equal to the other nion, false otherwise.
     * @details sorting is undefined for nions with degrees greater than 1. However, we can still compare
     *         nions with degrees greater than 1 by comparing the projections of the nions onto the real line.
     */
    constexpr inline bool operator>=(const nion <T> &other) const {
        return proj() >= other.proj();
    };

    /**
     * @brief overload the <= operator for nions.
     * @param other The nion to compare this nion to.
     * @return True if this nion is less than or equal to the other nion, false otherwise.
     * @details sorting is undefined for nions with degrees greater than 1. However, we can still compare
     *         nions with degrees greater than 1 by comparing the projections of the nions onto the real line.
     */
    constexpr inline bool operator<=(const nion <T> &other) const {
        return proj() <= other.proj();
    };

    /**
     * @brief overload the + operator for nions with scalars.
     * @tparam S The type of the scalar.
     * @param scalar The scalar to add this nion by.
     * @return The sum of this nion and the scalar.
     */
    template<typename S>
    constexpr inline nion<T, D> operator+(S scalar) const {
        nion<T, D> sum = *this;
        sum.components[0] += static_cast<T>(scalar);
        return sum;
    };

    /**
     * @brief overload the - operator for nions with scalars.
     * @tparam S The type of the scalar.
     * @param scalar The scalar to subtract this nion by.
     * @return The subtraction of this nion and the scalar.
     */
    template<typename S>
    constexpr inline nion<T, D> operator-(S scalar) const {
        return *this + (static_cast<T>(-scalar));
    };

    /**
     * @brief overload the * operator for nions with scalars.
     * @tparam S The type of the scalar.
     * @param scalar The scalar to multiply this nion by.
     * @return The product of this nion and the scalar.
     */
    template<typename S>
    constexpr inline nion<T, D> operator*(S scalar) const {
        nion<T, D> product = *this;
        #pragma vector aligned
        for (D i = 0; i < degree; i++) {
            product[i] *= static_cast<T>(scalar);
        }
        return product;
    };

    /**
     * @brief overload the / operator for nions with scalars.
     * @tparam S The type of the scalar.
     * @param scalar The scalar to divide this nion by.
     * @return The division of this nion and the scalar.
     */
    template<typename S>
    constexpr inline nion<T, D> operator/(S scalar) const {
        return *this * (static_cast<T>(1 / scalar));
    };

    /**
     * @brief overload the += operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param other The scalar to add to this nion.
     * @return The sum of this nion and the scalar inplace.
     */
    template<typename S>
    constexpr inline void operator+=(S scalar) const {
        components[0] += static_cast<T>(scalar);
    };

    /**
     * @brief overload the -= operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param other The scalar to substract from this nion.
     * @return The subtraction of this nion and the scalar inplace.
     */
    template<typename S>
    constexpr inline void operator-=(S scalar) const {
        components[0] -= static_cast<T>(scalar);
    };

    /**
     * @breif overload the *= operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param scalar The scalar to multiply this nion by.
     * @return The product of this nion and the scalar inplace.
     */
    template<typename S>
    constexpr inline void operator*=(S scalar) {
        #pragma vector aligned
        for (D i = 0; i < degree; i++) {
            components[i] *= static_cast<T>(scalar);
        }
    };

    /**
     * @breif overload the /= operator for nions with scalars.
     * @tparam S type of the scalar.
     * @param scalar The scalar to divide this nion by.
     * @return The division of this nion and the scalar inplace.
     */
    template<typename S>
    constexpr inline void operator/=(S scalar) {
        # pragma vector aligned
        for (D i = 0; i < degree; i++) {
            components[i] /= static_cast<T>(scalar);
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
    template<typename S>
    constexpr inline bool operator==(S scalar) const {
        if (!value_is_similar(real(), scalar)) {
            return false;
        }
        #pragma vector aligned
        for (D i = 1; i < degree; i++) {
            if (!value_is_similar(components[i], 0)) {
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
    template<typename S>
    constexpr inline bool operator!=(S scalar) const {
        if (!value_is_similar(real(), scalar)) {
            return true;
        }
        #pragma vector aligned
        for (D i = 1; i < degree; i++) {
            if (!value_is_similar(components[i], 0)) {
                return true;
            }
        }
        return false;
    };

    /**
     * @brief comparison if nion is a real number.
     * @return True if the nion is a real number, false otherwise.
     * @details A nion is a real number if all components except the first are equal to zero.
     */
    constexpr inline bool is_real() const {
        #pragma vector aligned
        for (D i = 1; i < degree; i++) {
            if (!value_is_similar(components[i], 0)) {
                return false;
            }
        }
        return true;
    }
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
template<typename T, typename D = uint_fast16_t, typename S = T>
static constexpr inline nion<T, D> operator*(S scalar, const nion<T, D> &z) {
    return z * static_cast<T>(scalar);
}

/**
 * @brief overload the / operator for lhs scalars and rhs nions.
 * @tparam T type of the nion.
 * @tparam S type of the scalar.
 * @param scalar type of the scalar.
 * @param z The nion to divide the scalar by.
 * @return
 */
template<typename T, typename D = uint_fast16_t, typename S = T>
static constexpr inline nion<T, D> operator/(S scalar, const nion<T, D> &z) {
    return z.inv() * static_cast<T>(scalar);
}

/**
* @brief overload the + operator for lhs scalars and rhs nions.
 * @tparam S type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to add the scalar by.
*/
template<typename T, typename D = uint_fast16_t, typename S = T>
static constexpr inline nion<T, D> operator+(S scalar, const nion<T, D> &z) {
    return z + static_cast<T>(scalar);
}

/**
 * @brief overload the - operator for lhs scalars and rhs nions.
 * @tparam S The type of the scalar.
 * @tparam T type of the nion.
 * @param scalar type of the scalar.
 * @param z The nion to subtract the scalar by.
*/
template<typename T, typename D = uint_fast16_t, typename S = T>
static constexpr inline nion<T, D> operator-(S scalar, const nion<T, D> &z) {
    return -z + static_cast<T>(scalar);
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
template<typename T, typename D = uint_fast16_t, typename S = T>
static constexpr inline bool operator==(S scalar, const nion<T, D> &z) {
    return z == static_cast<T>(scalar);
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
template<typename T, typename D = uint_fast16_t, typename S = T>
static constexpr inline bool operator!=(S scalar, const nion<T, D> &z) {
    return z != static_cast<T>(scalar);
}

/**
     * @brief overload the << operator for nions.
     * @param os The output stream.
     * @param z The nion to print.
     * @return The output stream.
     */
    template<typename T, typename D = uint_fast16_t>
std::ostream &operator<<(std::ostream &os, const nion<T, D> &z) {
    T component = z.components[0];
    os << component;

    for (D i = 1; i < z.degree; i++) {
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
template<typename T, typename D = uint_fast16_t>
std::istream &operator>>(std::istream &is, nion<T, D> &z) {
    for (D i = 0; i < z.degree; i++) {
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
template<typename T, typename D = uint_fast16_t>
static constexpr inline T real(const nion<T, D> &z) {
    return z.real();
}

/**
 * @brief Calculate the imaginary part of a nion.
 * @tparam T type of the nion.
 * @param z The nion to calculate the imaginary part of.
 * @return The imaginary part of the nion.
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> imag(const nion<T, D> &z) {
    return z.imag();
}

/**
 * @brief compute the conjugate of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the conjugate of.
 * @return The conjugate of the nion.
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> conj(const nion<T, D> &z) {
    return z.conj();
}

/**
 * @brief compute the absolute value of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the absolute value of.
 * @return The absolute value of the nion.
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline T abs(const nion<T, D> &z) {
    return z.abs();
}

/**
 * @brief compute the norm of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the norm of.
 * @return The norm of the nion.
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline T norm(const nion<T, D> &z) {
    return z.norm();
}

/**
 * @brief compute the inverse of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse of.
 * @return The inverse of the nion.
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> inv(const nion<T, D> &z) {
    return z.inv();
}

/**
 * @brief compute the dot product of two nions.
 * @tparam T type of the nions.
 * @param lhs The left hand side nion.
 * @param rhs The right hand side nion.
 * @return The dot product of the nions.
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline T dot(const nion<T, D> &lhs, const nion<T, D> &rhs) {
    T dotProduct = 0;
    #pragma vector aligned
    for (D i = 0; i < std::min(lhs.degree, rhs.degree); i++) {
        dotProduct += lhs[i] * rhs[i];
    }
    return dotProduct;
}

/**
 * @brief compute the cross product of two nions.
 * @tparam T type of the nions.
 * @param lhs The left hand side nion.
 * @param rhs The right hand side nion.
 * @return The cross product of the nions.
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> cross(const nion<T, D> &lhs, const nion<T, D> &rhs); //TODO: implement cross product

/**
 * @brief compute the wedge product of two nions.
 * @tparam T type of the nions.
 * @param lhs The left hand side nion.
 * @param rhs The right hand side nion.
 * @return The wedge product of the nions.
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> wedge(const nion<T, D> &lhs, const nion<T, D> &rhs); //TODO: implement wedge product

/**
 * @brief compute the outer product of two nions.
 * @tparam T type of the nions.
 * @param lhs The left hand side nion.
 * @param rhs The right hand side nion.
 * @return The outer product of the nions.
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> outer(const nion<T, D> &lhs, const nion<T, D> &rhs); //TODO: implement outer product


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
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> exp(const nion<T, D> &z) noexcept {

    // get polar form of nion
    nion<T, D> v = z.imag();

    // make unit vector
    T v_norm = v.norm();
    if (v_norm <= std::numeric_limits<T>::epsilon())
        return nion<T, D>(std::exp(z[0]), z.degree);

    // compute exponential of nion
    return (v * (std::sin(v_norm) / v_norm) + std::cos(v_norm)) * std::exp(z[0]);
}

/**
 * @brief compute the principle logarithm of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the natural logarithm of.
 * @return The natural logarithm of the nion.
 * @details The natural logarithm of a nion is defined as ln(z) = ln(|z|) + v/|v| * atan(|v|/r).
 *          where a is the real component and v is the imaginary components.
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> log(const nion<T, D> &z) noexcept {

    // get polar form of nion
    T r = z.real();
    nion<T, D> v = z.imag();

    // compute norms
    T z_abs = z.abs();
    T z_norm = std::sqrt(z_abs);
    T v_norm = std::sqrt(z_abs - r*r);
    T theta = std::atan2(v_norm, r);

    if (v_norm <= std::numeric_limits<T>::epsilon())
        return nion<T, D>(std::log(z_norm), z.degree);

    // compute natural logarithm of nion
    return v * (theta / v_norm) + std::log(z_norm);
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
template<typename T, typename D = uint_fast16_t, typename S = T>
static constexpr inline nion<T, D> pow(const nion<T, D> &z, S power) noexcept {
    // get polar form of nion
    T r = z.real();
    nion<T, D> v = z.imag();
    T powa = static_cast<T>(power);

    // compute norms
    T z_abs = z.abs();
    T z_norm = std::sqrt(z_abs);
    T v_norm = std::sqrt(z_abs - r*r);
    T theta = std::atan2(v_norm, r);

    if (v_norm <= std::numeric_limits<T>::epsilon())
        return nion<T, D>(std::pow(r, powa), z.degree);

    // compute power of nion
    return pow(z_norm, powa) * (cos(theta * powa) + v * (sin(theta * powa) / v_norm));
}

/**
 * @brief compute the power of a nion with nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the power of.
 * @param power The power to raise the nion to.
 * @return The power of the nion.
 * @details The power of a nion is defined as z^p = e^(p * ln(z)).
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> pow(const nion<T, D> &z, const nion<T, D> &power) noexcept {
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
template<typename T, typename D = uint_fast16_t, typename S = T>
static constexpr inline nion<T, D> pow(S x, const nion<T, D> &z) noexcept {
    return exp(z * log(static_cast<T>(x)));
}

/**
 * @brief compute the square root of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the square root of.
 * @return The square root of the nion.
 * @details The square root of a nion is defined as sqrt(z) = z^(1/2).
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> sqrt(const nion<T, D> &z) noexcept{
    return pow(z, 0.5l);
}

/**
 * @brief compute the cube root of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the cube root of.
 * @return The cube root of the nion.
 * @details The cube root of a nion is defined as cbrt(z) = z^(1/3).
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> cbrt(const nion<T, D> &z) noexcept{
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
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> sinh(const nion<T, D> &z) noexcept{
    // get polar form of nion
    nion<T, D> v = z.imag();

    // make unit vector
    T v_norm = v.norm();
    if (v_norm <= std::numeric_limits<T>::epsilon())
        return nion<T, D>(std::exp(z[0]), z.degree);

    T e_z = std::exp(z[0]) / static_cast<T>(2);
    T e_mz = std::exp(-z[0]) / static_cast<T>(2);

    // compute exponential of nion
    nion<T, D> sin_nion = v * (( e_z + e_mz) * std::sin(v_norm) / v_norm);
    sin_nion += std::cos(v_norm) * (e_z - e_mz);
    return sin_nion;
}

/**
 * @brief compute the hyperbolic cosine of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic cosine of.
 * @return The hyperbolic cosine of the nion.
 * @details The hyperbolic cosine of a nion is defined as cosh(z) = (e^z + e^-z) / 2.
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> cosh(const nion<T, D> &z) noexcept{
    // get polar form of nion
    nion<T, D> v = z.imag();

    // make unit vector
    T v_norm = v.norm();
    if (v_norm <= std::numeric_limits<T>::epsilon())
        return nion<T, D>(std::exp(z[0]), z.degree);

    T e_z = std::exp(z[0]) / static_cast<T>(2);
    T e_mz = std::exp(-z[0]) / static_cast<T>(2);

    // compute exponential of nion
    nion<T, D> cos_nion = v * (( e_z - e_mz) * std::sin(v_norm) / v_norm);
    cos_nion += std::cos(v_norm) * (e_z + e_mz);
    return cos_nion;
}

/**
 * @brief compute the hyperbolic tangent of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic tangent of.
 * @return The hyperbolic tangent of the nion.
 * @details The hyperbolic tangent of a nion is defined as tanh(z) = sinh(z) / cosh(z).
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> tanh(const nion<T, D> &z) noexcept{
    return sinh(z) / cosh(z);
}

/**
 * @brief compute the hyperbolic cotangent of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic cotangent of.
 * @return The hyperbolic cotangent of the nion.
 * @details The hyperbolic cotangent of a nion is defined as coth(z) = 1 / tanh(z).
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> coth(const nion<T, D> &z) noexcept{
    return tanh(z).inv();
}

/**
 * @brief compute the hyperbolic secant of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic secant of.
 * @return The hyperbolic secant of the nion.
 * @details The hyperbolic secant of a nion is defined as sech(z) = 1 / cosh(z).
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> sech(const nion<T, D> &z) noexcept{
    return cosh(z).inv();
}

/**
 * @brief compute the hyperbolic cosecant of a nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the hyperbolic cosecant of.
 * @return The hyperbolic cosecant of the nion.
 * @details The hyperbolic cosecant of a nion is defined as csch(z) = 1 / sinh(z).
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> csch(const nion<T, D> &z) noexcept{
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
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> sin(const nion<T, D> &z) noexcept{
    // get the polar form of the nion
    T r = real(z);
    nion<T, D> v = imag(z);

    // make unit vector
    T v_norm = v.norm();
    if (v_norm <= std::numeric_limits<T>::epsilon())
        return nion<T, D>(sin(r), z.degree);

    // compute the sine of the nion
    return v * (sinh(v_norm) * cos(r) / v_norm) + sin(r) * cosh(v_norm);
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
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> cos(const nion<T, D> &z) noexcept{
    // get the polar form of the nion
    T r = real(z);
    nion<T, D> v = imag(z);

    // make unit vector
    T v_norm = v.norm();
    if (v_norm <= std::numeric_limits<T>::epsilon())
        return nion<T, D>(cos(r), z.degree);


    // compute the cosine of the nion
    return -v * (sinh(v_norm) * sin(r) / v_norm) + cos(r) * cosh(v_norm);
}

/**
 * @brief compute the tangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the tangent of.
 * @return The tangent of the nion.
 * @details The tangent of the nion is defined as
 * tan(z) = tan(a + bi) = (tan(a) + tanh(b)i) / (1 - tan(a) * tanh(b) i).
 * @see https://en.wikipedia.org/wiki/Proofs_of_trigonometric_identities#Angle_sum_identities
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> tan(const nion<T, D> &z) noexcept{
    // get the polar form of the nion
    T r = real(z);
    nion<T, D> v = imag(z);

    // make unit vector
    T v_norm = v.norm();
    if (v_norm <= std::numeric_limits<T>::epsilon())
        return nion<T, D>(tan(r), z.degree);

    // compute the tangent of the nion
    return (tan(r) + v * (tanh(v_norm) / v_norm)) / (1 - v*(tan(r) / v_norm * tanh(v_norm)));
}

/**
 * @brief compute the cotangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the cotangent of.
 * @return The cotangent of the nion.
 * @details The cotangent of the nion is defined as cot(z) = 1 / tan(z).
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> cot(const nion<T, D> &z) noexcept{
    return tan(z).inv();
}

/**
 * @brief compute the secant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the secant of.
 * @return The secant of the nion.
 * @details The secant of the nion is defined as sec(z) = 1 / cos(z).
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> sec(const nion<T, D> &z) noexcept{
    return cos(z).inv();
}

/**
 * @brief compute the cosecant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the cosecant of.
 * @return The cosecant of the nion.
 * @details The cosecant of the nion is defined as csc(z) = 1 / sin(z).
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> csc(const nion<T, D> &z) noexcept{
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
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> asinh(const nion<T, D> &z) noexcept{
    return log(z + sqrt(1 + pow(z,2)));
}

/**
 * @brief compute the inverse hyperbolic cosine of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse hyperbolic cosine of.
 * @return The inverse hyperbolic cosine of the nion.
 * @details The inverse hyperbolic cosine of the nion is defined as acosh(z) = ln(z + sqrt(z + 1)*sqrt(z - 1)).
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://mathworld.wolfram.com/InverseHyperbolicCosine.html
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> acosh(const nion<T, D> &z) noexcept{
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
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> atanh(const nion<T, D> &z) noexcept{
    return (log(1 + z) - log(1 - z)) * 0.5l;
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
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> acoth(const nion<T, D> &z) noexcept{
    return (log(1 + inv(z)) - log(1 - inv(z))) * 0.5l;
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
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> asech(const nion<T, D> &z) noexcept{
    return log(sqrt(pow(z,-2) - 1) + inv(z));
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
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> acsch(const nion<T, D> &z) noexcept{
    return log(sqrt(1 + pow(z,-2)) + inv(z));
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
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> asin(const nion<T, D> &z) noexcept {
    // get the polar form of the nion
    nion<T, D> v = z.imag();

    // make unit vector
    T v_norm = v.norm();
    if (v_norm == 0)
        return nion<T, D>(asin(z.real()), z.degree);
    v /= v_norm;


    // compute the inv sine of the nion
    return v * log(sqrt(1 - pow(z,2)) - v * z);
}

/**
 * @brief compute the inverse cosine of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse cosine of.
 * @return The inverse cosine of the nion.
 * @details The inverse cosine of the nion is defined as acos(z) = pi/2 - asin(z).
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Extension_to_complex_plane
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> acos(const nion<T, D> &z) noexcept {
    return M_PI_2l - asin(z);
}

/**
 * @brief compute the inverse tangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse tangent of.
 * @return The inverse tangent of the nion.
 * @details The inverse tangent of the nion is defined as atan(z) = v * atanh(-v*z) = -v/2 * log((v-z)/(v+z))
 * @note where r is the real part of z and v is the complex components of z in polar form.
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Logarithmic_forms
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> atan(const nion<T, D> &z) noexcept {
    // get the polar form of the nion
    nion<T, D> v = z.imag();

    // make unit vector
    T v_norm = v.norm();
    if (v_norm == 0)
        return nion<T, D>(atan(z.real()), z.degree);
    v /= v_norm;

    // compute the inv tangent of the nion
    return v * log((v - z)/(v + z)) * -0.5l;
}

/**
 * @brief compute the inverse cotangent of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse cotangent of.
 * @return The inverse cotangent of the nion.
 * @details The inverse cotangent of the nion is defined as acot(z) = pi/2 - atan(z).
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Extension_to_complex_plane
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> acot(const nion<T, D> &z) noexcept {
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
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> asec(const nion<T, D> &z) noexcept {
    return acos(inv(z));
}

/**
 * @brief compute the inverse cosecant of the nion.
 * @tparam T type of the nion.
 * @param z The nion to compute the inverse cosecant of.
 * @return The inverse cosecant of the nion.
 * @details The inverse cosecant of the nion is defined as acsc(z) = asin(1/z).
 * @see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Extension_to_complex_plane
 */
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> acsc(const nion<T, D> &z) noexcept {
    return asin(inv(z));
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
template<typename T, typename D = uint_fast16_t>
static constexpr inline nion<T, D> gamma(const nion<T, D> &z) noexcept {
    // compute the gamma function of the nion
    return exp(-z) * sqrt(inv(z)) * pow(1/(12 * z - inv(10 * z)) + z, z) * sqrt(2 * M_PIl);
}

#endif //NION_NION_HPP
