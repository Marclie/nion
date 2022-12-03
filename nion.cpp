#include <cmath>
#include <limits>
#include <type_traits>
#include <iostream>
#include <stdexcept>
#include "nion.hpp"

using Nion::nion;
namespace Nion {

    /***************************
    *  NION CONSTRUCTORS
    ***************************/

    template<typename T, typename D>
    constexpr size_t Nion::nion<T, D>::alignSize(D size) {
        // checks if degree is a power of 2
        if (size & (size - 1)) {
            return size;
        }

        // if degree is not a power of 2, align it to the next power of 2
        size_t alignedDegree = 1;
        while (alignedDegree <= size) {
            alignedDegree <<= 1;
        }
        return alignedDegree;
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(T* vals, D degree) : degree(degree) {

        // check if the degree is greater than zero
        if (degree <= 0) {
            throw std::invalid_argument("The degree of the nion must be greater than zero.");
        }
        this->components = (T *) malloc(degree * sizeof(T));
        memcpy(this->components, vals, degree * sizeof(T));
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(const std::initializer_list<T> &vals) : degree(vals.size()) {
        if (degree <= 0) {
            throw std::invalid_argument("The degree of the nion must be greater than zero.");
        }
        this->components = (T *) malloc(degree * sizeof(T));
        memcpy(this->components, vals.begin(), degree * sizeof(T));
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(D degree) : degree(degree) {
        // check if the degree is greater than zero
        if (degree <= 0) {
            throw std::invalid_argument("The degree of the nions must be greater than zero.");
        }
        this->components = (T *) calloc(this->degree, sizeof(T));
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(const nion<T, D> &other) : degree(other.degree) {
        this->components = (T *) malloc(other.degree * sizeof(T));
        memcpy(this->components, other.components, other.degree * sizeof(T));
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>(nion<T, D> &&other) noexcept{
        this->degree = other.degree;
        this->components = other.components;
        other.components = nullptr;
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline nion<T, D>::nion(S realVal, D degree) : degree(degree) {
        // check if the degree is greater than zero
        if (degree <= 0) {
            throw std::invalid_argument("The degree of the nions must be greater than zero.");
        }
        this->components = (T *) malloc(degree * sizeof(T));
        components[0] = static_cast<T>(realVal);
        if (degree > 1) {
            memset(components + 1, 0, (degree - 1) * sizeof(T));
        }
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(const nion<T, D> &a, const nion<T, D> &b) : degree(a.degree + b.degree) {
        this->components = (T *) malloc(degree * sizeof(T));

        memcpy(this->components, a.components, a.degree * sizeof(T));
        memcpy(this->components + a.degree, b.components, b.degree * sizeof(T));
    }

    /************************************
    *         ASSIGNMENT OPERATORS
    *************************************/
    
    template<typename T, typename D>
    constexpr inline nion<T, D> &nion<T, D>::operator=(const nion<T, D> &other) noexcept {
        if (&other == this) {
            return *this;
        }
        if (this->degree != other.degree) {
            this->degree = other.degree;
            if (this->components == nullptr) {
                this->components = (T *) malloc(other.degree * sizeof(T));
            } else {
                this->components = (T *) realloc(this->components, other.degree * sizeof(T));
            }
        }

        memcpy(this->components, other.components, degree * sizeof(T));
        return *this;
    }

    template<typename T, typename D>
    constexpr inline nion<T, D> &nion<T, D>::operator=(nion<T, D> &&other) noexcept {
        if (&other != this) {
            free(this->components);
            this->degree = other.degree;
            this->components = other.components;
            other.components = nullptr;
        }
        return *this;
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline nion<T, D> &nion<T, D>::operator=(S scalar) {
        free(this->components);
        this->components = nullptr;
        *this = nion<T, D>(static_cast<T>(scalar), 1);
        return *this;
    }

    /************************************
    *  ASSIGNMENT AND ADDITION OPERATORS
    *************************************/

    template<typename T, typename D>
    constexpr inline void nion<T, D>::resize(D newDegree) {
        components = (T *) realloc(components, newDegree * sizeof(T));
        if (newDegree > degree)
            memset(components + degree, 0, (newDegree - degree) * sizeof(T));
        degree = newDegree;
    }
    
    template<typename T, typename D>
    constexpr inline void nion<T, D>::operator+=(const nion<T, D> &other) {
        // if the degree is less than the other nion, resize this nion.
        if (degree < other.degree)
            resize(other.degree);

        // add the components of the other nion to this nion.
        #pragma vector aligned
        for (D i = 0; i < other.degree; i++) {
            components[i] += other[i];
        }
    }

    template<typename T, typename D>
    constexpr inline void nion<T, D>::operator-=(const nion<T, D> &other) {
        // if the degree is less than the other nion, resize this nion.
        if (this->degree < other.degree)
            resize(other.degree);

        // substract the components of the other nion from this nion.
        #pragma vector aligned
        for (D i = 0; i < other.degree; i++) {
            components[i] -= other[i];
        }
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline void nion<T, D>::operator+=(S scalar) const {
        components[0] += static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline void nion<T, D>::operator-=(S scalar) const {
        components[0] -= static_cast<T>(scalar);
    }

    /************************************
    *        ADDITION OPERATORS
    *************************************/

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator+(const nion<T, D> &other) const {
        nion<T, D> sum = *this;
        sum += other;
        return sum;
    }

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator-(const nion<T, D> &other) const {
        nion<T, D> difference = *this;
        difference -= other;
        return difference;
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline nion<T, D> nion<T, D>::operator+(S scalar) const {
        nion<T, D> sum = *this;
        sum.components[0] += static_cast<T>(scalar);
        return sum;
    }

    template<typename T, typename D, typename S>
    static constexpr inline nion<T, D> operator+(S scalar, const nion<T, D> &z) {
        return z + static_cast<T>(scalar);
    }
    
    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator++() {
        components[0]++;
        return *this;
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline nion<T, D> nion<T, D>::operator-(S scalar) const {
        return *this + (static_cast<T>(-scalar));
    }

    template<typename T, typename D, typename S>
    static constexpr inline nion<T, D> operator-(S scalar, const nion<T, D> &z) {
        return -z + static_cast<T>(scalar);
    }
    
    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator--() {
        components[0]--;
        return *this;
    }
    
    /******************************************
    *  ASSIGNMENT AND MULTIPLICATION OPERATORS
    *******************************************/

    template<typename T, typename D>
    constexpr inline void nion<T, D>::operator*=(const nion<T, D> &other) {
        *this = *this * other;
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline void nion<T, D>::operator*=(S scalar) {
        #pragma vector aligned
        for (D i = 0; i < degree; i++) {
            components[i] *= static_cast<T>(scalar);
        }
    }
    
    template<typename T, typename D>
    constexpr inline void nion<T, D>::operator/=(const nion<T, D> &other) {
        *this = *this / other;
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline void nion<T, D>::operator/=(S scalar) {
        # pragma vector aligned
        for (D i = 0; i < degree; i++) {
            components[i] /= static_cast<T>(scalar);
        }
    }

    /******************************************
    *        MULTIPLICATION OPERATORS
    *******************************************/

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::conj() const {
        nion<T, D> conjugate = *this;

        // negate all components except the first
        #pragma vector aligned
        for (D i = 1; i < degree; i++) {
            conjugate.components[i] = -conjugate.components[i];
        }

        return conjugate;
    }

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator-() const {
        nion<T, D> negated = *this;

        #pragma vector aligned
        for (D i = 0; i < degree; i++) {
            negated.components[i] = -negated.components[i];
        }

        return negated;
    }
    
    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator*(const nion<T, D> &other) const {

        switch (degree) {
            case 1:
                // if the degree is 1, then the product is just the scalar product.
                if (other.degree == 1) {
                    nion<T, D> product;
                    product.degree = 1;
                    product.components = (T *) aligned_alloc(1, sizeof(T));
                    product.components[0] = components[0] * other.components[0];
                    return std::move(product);
                }
                else
                    return other * components[0];
            case 2:
                if (other.degree == 2) { // hard-coded complex product
                    nion<T, D> product;
                    product.degree = 2;
                    product.components = (T *) aligned_alloc(2, sizeof(T) * 2);
                    product.components[0] = components[0] * other[0] - components[1] * other[1];
                    product.components[1] = components[1] * other[0] + components[0] * other[1];
                    return product;
                }
                break;
            case 4:
                if (other.degree == 4) { // hard-coded quaternion product
                    nion<T, D> product;
                    product.degree = 4;
                    product.components = (T *) aligned_alloc(4, sizeof(T) * 4);
                    product.components[0] = components[0] * other[0] - components[1] * other[1] - components[2] * other[2] - components[3] * other[3];
                    product.components[1] = components[1] * other[0] + components[0] * other[1] - components[3] * other[2] + components[2] * other[3];
                    product.components[2] = components[2] * other[0] + components[3] * other[1] + components[0] * other[2] - components[1] * other[3];
                    product.components[3] = components[3] * other[0] - components[2] * other[1] + components[1] * other[2] + components[0] * other[3];
                    return product;
                }
                break;
            case 8:
                if (other.degree == 8) { // hard-coded octonion product
                    nion<T, D> product;
                    product.degree = 8;
                    product.components = (T *) aligned_alloc(8, sizeof(T) * 8);
                    product.components[0] = components[0] * other[0] - components[1] * other[1] - components[2] * other[2] - components[3] * other[3] - components[4] * other[4] - components[5] * other[5] - components[6] * other[6] - components[7] * other[7];
                    product.components[1] = components[1] * other[0] + components[0] * other[1] - components[3] * other[2] + components[2] * other[3] - components[5] * other[4] + components[4] * other[5] + components[7] * other[6] - components[6] * other[7];
                    product.components[2] = components[2] * other[0] + components[3] * other[1] + components[0] * other[2] - components[1] * other[3] - components[6] * other[4] - components[7] * other[5] + components[4] * other[6] + components[5] * other[7];
                    product.components[3] = components[3] * other[0] - components[2] * other[1] + components[1] * other[2] + components[0] * other[3] - components[7] * other[4] + components[6] * other[5] - components[5] * other[6] + components[4] * other[7];
                    product.components[4] = components[4] * other[0] + components[5] * other[1] + components[6] * other[2] + components[7] * other[3] + components[0] * other[4] - components[1] * other[5] - components[2] * other[6] - components[3] * other[7];
                    product.components[5] = components[5] * other[0] - components[4] * other[1] + components[7] * other[2] - components[6] * other[3] + components[1] * other[4] + components[0] * other[5] + components[3] * other[6] - components[2] * other[7];
                    product.components[6] = components[6] * other[0] - components[7] * other[1] - components[4] * other[2] + components[5] * other[3] + components[2] * other[4] - components[3] * other[5] + components[0] * other[6] + components[1] * other[7];
                    product.components[7] = components[7] * other[0] + components[6] * other[1] - components[5] * other[2] - components[4] * other[3] + components[3] * other[4] + components[2] * other[5] - components[1] * other[6] + components[0] * other[7];
                    return product;
                }
                break;
            case 16:
                if (other.degree == 16) { // hard-coded sedenion product
                    nion<T, D> product;
                    product.degree = 16;
                    product.components = (T *) aligned_alloc(16, sizeof(T) * 16);

                    product.components[0] = components[0] *  other[0] - components[1] *  other[1] - components[2] *  other[2] - components[3] *  other[3] - components[4] *  other[4] - components[5] *  other[5] - components[6] *  other[6] - components[7] *  other[7] - components[8] *  other[8] - components[9] *  other[9] - components[10] * other[10] - components[11] * other[11] - components[12] * other[12] - components[13] * other[13] - components[14] * other[14] - components[15] * other[15];
                    product.components[1] = components[1] *  other[0] + components[0] *  other[1] - components[3] *  other[2] + components[2] *  other[3] - components[5] *  other[4] + components[4] *  other[5] - components[7] *  other[6] + components[6] *  other[7] - components[9] *  other[8] + components[8] *  other[9] - components[11] * other[10] + components[10] * other[11] - components[13] * other[12] + components[12] * other[13] - components[15] * other[14] + components[14] * other[15];
                    product.components[2] = components[0] *  other[2] - components[1] *  other[3] + components[2] *  other[0] + components[3] *  other[1] + components[4] *  other[6] + components[5] *  other[7] - components[6] *  other[4] - components[7] *  other[5] + components[8] * other[10] + components[9] * other[11] - components[10] *  other[8] - components[11] *  other[9] - components[12] * other[14] - components[13] * other[15] + components[14] * other[12] + components[15] * other[13];
                    product.components[3] = components[0] *  other[3] + components[1] *  other[2] - components[2] *  other[1] + components[3] *  other[0] + components[4] *  other[7] - components[5] *  other[6] + components[6] *  other[5] - components[7] *  other[4] + components[8] * other[11] - components[9] * other[10] + components[10] *  other[9] - components[11] *  other[8] - components[12] * other[15] + components[13] * other[14] - components[14] * other[13] + components[15] * other[12];
                    product.components[4] = components[0] *  other[4] - components[1] *  other[5] - components[2] *  other[6] - components[3] *  other[7] + components[4] *  other[0] + components[5] *  other[1] + components[6] *  other[2] + components[7] *  other[3] + components[8] * other[12] + components[9] * other[13] + components[10] * other[14] + components[11] * other[15] - components[12] *  other[8] - components[13] *  other[9] - components[14] * other[10] - components[15] * other[11];
                    product.components[5] = components[0] *  other[5] + components[1] *  other[4] - components[2] *  other[7] + components[3] *  other[6] - components[4] *  other[1] + components[5] *  other[0] - components[6] *  other[3] + components[7] *  other[2] + components[8] * other[13] - components[9] * other[12] + components[10] * other[15] - components[11] * other[14] + components[12] *  other[9] - components[13] *  other[8] + components[14] * other[11] - components[15] * other[10];
                    product.components[6] = components[0] *  other[6] + components[1] *  other[7] + components[2] *  other[4] - components[3] *  other[5] - components[4] *  other[2] + components[5] *  other[3] + components[6] *  other[0] - components[7] *  other[1] + components[8] * other[14] - components[9] * other[15] - components[10] * other[12] + components[11] * other[13] + components[12] * other[10] - components[13] * other[11] - components[14] *  other[8] + components[15] * other[9];
                    product.components[7] = components[0] *  other[7] - components[1] *  other[6] + components[2] *  other[5] + components[3] *  other[4] - components[4] *  other[3] - components[5] *  other[2] + components[6] *  other[1] + components[7] *  other[0] + components[8] * other[15] + components[9] * other[14] - components[10] * other[13] - components[11] * other[12] + components[12] * other[11] + components[13] * other[10] - components[14] *  other[9] - components[15] * other[8];
                    product.components[8] = components[0] *  other[8] - components[1] *  other[9] - components[2] * other[10] - components[3] * other[11] - components[4] * other[12] - components[5] * other[13] - components[6] * other[14] - components[7] * other[15] + components[8] *  other[0] + components[9] *  other[1] + components[10] *  other[2] + components[11] *  other[3] + components[12] *  other[4] + components[13] *  other[5] + components[14] *  other[6] + components[15] * other[7];
                    product.components[9] = components[0] *  other[9] + components[1] *  other[8] - components[2] * other[11] + components[3] * other[10] - components[4] * other[13] + components[5] * other[12] + components[6] * other[15] - components[7] * other[14] - components[8] *  other[1] + components[9] *  other[0] - components[10] *  other[3] + components[11] *  other[2] - components[12] *  other[5] + components[13] *  other[4] + components[14] *  other[7] - components[15] * other[6];
                    product.components[10] = components[0] * other[10] + components[1] * other[11] + components[2] *  other[8] - components[3] *  other[9] - components[4] * other[14] - components[5] * other[15] + components[6] * other[12] + components[7] * other[13] - components[8] *  other[2] + components[9] *  other[3] + components[10] *  other[0] - components[11] *  other[1] - components[12] *  other[6] - components[13] *  other[7] + components[14] *  other[4] + components[15] * other[5];
                    product.components[11] = components[0] * other[11] - components[1] * other[10] + components[2] *  other[9] + components[3] *  other[8] - components[4] * other[15] + components[5] * other[14] - components[6] * other[13] + components[7] * other[12] - components[8] *  other[3] - components[9] *  other[2] + components[10] *  other[1] + components[11] *  other[0] - components[12] *  other[7] + components[13] *  other[6] - components[14] *  other[5] + components[15] * other[4];
                    product.components[12] = components[0] * other[12] + components[1] * other[13] + components[2] * other[14] + components[3] * other[15] + components[4] *  other[8] - components[5] *  other[9] - components[6] * other[10] - components[7] * other[11] - components[8] *  other[4] + components[9] *  other[5] + components[10] *  other[6] + components[11] *  other[7] + components[12] *  other[0] - components[13] *  other[1] - components[14] *  other[2] - components[15] * other[3];
                    product.components[13] = components[0] * other[13] - components[1] * other[12] + components[2] * other[15] - components[3] * other[14] + components[4] *  other[9] + components[5] *  other[8] + components[6] * other[11] - components[7] * other[10] - components[8] *  other[5] - components[9] *  other[4] + components[10] *  other[7] - components[11] *  other[6] + components[12] *  other[1] + components[13] *  other[0] + components[14] *  other[3] - components[15] * other[2];
                    product.components[14] = components[0] * other[14] - components[1] * other[15] - components[2] * other[12] + components[3] * other[13] + components[4] * other[10] - components[5] * other[11] + components[6] *  other[8] + components[7] *  other[9] - components[8] *  other[6] - components[9] *  other[7] - components[10] *  other[4] + components[11] *  other[5] + components[12] *  other[2] - components[13] *  other[3] + components[14] *  other[0] + components[15] * other[1];
                    product.components[15] = components[0] * other[15] + components[1] * other[14] - components[2] * other[13] - components[3] * other[12] + components[4] * other[11] + components[5] * other[10] - components[6] *  other[9] + components[7] *  other[8] - components[8] *  other[7] + components[9] *  other[6] - components[10] *  other[5] - components[11] *  other[4] + components[12] *  other[3] + components[13] *  other[2] - components[14] *  other[1] + components[15] * other[0];

                    return product;
                }
                break;
            //case 32: TODO: No way. You can't make me. 32x32 products are just too big. I could do it, but it wouldn't be pretty. I'll leave it to you.
        }

        // if the other degree is 1, then the product is just the scalar product.
        if (other.degree == 1)
            return *this * other[0];


        ///*** if the degree is not 1, 2, 4, or 8, then we need to do the recursive product ***

        // the degree of the first half of the nion
        D this_half_degree = degree - degree / 2;
        D other_half_degree = other.degree - other.degree / 2;

        nion<T, D> a, b, c, d;

        a.degree = this_half_degree;
        b.degree = degree - this_half_degree;
        c.degree = other_half_degree;
        d.degree = other.degree - other_half_degree;

        a.components = this->components;
        b.components = this->components + this_half_degree;
        c.components = other.components;
        d.components = other.components + other_half_degree;

        // calculate the product
        nion<T, D> product(
                (a * c) - (d.conj() * b), // consider adding involution parameter for sign
                (d * a) + (b * c.conj())
        );

        a.components = nullptr;
        b.components = nullptr;
        c.components = nullptr;
        d.components = nullptr;
        return product;

    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline nion<T, D> nion<T, D>::operator*(S scalar) const {
        nion<T, D> product;
        product.degree = degree;
        product.components = (T *) malloc(degree * sizeof(T));

        #pragma vector aligned
        for (D i = 0; i < degree; i++) {
            product.components[i] = components[i] * static_cast<T>(scalar);
        }
        return product;
    }
    
    template<typename T, typename D, typename S>
    static constexpr inline nion<T, D> operator*(S scalar, const nion<T, D> &z) {
        static_assert(std::is_arithmetic_v<S>, "scalar must be arithmetic");
        return z * static_cast<T>(scalar);
    }
    
    template<typename T, typename D>
    constexpr inline T nion<T, D>::abs() const {
        T absVal = 0;

        #pragma vector aligned
        for (D i = 0; i < degree; i++) {
            absVal += components[i] * components[i];
        }

        return absVal;
    }

    template<typename T, typename D>
    constexpr inline T nion<T, D>::norm() const {
        return std::sqrt(abs());
    }
    
    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::inv() const {
        constexpr T epsilon = std::numeric_limits<T>::epsilon();

        nion<T, D> inverse = *this;
        T absolute = abs();
        if (absolute < epsilon) {
            // if the absolute value is zero, then use the product definition of the absolute value.
            // zero divisors are possible in nions with degree >= 16, so we need to check for them.
            absolute = (*this * this->conj()).real();
        }

        inverse[0] /= absolute;
        #pragma vector aligned
        for (D i = 1; i < degree; i++) {
            inverse.components[i] /= -absolute;
        }
        return inverse;
    }

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator/(const nion<T, D> &other) const {
        return *this * other.inv();
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline nion<T, D> nion<T, D>::operator/(S scalar) const {
        return *this * (static_cast<T>(1 / scalar));
    }

    template<typename T, typename D, typename S>
    static constexpr inline nion<T, D> operator/(S scalar, const nion<T, D> &z) {
        return z.inv() * static_cast<T>(scalar);
    }


    /******************************************
    *            ACCESSOR FUNCTIONS
    *******************************************/

    template<typename T, typename D>
    constexpr inline T &nion<T, D>::operator[](D index) const {
        return this->components[index];
    }
    
    template<typename T, typename D>
    constexpr inline T nion<T, D>::real() const { return components[0]; }

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::imag() const {
        nion<T, D> imag;

        T *imag_components = (T *) malloc(degree * sizeof(T));
        imag.components = imag_components;
        imag.degree = degree;

        imag[0] = 0;
        if (degree > 1)
            memcpy(imag_components + 1, components + 1, (degree - 1) * sizeof(T));

        return imag;
    }

    /******************************************
    *            COMPARISON OPERATORS
    *******************************************/

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::operator==(const nion<T, D> &other) const {
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
    }

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::operator!=(const nion<T, D> &other) const {
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
    }

    template<typename T, typename D>
    constexpr inline T nion<T, D>::proj() const { return copysign(norm(), real()); }
    
    template<typename T, typename D>
    constexpr inline bool nion<T, D>::operator>(const nion<T, D> &other) const {
        return proj() > other.proj();
    }

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::operator<(const nion<T, D> &other) const {
        return proj() < other.proj();
    }

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::operator>=(const nion<T, D> &other) const {
        if (proj() > other.proj())
            return true;
        return *this == *other;
    }

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::operator<=(const nion<T, D> &other) const {
        if (proj() < other.proj())
            return true;
        return *this == *other;
    }

    template<typename T, typename D, typename S>
    constexpr inline bool value_is_similar(const T a, const S b, const T epsilon){
        return std::fabs(a - static_cast<T>(b)) <= epsilon;
    }

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::is_real() const {
        #pragma vector aligned
        for (D i = 1; i < degree; i++) {
            if (!value_is_similar(components[i], 0)) {
                return false;
            }
        }
        return true;
    }
    
    template<typename T, typename D>
    template<typename S>
    constexpr inline bool nion<T, D>::operator==(S scalar) const {
        if (!value_is_similar(real(), static_cast<T>(scalar))) {
            return false;
        }
        return is_real();
    }

    template<typename T, typename D, typename S>
    static constexpr inline bool operator==(S scalar, const nion<T, D> &z) {
        return z == static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline bool nion<T, D>::operator!=(S scalar) const {
        if (!value_is_similar(real(), static_cast<T>(scalar))) {
            return true;
        }
        return !is_real();
    }

    template<typename T, typename D, typename S>
    static constexpr inline bool operator!=(S scalar, const nion<T, D> &z) {
        return z != static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline bool nion<T,D>::operator>(S scalar) const{
        return proj() > static_cast<T>(scalar);
    }

    template<typename T, typename D, typename S>
    static constexpr inline bool operator>(S scalar, const nion<T, D> &z) {
        return z < static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline bool nion<T,D>::operator<(S scalar) const{
        return proj() < static_cast<T>(scalar);
    }

    template<typename T, typename D, typename S>
    static constexpr inline bool operator<(S scalar, const nion<T, D> &z) {
        return z > static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline bool nion<T,D>::operator>=(S scalar) const{
        if (*this == nion<T, D>(static_cast<T>(scalar), this->degree))
            return true;
        return proj() > static_cast<T>(scalar);
    }

    template<typename T, typename D, typename S>
    static constexpr inline bool operator>=(S scalar, const nion<T, D> &z) {
        return z <= static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline bool nion<T,D>::operator<=(S scalar) const{
        if (*this == nion<T, D>(static_cast<T>(scalar), this->degree))
            return true;
        return proj() < static_cast<T>(scalar);
    }
    
    template<typename T, typename D, typename S>
    static constexpr inline bool operator<=(S scalar, const nion<T, D> &z) {
        return z >= static_cast<T>(scalar);
    }

    /******************************************
    *            STREAMING OPERATORS
    *******************************************/
    
    template<typename T, typename D>
    static std::ostream &operator<<(std::ostream &os, const nion<T, D> &z) {
        T component = z.components[0];
        os << "(" << component;

        for (D i = 1; i < z.degree; i++) {
            component = z.components[i];
            os << "," << component;
        }
        os << ")";
        return os;
    }

    template<typename T, typename D>
    static std::istream &operator>>(std::istream &is, nion<T, D> &z) {
        for (D i = 0; i < z.degree; i++) {
            is >> z.components[i];
        }
        return is;
    }

    /**
    * @brief Converts a nion to a string.
    * @return The string representation of the nion.
    */
    template<typename T, typename D>
    inline std::string nion<T, D>::to_string() const {
        std::string nion_string = "(" + std::to_string(components[0]);
        for (D i = 1; i < degree; i++) {
            nion_string += "," + std::to_string(components[i]);
        }
        nion_string += ")";
        return nion_string;
    }

/*********************************
*  NION FUNCTION IMPLEMENTATIONS *
**********************************/

    template<typename T, typename D>
    static constexpr inline T real(const nion<T, D> &z) {
        return z.real();
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> imag(const nion<T, D> &z) {
        return z.imag();
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> conj(const nion<T, D> &z) {
        return z.conj();
    }

    template<typename T, typename D>
    static constexpr inline T abs(const nion<T, D> &z) {
        return z.abs();
    }

    template<typename T, typename D>
    static constexpr inline T norm(const nion<T, D> &z) {
        return z.norm();
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> inv(const nion<T, D> &z) {
        return z.inv();
    }

    template<typename T, typename D>
    static constexpr inline T dot(const nion<T, D> &lhs, const nion<T, D> &rhs) {
        T dotProduct = 0;
        #pragma vector aligned
        for (D i = 0; i < std::min(lhs.degree, rhs.degree); i++) {
            dotProduct += lhs[i] * rhs[i];
        }
        return dotProduct;
    }

    /*
    template<typename T, typename D>
    static constexpr inline nion<T, D>
    cross(const nion<T, D> &lhs, const nion<T, D> &rhs){} //TODO: implement cross product

    template<typename T, typename D>
    static constexpr inline nion<T, D>
    wedge(const nion<T, D> &lhs, const nion<T, D> &rhs){} //TODO: implement wedge product

    template<typename T, typename D>
    static constexpr inline nion<T, D>
    outer(const nion<T, D> &lhs, const nion<T, D> &rhs){} //TODO: implement outer product
     */


    /****************************
    *  NION ALGEBRAIC FUNCTIONS *
    *****************************/

    template<typename T, typename D>
    static constexpr inline nion<T, D> exp(const nion<T, D> &z) noexcept {

        // get polar form of nion
        T r = z.real();
        nion<T, D> i = z.imag();

        // make unit vector
        T i_norm = i.norm();

        // compute exponential of nion
        T cos_theta;
        T sin_theta;
        T exp_r = std::exp(r);

        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (std::abs(i_norm) < epsilon)
            return nion<T, D>(exp_r, z.degree);

        // compute exponential of nion
        cos_theta = std::cos(i_norm);
        sin_theta = std::sin(i_norm);

        return i*(exp_r * sin_theta / i_norm) + exp_r * cos_theta;

    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> log(const nion<T, D> &z) noexcept {

        // get polar form of nion
        T r = z.real();
        nion<T, D> i = z.imag();

        // compute norms
        T z_abs = z.abs();
        T z_norm = std::sqrt(z_abs);
        T i_norm = std::sqrt(z_abs - r * r);
        T theta = std::atan2(i_norm, r);

        // compute natural logarithm of nion
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_norm <= epsilon) // TODO: compare to abs not norm for all
            return std::log(z_norm) + i * theta;
        else
            return i * (theta / i_norm) + std::log(z_norm);
    }

    template<typename T, typename D, typename S>
    static constexpr inline nion<T, D> pow(const nion<T, D> &base, S power) noexcept {

        // get polar form of nion
        T r = base.real();
        nion<T, D> i = base.imag();

        // compute norms
        T base_abs = base.abs();
        T base_norm = std::sqrt(base_abs);
        T i_norm = std::sqrt(base_abs - r * r);

        T power_t;
        T theta;

        T cos_ptheta;
        T sin_ptheta;

        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_norm <= epsilon) { // real nion
            if (std::signbit(r)) { // if base is negative
                power_t = static_cast<T>(power);
                theta = power_t * std::atan2(i_norm, r);
                cos_ptheta = std::cos(theta);
                sin_ptheta = std::sin(theta);

                return std::pow(base_norm, power_t) * (cos_ptheta + i * (sin_ptheta));
            } else {
                return nion<T, D>(std::pow(base_norm, static_cast<T>(power)), base.degree); // if base is positive
            }
        }
        
        if constexpr (std::is_integral_v<S>) { // if power is integer, use faster algorithm

            nion<T, D> z = base;
            if (std::signbit(power)) {
                z = inv(z);
                power = -power;
            }

            switch (power) {
                case 0:
                    return nion<T, D>(1, z.degree);
                case 1:
                    return z;
                case 2:
                    return sqr(z);
                default:
                    power_t = static_cast<T>(power);
                    theta = power_t * std::atan2(i_norm, r);
                    cos_ptheta = std::cos(theta);
                    sin_ptheta = std::sin(theta);
                    break;
            }
        } else {
            power_t = static_cast<T>(power);
            theta = power_t * std::atan2(i_norm, r);
            cos_ptheta = std::cos(theta); //cos(atan(y/x)) = 1/sqrt(1+y^2/x^2) --> cos(p*atan(y/x)) = ???
            sin_ptheta = std::sin(theta); //sin(atan(y/x)) = y/(x*sqrt(1 + y^2/x^2)) --> sin(p*atan(y/x)) = ???
        }

        // compute power of nion
        return std::pow(base_norm, power_t) * (cos_ptheta + i * (sin_ptheta / i_norm));
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> pow(const nion<T, D> &base, const nion<T, D> &power) noexcept {
        return exp(power * log(base));
    }

    template<typename T, typename D, typename S>
    static constexpr inline nion<T, D> pow(S base, const nion<T, D> &power) noexcept {
        return pow(nion<T, D>(static_cast<T>(base), power.degree), power);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> sqr(const nion<T, D> &base) noexcept {
        // get polar form of nion
        T r = base.real();
        nion<T, D> i = base.imag();

        // compute norms
        T base_abs = base.abs();
        T base_norm = std::sqrt(base_abs);
        T i_norm = std::sqrt(base_abs - r * r);

        T power_t = static_cast<T>(2);
//        T theta = 2 * std::atan2(i_norm, r);

        T x2 = r * r;
        T y2 = i_norm * i_norm;
        T denom = static_cast<T>(1) / (x2 + y2);
        T cos_2theta = (x2 - y2) * denom;
        T sin_2theta = static_cast<T>(2) * r * i_norm * denom;

        // compute power of nion
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_norm <= epsilon)
            return std::pow(base_norm, power_t) * (cos_2theta + i * sin_2theta);
        else
            return std::pow(base_norm, power_t) * (cos_2theta + i * (sin_2theta / i_norm));
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> sqrt(const nion<T, D> &z) noexcept {
        return pow(z, static_cast<T>(0.5));
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> cbrt(const nion<T, D> &z) noexcept {
        return pow(z, static_cast<T>(1.0) / static_cast<T>(3.0));
    }

    /*******************************************
    *  NION HYPERBOLIC TRIGONOMETRIC FUNCTIONS *
    ********************************************/

    template<typename T, typename D>
    static constexpr inline nion<T, D> sinh(const nion<T, D> &z) noexcept {
        // get polar form of nion
        nion<T, D> i = z.imag();

        // make unit vector
        T i_norm = i.norm();

        // calculate scalars
        T e_z = std::exp(z[0]) / static_cast<T>(2);
        T e_mz = std::exp(-z[0]) / static_cast<T>(2);

        // compute exponential of nion
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_norm <= epsilon) {
            nion<T, D> sin_nion = i * ((e_z + e_mz) * std::sin(i_norm));
            sin_nion += std::cos(i_norm) * (e_z - e_mz);
            return sin_nion;
        } else {
            nion<T, D> sin_nion = i * ((e_z + e_mz) * std::sin(i_norm) / i_norm);
            sin_nion += std::cos(i_norm) * (e_z - e_mz);
            return sin_nion;
        }

    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> cosh(const nion<T, D> &z) noexcept {
        // get polar form of nion
        nion<T, D> i = z.imag();

        // make unit vector
        T i_norm = i.norm();

        // calculate scalars
        T e_z = std::exp(z[0]) / static_cast<T>(2);
        T e_mz = std::exp(-z[0]) / static_cast<T>(2);

        // compute exponential of nion
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_norm <= epsilon) {
            nion<T, D> cos_nion = i * ((e_z - e_mz) * std::sin(i_norm));
            cos_nion += std::cos(i_norm) * (e_z + e_mz);
            return cos_nion;
        } else {
            nion<T, D> cos_nion = i * ((e_z - e_mz) * std::sin(i_norm) / i_norm);
            cos_nion += std::cos(i_norm) * (e_z + e_mz);
            return cos_nion;
        }
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> tanh(const nion<T, D> &z) noexcept {
        return (exp(z*2) - 1) / (exp(z*2) + 1);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> coth(const nion<T, D> &z) noexcept {
        return tanh(z).inv();
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> sech(const nion<T, D> &z) noexcept {
        return cosh(z).inv();
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> csch(const nion<T, D> &z) noexcept {
        return sinh(z).inv();
    }

    /********************************
    *  NION TRIGONOMETRIC FUNCTIONS *
    *********************************/


    template<typename T, typename D>
    static constexpr inline nion<T, D> sin(const nion<T, D> &z) noexcept {
        // get the polar form of the nion
        T r = real(z);
        nion<T, D> i = imag(z);

        // make unit vector
        T i_norm = i.norm();

        // compute the sine of the nion
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_norm <= epsilon)
            return i * (std::sinh(i_norm) * std::cos(r)) + std::sin(r) * std::cosh(i_norm);
        else
            return i * (std::sinh(i_norm) * std::cos(r) / i_norm) + std::sin(r) * std::cosh(i_norm);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> cos(const nion<T, D> &z) noexcept {
        // get the polar form of the nion
        T r = real(z);
        nion<T, D> i = imag(z);

        // make unit vector
        T i_norm = i.norm();

        // compute the cosine of the nion
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_norm <= epsilon)
            return -i * (std::sinh(i_norm) * std::sin(r)) + std::cos(r) * std::cosh(i_norm);
        else
            return -i * (std::sinh(i_norm) * std::sin(r) / i_norm) + std::cos(r) * std::cosh(i_norm);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> tan(const nion<T, D> &z) noexcept {
        // get the polar form of the nion
        T r = real(z);
        nion<T, D> i = imag(z);

        // make unit vector
        T i_norm = i.norm();

        // compute the tangent of the nion
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_norm <= epsilon)
            return (std::tan(r) + i * std::tanh(i_norm)) / (static_cast<T>(1) - i * (std::tan(r) * std::tanh(i_norm)));
        else
            return (std::tan(r) + i * (std::tanh(i_norm) / i_norm)) / (static_cast<T>(1) - i * (std::tan(r) / i_norm * std::tanh(i_norm)));
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> cot(const nion<T, D> &z) noexcept {
        return tan(z).inv();
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> sec(const nion<T, D> &z) noexcept {
        return cos(z).inv();
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> csc(const nion<T, D> &z) noexcept {
        return sin(z).inv();
    }

    /***************************************************
    *  NION INVERSE HYPERBOLIC TRIGONOMETRIC FUNCTIONS *
    ****************************************************/

    template<typename T, typename D>
    static constexpr inline nion<T, D> asinh(const nion<T, D> &z) noexcept {
        return log(z + sqrt(sqr(z) + static_cast<T>(1)));
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> acosh(const nion<T, D> &z) noexcept {
        // compute the inverse hyperbolic cosine of the nion
        return 2 * log(sqrt((z-1) * static_cast<T>(0.5)) + sqrt((z+1) * static_cast<T>(0.5)));

    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> atanh(const nion<T, D> &z) noexcept {
        return (log(static_cast<T>(1) + z) - log(static_cast<T>(1) - z)) * static_cast<T>(0.5);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> acoth(const nion<T, D> &z) noexcept {
        return (log(static_cast<T>(1) + inv(z)) - log(static_cast<T>(1) - inv(z))) * static_cast<T>(0.5);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> asech(const nion<T, D> &z) noexcept {
        return log(sqrt(1/sqr(z) - static_cast<T>(1)) + inv(z));
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> acsch(const nion<T, D> &z) noexcept {
        return log(sqrt(static_cast<T>(1) + 1/sqr(z)) + inv(z));
    }

    /****************************************
    *  NION INVERSE TRIGONOMETRIC FUNCTIONS *
    *****************************************/

    template<typename T, typename D>
    static constexpr inline nion<T, D> asin(const nion<T, D> &z) noexcept {

        // get the polar form of the nion
        T r = real(z);
        nion<T, D> i = z.imag();

        // make unit vector
        T i_norm = i.norm();
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_norm > epsilon)
            i /= i_norm;

        // compute the inv sine of the nion
        return -i * log(sqrt(static_cast<T>(1) - sqr(z)) + (i * r) - i_norm);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> acos(const nion<T, D> &z) noexcept {
        return static_cast<T>(M_PI_2l) - asin(z);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> atan(const nion<T, D> &z) noexcept {
        // get the polar form of the nion
        T r = real(z);
        nion<T, D> i = z.imag();

        // make unit vector
        T i_norm = i.norm();
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_norm > epsilon)
            i /= i_norm;

        // compute the inv tangent of the nion:
        T one = static_cast<T>(1);
        return static_cast<T>(0.5) * i * (-log((one - i_norm) + (i * r) ) + log((one + i_norm) + (i * -r) ));
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> acot(const nion<T, D> &z) noexcept {
        return static_cast<T>(M_PI_2l) - atan(z);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> asec(const nion<T, D> &z) noexcept {
        return acos(inv(z));
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> acsc(const nion<T, D> &z) noexcept {
        return asin(inv(z));
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> atan2(const nion<T, D> &y, const nion<T, D> &x) noexcept {
        return atan(y / x);
    }

    /***************************
     *   NION GAMMA FUNCTION   *
     **************************/

    template<typename T, typename D>
    static constexpr inline nion<T, D> gamma(const nion<T, D> &z) noexcept {
        // compute the gamma function of the nion
        return exp(-z) * sqrt(inv(z))
               * pow(static_cast<T>(1) / (static_cast<T>(12) * z - inv(static_cast<T>(10) * z)) + z, z)
               * std::sqrt(static_cast<T>(2) * static_cast<T>(M_PIl));
    }
}
