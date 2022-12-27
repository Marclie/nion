#include <cmath>
#include <limits>
#include <type_traits>
#include <iostream>
#include <stdexcept>
#include "nion.hpp"

// max size of stack array for storing elements of nion. Defaults to 2^10 (nion degree of 10)
#ifndef NION_MAX_SIZE
#define NION_MAX_SIZE 1024
#endif

#define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << message << ": file=" << __FILE__ \
                      << ", line=" << __LINE__ << std::endl; \
            std::terminate(); \
        } \
    } while (false)

using Nion::nion;
namespace Nion {

    /***************************
    *  NION CONSTRUCTORS
    ***************************/

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(T* vals, D size) : size_(size) {

        /// check if the degree is greater than zero and less than the maximum size
        ASSERT(size_ > 0, "The size of the nion must be greater than zero.");
        ASSERT(size_ <= NION_MAX_SIZE, "The size of the nion is too large. "
                                          "consider increasing CMake variable NION_MAX_SIZE.");

        /// copy the values into the nion
        memcpy(elem_, vals, size_ * sizeof(T));
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(const std::initializer_list<T> &vals) : size_(vals.size()) {
        /// check if the degree is greater than zero
        ASSERT(size_ > 0, "The size of the nion must be greater than zero.");
        ASSERT(size_ <= NION_MAX_SIZE, "The size of the nion is too large. "
                                          "consider increasing CMake variable NION_MAX_SIZE.");

        /// copy the values into the nion
        memcpy(elem_, vals.begin(), size_ * sizeof(T));
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(D size) : size_(size) {
        /// check if the degree is greater than zero
        ASSERT(size_ > 0, "The size of the nion must be greater than zero.");
        ASSERT(size_ <= NION_MAX_SIZE, "The size of the nion is too large. "
                                          "consider increasing CMake variable NION_MAX_SIZE.");

        /// initialize the values to zero
        zero();
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(const nion<T, D> &other) : size_(other.size_) {
        /// copy the values into the nion
        memcpy(elem_, other.elem_, size_ * sizeof(T));
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(nion<T, D> &&other) noexcept: size_(other.size_) {
        /// copy the values into the nion
        memcpy(elem_, other.elem_, size_ * sizeof(T));
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline nion<T, D>::nion(S realVal, D size) : size_(size) {
        // check if the degree is greater than zero
        ASSERT(size_ > 0, "The degree of the nion must be greater than zero.");
        ASSERT(size_ <= NION_MAX_SIZE, "The size of the nion is too large. "
                                          "consider increasing CMake variable NION_MAX_SIZE.");

        // initialize the values to zero
        zero();

        // set the real part
        elem_[0] = static_cast<T>(realVal);
    }

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::make_pair(const nion<T, D> &a, const nion<T, D> &b) {

        /// initialize the nion pair
        nion<T, D> pair;

        /// set the size of the nion
        pair.size_ = a.size_ + b.size_;

        ASSERT(a.size_ > 0 && b.size_ > 0, "The sizes of the nion pair (a, b) must both be greater than zero.");
        ASSERT(pair.size_ <= NION_MAX_SIZE, "The size of the nion is too large. "
                                          "consider increasing CMake variable NION_MAX_SIZE.");
        
        /// copy the values into the nion
        memcpy(pair.elem_, a.elem_, a.size_ * sizeof(T));
        memcpy(pair.elem_ + a.size_, b.elem_, b.size_ * sizeof(T));

        return pair;
    }

    template<typename T, typename D>
    constexpr inline void nion<T, D>::resize(int size) {
        ASSERT(size > 0, "new nion size must be greater than zero");
        size_ = size;

        // set the new elements to zero
        memset(elem_ + size_, 0, (size - size_) * sizeof(T));
    }

    /************************************
    *         ASSIGNMENT OPERATORS
    *************************************/
    
    template<typename T, typename D>
    constexpr inline nion<T, D> &nion<T, D>::operator=(const nion<T, D> &other) {
        /// check if the nions are the same
        if (&other == this) {
            return *this; // return the nion
        }

        /// set the size
        size_ = other.size_;

        /// copy the values into the nion
        memcpy(elem_, other.elem_, size_ * sizeof(T));
        return *this; // return the nion
    }

    template<typename T, typename D>
    constexpr nion<T, D> &Nion::nion<T, D>::operator=(const std::initializer_list<T> &vals) {

        /// set the size of the nion
        size_ = vals.size();

        /// copy the values into the nion
        memcpy(elem_, vals.begin(), size_ * sizeof(T));
        return *this; // return the nion
    }

    template<typename T, typename D>
    constexpr inline nion<T, D> &nion<T, D>::operator=(nion<T, D> &&other)  noexcept {
        /// check if the nions are the same
        if (&other == this) {
            return *this; // return the nion
        }

        /// set the size
        size_ = other.size_;

        /// copy the values into the nion
        memcpy(elem_, other.elem_, size_ * sizeof(T));
        return *this; // return the nion
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline nion<T, D> &nion<T, D>::operator=(S scalar) {
        /// check if the nion is initialized
        if (size_ <= 0) size_ = 1; // set the degree

        zero(); // set the nion to zero
        elem_[0] = static_cast<T>(scalar); // set the real component
        return *this; // return the nion
    }

    /************************************
    *  ASSIGNMENT AND ADDITION OPERATORS
    *************************************/

    template<typename T, typename D>
    constexpr inline void nion<T, D>::operator+=(const nion<T, D> &other) {

        // add the components of the other nion to this nion.
        if (size_ >= other.size_) {
            #pragma unroll
            for (D i = 0; i < other.size_; i++) elem_[i] += other.elem_[i];
        } else {
            #pragma unroll
            for (D i = 0; i < size_; i++) elem_[i] += other.elem_[i];

            // copy the remaining values
            memcpy(elem_ + size_, other.elem_ + size_, (other.size_ - size_) * sizeof(T));

            // set the new size of the nion
            size_ = other.size_;
        }

    }

    template<typename T, typename D>
    constexpr inline void nion<T, D>::operator-=(const nion<T, D> &other) {

        // subtract the components of the other nion to this nion.
        if (size_ >= other.size_) {
            #pragma unroll
            for (D i = 0; i < other.size_; i++) elem_[i] -= other.elem_[i];
        } else {
            #pragma unroll
            for (D i = 0; i < size_; i++) elem_[i] -= other.elem_[i];

            // copy the remaining values and negate them
            #pragma unroll
            for (D i = size_; i < other.size_; i++) elem_[i] = -other.elem_[i];

            // set the new size of the nion
            size_ = other.size_;
        }
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline void nion<T, D>::operator+=(S scalar) {
        elem_[0] += static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline void nion<T, D>::operator-=(S scalar) {
        elem_[0] -= static_cast<T>(scalar);
    }

    /************************************
    *        ADDITION OPERATORS
    *************************************/

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator+(const nion<T, D> &other) const {

        nion<T, D> sum; // create a nion to store the sum

        if (size_ >= other.size_) {
            // set the size of the sum
            sum.size_ = size_;

            // add the components of the other nion and this nion.
            #pragma unroll
            for (D i = 0; i < other.size_; i++) sum.elem_[i] = elem_[i] + other.elem_[i];

            // copy the remaining values of this nion
            memcpy(sum.elem_ + other.size_, elem_ + other.size_, (size_ - other.size_) * sizeof(T));
        } else {
            // set the size of the sum
            sum.size_ = other.size_;

            // add the components of the other nion and this nion.
            #pragma unroll
            for (D i = 0; i < size_; i++) sum.elem_[i] = elem_[i] + other.elem_[i];

            // copy the remaining values of the other nion
            memcpy(sum.elem_ + size_, other.elem_ + size_, (other.size_ - size_) * sizeof(T));
        }

        return sum;
    }

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator-(const nion<T, D> &other) const {

        nion<T, D> diff; // create a nion to store the difference

        if (size_ >= other.size_) {
            // set the size of the difference
            diff.size_ = size_;

            // subtract the components of the other nion and this nion.
            #pragma unroll
            for (D i = 0; i < other.size_; i++) diff.elem_[i] = elem_[i] - other.elem_[i];

            // copy the remaining values of this nion
            memcpy(diff.elem_ + other.size_, elem_ + other.size_, (size_ - other.size_) * sizeof(T));
        } else {
            // set the size of the difference
            diff.size_ = other.size_;

            // subtract the components of the other nion and this nion.
            #pragma unroll
            for (D i = 0; i < size_; i++) diff.elem_[i] = elem_[i] - other.elem_[i];

            // copy the remaining values of the other nion and negate them
            #pragma unroll
            for (D i = size_; i < other.size_; i++) diff.elem_[i] = -other.elem_[i];
        }

        return diff;
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline nion<T, D> nion<T, D>::operator+(S scalar) const {
        nion<T, D> sum(size_); // create a nion to store the sum

        /// copy the components of this nion to the sum nion.
        memcpy(sum.elem_, elem_, size_ * sizeof(T));

        /// add the scalar to the real component of the sum nion.
        sum.elem_[0] += static_cast<T>(scalar);
        return sum;
    }

    template<typename T, typename D, typename S>
    extern constexpr inline nion<T, D> operator+(S scalar, const nion<T, D> &z) {
        return z + static_cast<T>(scalar);
    }

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator++() {
        elem_[0]++;
        return *this;
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline nion<T, D> nion<T, D>::operator-(S scalar) const {
        return *this + (static_cast<T>(-scalar));
    }

    template<typename T, typename D, typename S>
    extern constexpr inline nion<T, D> operator-(S scalar, const nion<T, D> &z) {
        return -z + static_cast<T>(scalar);
    }

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator--() {
        elem_[0]--;
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
        T product_scalar = static_cast<T>(scalar);
        #pragma unroll
        for (D i = 0; i < size_; i++) elem_[i] *= product_scalar;

    }

    template<typename T, typename D>
    constexpr inline void nion<T, D>::operator/=(const nion<T, D> &other) {
        *this = *this / other;
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline void nion<T, D>::operator/=(S scalar) {
        T quotient_scalar = static_cast<T>(scalar);
        #pragma unroll
        for (D i = 0; i < size_; i++) elem_[i] /= quotient_scalar;

    }

    /******************************************
    *        MULTIPLICATION OPERATORS
    *******************************************/

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::conj() const {
        nion<T, D> conjugate = *this; // copy this nion

        // negate all components except the first
        #pragma unroll
        for (D i = 1; i < size_; i++) conjugate.elem_[i] *= -1;

        return conjugate;
    }

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator-() const {
        nion<T, D> negated = *this; // copy this nion

        // negate all components
        #pragma unroll
        for (D i = 0; i < size_; i++) negated.elem_[i] *= -1;

        return negated;
    }

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator*(const nion<T, D> &other) const {

        if (size_ == other.size_) {
            nion<T, D> product; // create a nion to store the product
            product.size_ = size_;
            
            switch (size_) {
                case 1: // if this size is 1, then the product is just the scalar product.
                    product.elem_[0] = elem_[0] * other.elem_[0];
                    return product;
                    
// if not in debug mode, use the optimized code
#ifndef NION_DEBUG
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
        } else { // If the matrix is not square, we need to do the recursive product.

            // if any of the sizes are 1, then the product is just the scalar product.
            if (size_ == 1) return other * elem_[0];
            if (other.size_ == 1) return *this * other.elem_[0];
        }

        /// ********** Recursive Product ********** ///

        nion<T> a, b, c, d; // the four halves of the nions

        // the elements of the first half of the nion
        std::size_t this_half_size = size_ - size_ / 2;
        std::size_t other_half_size = other.size_ - other.size_ / 2;

        a.size_ = this_half_size; // the size of the first half of the nion
        b.size_ = size_ - this_half_size; // the size of the second half of the nion (can be one less than the first half)
        c.size_ = other_half_size; // same as a for the other nion
        d.size_ = other.size_ - other_half_size; // same as b for the other nion

        /// copy the elements of the nions into the halves ///
        memcpy(a.elem_, elem_, this_half_size * sizeof(T));
        memcpy(b.elem_, elem_ + this_half_size, (size_ - this_half_size) * sizeof(T));
        memcpy(c.elem_, other.elem_, other_half_size * sizeof(T));
        memcpy(d.elem_, other.elem_ + other_half_size, (other.size_ - other_half_size) * sizeof(T));

        /// calculate the cayley-dickson product
        return make_pair(
                (a * c) - (d.conj() * b), // add involution parameter for sign with macro? (split hypercomplex numbers)
                (d * a) + (b * c.conj())
        );
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline nion<T, D> nion<T, D>::operator*(S scalar) const {
        nion<T, D> product;
        product.size_ = size_;
        T product_scalar = static_cast<T>(scalar);
        
        // compute the product of each element of the nion with the scalar
        #pragma unroll
        for (D i = 0; i < size_; i++) product.elem_[i] = elem_[i] * product_scalar;
        
        return product;
    }

    template<typename T, typename D, typename S>
    extern constexpr inline nion<T, D> operator*(S scalar, const nion<T, D> &z) {
        static_assert(std::is_arithmetic_v<S>, "scalar must be arithmetic");
        return z * static_cast<T>(scalar);
    }

    template<typename T, typename D>
    constexpr inline T nion<T, D>::abs() const {
        T absVal = 0;
        
        #pragma unroll
        for (D i = 0; i < size_; i++) absVal += elem_[i] * elem_[i];
        
        return absVal;
    }

    template<typename T, typename D>
    constexpr inline T nion<T, D>::norm() const {
        return std::sqrt(abs());
    }

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::inv() const {
        constexpr T epsilon = std::numeric_limits<T>::epsilon();

        T absolute = abs();
        if (absolute < epsilon) {
            // if the absolute value is zero, then use the product definition of the absolute value.
            // zero divisors are possible in nions with degree >= 16, so we need to check for them.
            // this is a bit of a hack, and probably not useful. (like the rest of this library)

            // q* / |q| = q* / sqrt(q * q*)
            nion<T, D> inverse = this->conj() * pow((*this * this->conj()), static_cast<T>(-0.5)).real();
            return inverse;
        }

        nion<T, D> inverse;
        inverse.size_ = size_;
        T inverse_scalar = static_cast<T>(1) / absolute;
        
        #pragma unroll
        for (D i = 0; i < size_; i++) inverse.elem_[i] = -elem_[i] * inverse_scalar;
        
        inverse.elem_[0] = elem_[0] * inverse_scalar;
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
    extern constexpr inline nion<T, D> operator/(S scalar, const nion<T, D> &z) {
        return z.inv() * static_cast<T>(scalar);
    }


    /******************************************
    *            ACCESSOR FUNCTIONS
    *******************************************/

    template<typename T, typename D>
    constexpr inline T nion<T, D>::operator[](D index) const {
        return this->elem_[index];
    }

    template<typename T, typename D>
    constexpr inline T nion<T, D>::real() const { return elem_[0]; }

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::imag() const {
        nion<T> imag = *this;
        imag.elem_[0] = 0;
        return imag;
    }

    /******************************************
    *            COMPARISON OPERATORS
    *******************************************/

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::operator==(const nion<T, D> &other) const {
        if (this == &other) return true;
        if (size_ != other.size_) return false;
        
        #pragma unroll
        for (D i = 0; i < size_; i++) 
            if (!value_is_similar(elem_[i], other.elem_[i])) return false;
        return true;
    }

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::operator!=(const nion<T, D> &other) const {
        if (this == &other) return false;
        if (size_ != other.size_) return true;
        
        #pragma unroll
        for (D i = 0; i < size_; i++)
            if (!value_is_similar(elem_[i], other.elem_[i])) return true;
        return false;
    }

    template<typename T, typename D>
    constexpr inline T nion<T, D>::rotate_real() const {
        // this yields the shortest rotation of the nion onto the real axis while preserving the norm and the sign of the real component
        // this is useful for comparing nions with different degrees
        return copysign(norm(), real());
    }

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::operator>(const nion<T, D> &other) const {
        // nions with degree > 1 are not ordered, but we can arbitrarily order them by their rotation to the real line
        // this is not a good idea, but it's better than nothing (and it's what I'm doing for now)
        return rotate_real() > other.rotate_real();
    }

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::operator<(const nion<T, D> &other) const {
        // nions with degree > 1 are not ordered, but we can arbitrarily order them by their rotation to the real line
        // this is not a good idea, but it's better than nothing (and it's what I'm doing for now)
        return rotate_real() < other.rotate_real();
    }

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::operator>=(const nion<T, D> &other) const {
        if (rotate_real() > other.rotate_real())
            return true;
        return *this == *other;
    }

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::operator<=(const nion<T, D> &other) const {
        if (rotate_real() < other.rotate_real())
            return true;
        return *this == *other;
    }

    template<typename T, typename D, typename S>
    constexpr inline bool value_is_similar(const T a, const S b){
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        return std::fabs(a - static_cast<T>(b)) <= epsilon;
    }

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::is_real() const {
        #pragma unroll
        for (D i = 1; i < size_; i++)
            if (!value_is_similar(elem_[i], 0)) return false;
        return true;
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline bool nion<T, D>::operator==(S scalar) const {
        if (!value_is_similar(real(), static_cast<T>(scalar))) return false;
        return is_real();
    }

    template<typename T, typename D, typename S>
    constexpr inline bool operator==(S scalar, const nion<T, D> &z) {
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
    constexpr inline bool operator!=(S scalar, const nion<T, D> &z) {
        return z != static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline bool nion<T,D>::operator>(S scalar) const{
        return rotate_real() > static_cast<T>(scalar);
    }

    template<typename T, typename D, typename S>
    constexpr inline bool operator>(S scalar, const nion<T, D> &z) {
        return z < static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline bool nion<T,D>::operator<(S scalar) const{
        return rotate_real() < static_cast<T>(scalar);
    }

    template<typename T, typename D, typename S>
    constexpr inline bool operator<(S scalar, const nion<T, D> &z) {
        return z > static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline bool nion<T,D>::operator>=(S scalar) const{
        if (*this == nion<T, D>(static_cast<T>(scalar), this->size_))
            return true;
        return rotate_real() > static_cast<T>(scalar);
    }

    template<typename T, typename D, typename S>
    constexpr inline bool operator>=(S scalar, const nion<T, D> &z) {
        return z <= static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline bool nion<T,D>::operator<=(S scalar) const{
        if (*this == nion<T, D>(static_cast<T>(scalar), this->size_))
            return true;
        return rotate_real() < static_cast<T>(scalar);
    }

    template<typename T, typename D, typename S>
    constexpr inline bool operator<=(S scalar, const nion<T, D> &z) {
        return z >= static_cast<T>(scalar);
    }

    /******************************************
    *            STREAMING OPERATORS
    *******************************************/

    template<typename T, typename D>
    std::ostream &operator<<(std::ostream &os, const nion<T, D> &z) {
        T component = z.elem_[0];
        os << "(" << component;

        for (D i = 1; i < z.size_; i++) {
            component = z.elem_[i];
            os << "," << component;
        }
        os << ")";
        return os;
    }

    template<typename T, typename D>
    std::istream &operator>>(std::istream &is, nion<T, D> &z) {
        for (D i = 0; i < z.size_; i++) {
            is >> z.elem_[i];
        }
        return is;
    }

    /**
    * @brief Converts a nion to a string.
    * @return The string representation of the nion.
    */
    template<typename T, typename D>
    inline std::string nion<T, D>::to_string() const {
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

    template<typename T, typename D>
    extern constexpr inline T real(const nion<T, D> &z) {
        return z.real();
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> imag(const nion<T, D> &z) {
        return z.imag();
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> conj(const nion<T, D> &z) {
        return z.conj();
    }

    template<typename T, typename D>
    extern constexpr inline T abs(const nion<T, D> &z) {
        return z.abs();
    }

    template<typename T, typename D>
    extern constexpr inline T norm(const nion<T, D> &z) {
        return z.norm();
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> inv(const nion<T, D> &z) {
        return z.inv();
    }

    template<typename T, typename D>
    extern constexpr inline T dot(const nion<T, D> &lhs, const nion<T, D> &rhs) {
        T dotProduct = 0;
        D minDegree = std::min(lhs.size_, rhs.size_);
        #pragma unroll
        for (D i = 0; i < minDegree; i++) {
            dotProduct += lhs.elem_[i] * rhs.elem_[i];
        }
        return dotProduct;
    }

    /*
    template<typename T, typename D>
    extern constexpr inline nion<T, D>
    cross(const nion<T, D> &lhs, const nion<T, D> &rhs){} //TODO: implement cross product

    template<typename T, typename D>
    extern constexpr inline nion<T, D>
    wedge(const nion<T, D> &lhs, const nion<T, D> &rhs){} //TODO: implement wedge product

    template<typename T, typename D>
    extern constexpr inline nion<T, D>
    outer(const nion<T, D> &lhs, const nion<T, D> &rhs){} //TODO: implement outer product
     */


    /****************************
    *  NION ALGEBRAIC FUNCTIONS *
    *****************************/

    template<typename T, typename D>
    extern constexpr inline nion<T, D> exp(const nion<T, D> &z) {

        // get polar form of nion
        T r = z.real();
        nion<T, D> i = z.imag();

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);

        // compute exponential of nion
        T cos_theta;
        T sin_theta;
        T exp_r = std::exp(r);

        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_abs < epsilon)
            return nion<T, D>(exp_r, z.size_);

        // compute exponential of nion
        cos_theta = std::cos(i_norm);
        sin_theta = std::sin(i_norm);

        return i*(exp_r * sin_theta / i_norm) + exp_r * cos_theta;

    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> log(const nion<T, D> &z) {

        // get polar form of nion
        T r = z.real();
        nion<T, D> i = z.imag();

        // compute norms
        T z_abs = z.abs();
        T z_norm = std::sqrt(z_abs);
        T i_abs = z_abs - r * r;
        T i_norm = std::sqrt(z_abs - r * r);
        T theta = std::atan2(i_norm, r);

        // compute natural logarithm of nion
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_abs <= epsilon)
            return std::log(z_norm) + i * theta;
        else
            return i * (theta / i_norm) + std::log(z_norm);
    }

    template<typename T, typename D, typename S>
    extern constexpr inline nion<T, D> pow(const nion<T, D> &base, S power) {

        // get polar form of nion
        T r = base.real();
        nion<T, D> i = base.imag();

        // compute norms
        T base_abs = base.abs();
        T base_norm = std::sqrt(base_abs);

        // make unit vector
        T i_abs = base_abs - r * r;
        T i_norm = std::sqrt(i_abs);

        T power_t;
        T theta;

        T cos_ptheta;
        T sin_ptheta;

        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_abs <= epsilon) { // real nion
            if (std::signbit(r)) { // if base is negative
                power_t = static_cast<T>(power);
                theta = power_t * std::atan2(i_norm, r);
                cos_ptheta = std::cos(theta);
                sin_ptheta = std::sin(theta);

                return std::pow(base_norm, power_t) * (cos_ptheta + i * (sin_ptheta));
            } else {
                return nion<T, D>(std::pow(base_norm, static_cast<T>(power)), base.size_); // if base is positive
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
                    return nion<T, D>(1, z.size_);
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
    extern constexpr inline nion<T, D> pow(const nion<T, D> &base, const nion<T, D> &power) {
        return exp(power * log(base));
    }

    template<typename T, typename D, typename S>
    extern constexpr inline nion<T, D> pow(S base, const nion<T, D> &power) {
        return pow(nion<T, D>(static_cast<T>(base), power.size_), power);
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> sqr(const nion<T, D> &base) {
        // get polar form of nion
        T r = base.real();
        nion<T, D> i = base.imag();

        // compute norms
        T base_abs = base.abs();
        T base_norm = std::sqrt(base_abs);

        // make unit vector
        T i_abs = base_abs - r * r;
        T i_norm = std::sqrt(i_abs);

        T power_t = static_cast<T>(2);
        constexpr T epsilon = std::numeric_limits<T>::epsilon();

        T x2 = r * r;
        T y2 = i_norm * i_norm;
        if (x2 + y2 <= epsilon)
            return nion<T>(base_norm * base_norm, base.size_); // if base is zero return zero (0^2 = 0)

        T denom = static_cast<T>(1) / (x2 + y2);
        T cos_2theta = (x2 - y2) * denom;
        T sin_2theta = static_cast<T>(2) * r * i_norm * denom;

        // compute power of nion
        if (i_abs <= epsilon)
            return std::pow(base_norm, power_t) * (cos_2theta + i * sin_2theta);
        else
            return std::pow(base_norm, power_t) * (cos_2theta + i * (sin_2theta / i_norm));
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> sqrt(const nion<T, D> &z) {
        return pow(z, static_cast<T>(0.5));
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> cbrt(const nion<T, D> &z) {
        return pow(z, static_cast<T>(1.0) / static_cast<T>(3.0));
    }

    /*******************************************
    *  NION HYPERBOLIC TRIGONOMETRIC FUNCTIONS *
    ********************************************/

    template<typename T, typename D>
    extern constexpr inline nion<T, D> sinh(const nion<T, D> &z) {
        // get polar form of nion
        nion<T, D> i = z.imag();

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);

        // calculate scalars
        T e_z = std::exp(z.elem_[0]) / static_cast<T>(2);
        T e_mz = std::exp(-z.elem_[0]) / static_cast<T>(2);

        // compute exponential of nion
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_abs <= epsilon) {
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
    extern constexpr inline nion<T, D> cosh(const nion<T, D> &z) {
        // get polar form of nion
        nion<T, D> i = z.imag();

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);

        // calculate scalars
        T e_z = std::exp(z.elem_[0]) / static_cast<T>(2);
        T e_mz = std::exp(-z.elem_[0]) / static_cast<T>(2);

        // compute exponential of nion
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_abs <= epsilon) {
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
    extern constexpr inline nion<T, D> tanh(const nion<T, D> &z) {
        return (exp(z*2) - 1) / (exp(z*2) + 1);
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> coth(const nion<T, D> &z) {
        return tanh(z).inv();
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> sech(const nion<T, D> &z) {
        return cosh(z).inv();
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> csch(const nion<T, D> &z) {
        return sinh(z).inv();
    }

    /********************************
    *  NION TRIGONOMETRIC FUNCTIONS *
    *********************************/


    template<typename T, typename D>
    extern constexpr inline nion<T, D> sin(const nion<T, D> &z) {
        // get the polar form of the nion
        T r = real(z);
        nion<T, D> i = imag(z);

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);

        // compute the sine of the nion
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_abs <= epsilon)
            return i * (std::sinh(i_norm) * std::cos(r)) + std::sin(r) * std::cosh(i_norm);
        else
            return i * (std::sinh(i_norm) * std::cos(r) / i_norm) + std::sin(r) * std::cosh(i_norm);
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> cos(const nion<T, D> &z) {
        // get the polar form of the nion
        T r = real(z);
        nion<T, D> i = imag(z);

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);

        // compute the cosine of the nion
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_abs <= epsilon)
            return -i * (std::sinh(i_norm) * std::sin(r)) + std::cos(r) * std::cosh(i_norm);
        else
            return -i * (std::sinh(i_norm) * std::sin(r) / i_norm) + std::cos(r) * std::cosh(i_norm);
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> tan(const nion<T, D> &z) {
        // get the polar form of the nion
        T r = real(z);
        nion<T, D> i = imag(z);

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);

        // compute the tangent of the nion
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_norm <= epsilon)
            return (std::tan(r) + i * std::tanh(i_norm)) / (static_cast<T>(1) - i * (std::tan(r) * std::tanh(i_norm)));
        else
            return (std::tan(r) + i * (std::tanh(i_norm) / i_norm)) / (static_cast<T>(1) - i * (std::tan(r) / i_norm * std::tanh(i_norm)));
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> cot(const nion<T, D> &z) {
        return tan(z).inv();
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> sec(const nion<T, D> &z) {
        return cos(z).inv();
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> csc(const nion<T, D> &z) {
        return sin(z).inv();
    }

    /***************************************************
    *  NION INVERSE HYPERBOLIC TRIGONOMETRIC FUNCTIONS *
    ****************************************************/

    template<typename T, typename D>
    extern constexpr inline nion<T, D> asinh(const nion<T, D> &z) {
        return log(z + sqrt(sqr(z) + static_cast<T>(1)));
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> acosh(const nion<T, D> &z) {
        // compute the inverse hyperbolic cosine of the nion
        return 2 * log(sqrt((z-1) * static_cast<T>(0.5)) + sqrt((z+1) * static_cast<T>(0.5)));

    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> atanh(const nion<T, D> &z) {
        return (log(static_cast<T>(1) + z) - log(static_cast<T>(1) - z)) * static_cast<T>(0.5);
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> acoth(const nion<T, D> &z) {
        return (log(static_cast<T>(1) + inv(z)) - log(static_cast<T>(1) - inv(z))) * static_cast<T>(0.5);
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> asech(const nion<T, D> &z) {
        return log(sqrt(1/sqr(z) - static_cast<T>(1)) + inv(z));
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> acsch(const nion<T, D> &z) {
        return log(sqrt(static_cast<T>(1) + 1/sqr(z)) + inv(z));
    }

    /****************************************
    *  NION INVERSE TRIGONOMETRIC FUNCTIONS *
    *****************************************/

    template<typename T, typename D>
    extern constexpr inline nion<T, D> asin(const nion<T, D> &z) {

        // get the polar form of the nion
        T r = real(z);
        nion<T, D> i = z.imag();

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_abs > epsilon)
            i /= i_norm;

        // compute the inv sine of the nion
        return -i * log(sqrt(static_cast<T>(1) - sqr(z)) + (i * r) - i_norm);
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> acos(const nion<T, D> &z) {
        return static_cast<T>(M_PI_2l) - asin(z);
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> atan(const nion<T, D> &z) {
        // get the polar form of the nion
        T r = real(z);
        nion<T, D> i = z.imag();

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        if (i_abs > epsilon)
            i /= i_norm;

        // compute the inv tangent of the nion:
        T one = static_cast<T>(1);
        return static_cast<T>(0.5) * i * (-log((one - i_norm) + (i * r) ) + log((one + i_norm) + (i * -r) ));
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> acot(const nion<T, D> &z) {
        return static_cast<T>(M_PI_2l) - atan(z);
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> asec(const nion<T, D> &z) {
        return acos(inv(z));
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> acsc(const nion<T, D> &z) {
        return asin(inv(z));
    }

    template<typename T, typename D>
    extern constexpr inline nion<T, D> atan2(const nion<T, D> &y, const nion<T, D> &x) {
        return atan(y / x);
    }

    /***************************
     *   NION GAMMA FUNCTION   *
     **************************/

    template<typename T, typename D>
    extern constexpr inline nion<T, D> gamma(const nion<T, D> &z) {
        // compute the gamma function of the nion
        return exp(-z) * sqrt(inv(z))
               * pow(static_cast<T>(1) / (static_cast<T>(12) * z - inv(static_cast<T>(10) * z)) + z, z)
               * std::sqrt(static_cast<T>(2) * static_cast<T>(M_PIl));
    }
}
