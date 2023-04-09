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

    template<typename T, unsigned long int N, typename D>
    constexpr inline nion<T,N,D>::nion(T* vals, D size) : size_(size) {

        /// check if the degree is greater than zero and less than the maximum size
        ASSERT(size_ > 0, "The size of the nion must be greater than zero.");
        ASSERT(size_ <= N, "The size of the nion is too large. "
                                          "consider increasing template parameter, N.");

        /// copy the values into the nion
        memcpy(elem_, vals, size_ * sizeof(T));
    }

    template<typename T, unsigned long int N, typename D>
    constexpr inline nion<T,N,D>::nion(const std::initializer_list<T> &vals) : size_(vals.size()) {
        /// check if the degree is greater than zero
        ASSERT(size_ > 0, "The size of the nion must be greater than zero.");
        ASSERT(size_ <= N, "The size of the nion is too large. "
                                          "consider increasing template parameter, N.");

        /// copy the values into the nion
        memcpy(elem_, vals.begin(), size_ * sizeof(T));
    }

    template<typename T, unsigned long int N, typename D>
    constexpr inline nion<T,N,D>::nion(D size) : size_(size) {
        /// check if the degree is greater than zero
        ASSERT(size_ > 0, "The size of the nion must be greater than zero.");
        ASSERT(size_ <= N, "The size of the nion is too large. "
                                          "consider increasing template parameter, N.");

        /// initialize the values to zero
        zero();
    }

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline nion<T,N,D>::nion(const nion<T,M,D> &other) : size_(other.size_) {
        /// copy the values into the nion
        memcpy(elem_, other.elem_, size_ * sizeof(T));
    }

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline nion<T,N,D>::nion(nion<T,M,D> &&other) noexcept: size_(other.size_) {
        /// copy the values into the nion
        memcpy(elem_, other.elem_, size_ * sizeof(T));
    }

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline nion<T,N,D>::nion(S realVal, D size) : size_(size) {
        // check if the degree is greater than zero
        ASSERT(size_ > 0, "The degree of the nion must be greater than zero.");
        ASSERT(size_ <= N, "The size of the nion is too large. "
                                          "consider increasing template parameter, N.");

        // initialize the values to zero
        zero();

        // set the real part
        elem_[0] = static_cast<T>(realVal);
    }

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M, unsigned long int P>
    constexpr inline nion<T,N,D> nion<T,N,D>::make_pair(const nion<T,M,D> &a, const nion<T,P,D> &b) {

        /// initialize the nion pair
        nion<T,N,D> pair;

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

    template<typename T, unsigned long int N, typename D>
    constexpr inline void nion<T,N,D>::resize(int size) {
        ASSERT(size > 0, "new nion size must be greater than zero");
        size_ = size;

        // set the new elements to zero
        memset(elem_ + size_, 0, (size - size_) * sizeof(T));
    }

    /************************************
    *         ASSIGNMENT OPERATORS
    *************************************/
    
    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline nion<T,N,D> &nion<T,N,D>::operator=(const nion<T,M,D> &other) {
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

    template<typename T, unsigned long int N, typename D>
    constexpr nion<T,N,D> &nion<T,N,D>::operator=(const std::initializer_list<T> &vals) {

        /// set the size of the nion
        size_ = vals.size();

        /// copy the values into the nion
        memcpy(elem_, vals.begin(), size_ * sizeof(T));
        return *this; // return the nion
    }

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline nion<T,N,D> &nion<T,N,D>::operator=(nion<T,M,D> &&other)  noexcept {
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

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline nion<T,N,D> &nion<T,N,D>::operator=(S scalar) {
        /// check if the nion is initialized
        if (size_ <= 0) size_ = 1; // set the degree

        zero(); // set the nion to zero
        elem_[0] = static_cast<T>(scalar); // set the real component
        return *this; // return the nion
    }

    /************************************
    *  ASSIGNMENT AND ADDITION OPERATORS
    *************************************/

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline void nion<T,N,D>::operator+=(const nion<T,M,D> &other) {

        // add the components of the other nion to this nion.
        if (size_ >= other.size_) {
            for (D i = 0; i < other.size_; i++) elem_[i] += other.elem_[i];
        } else {
            for (D i = 0; i < size_; i++) elem_[i] += other.elem_[i];

            // copy the remaining values
            memcpy(elem_ + size_, other.elem_ + size_, (other.size_ - size_) * sizeof(T));

            // set the new size of the nion
            size_ = other.size_;
        }

    }

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline void nion<T,N,D>::operator-=(const nion<T,M,D> &other) {

        // subtract the components of the other nion to this nion.
        if (size_ >= other.size_) {
            for (D i = 0; i < other.size_; i++) elem_[i] -= other.elem_[i];
        } else {
            for (D i = 0; i < size_; i++) elem_[i] -= other.elem_[i];

            // copy the remaining values and negate them
            for (D i = size_; i < other.size_; i++) elem_[i] = -other.elem_[i];

            // set the new size of the nion
            size_ = other.size_;
        }
    }

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline void nion<T,N,D>::operator+=(S scalar) {
        elem_[0] += static_cast<T>(scalar);
    }

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline void nion<T,N,D>::operator-=(S scalar) {
        elem_[0] -= static_cast<T>(scalar);
    }

    /************************************
    *        ADDITION OPERATORS
    *************************************/

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline nion<T,N,D> nion<T,N,D>::operator+(const nion<T,M,D> &other) const {

        nion<T,N,D> sum; // create a nion to store the sum

        if (size_ >= other.size_) {
            // set the size of the sum
            sum.size_ = size_;

            // add the components of the other nion and this nion.
            for (D i = 0; i < other.size_; i++) sum.elem_[i] = elem_[i] + other.elem_[i];

            // copy the remaining values of this nion
            memcpy(sum.elem_ + other.size_, elem_ + other.size_, (size_ - other.size_) * sizeof(T));
        } else {
            // set the size of the sum
            sum.size_ = other.size_;

            // add the components of the other nion and this nion.
            for (D i = 0; i < size_; i++) sum.elem_[i] = elem_[i] + other.elem_[i];

            // copy the remaining values of the other nion
            memcpy(sum.elem_ + size_, other.elem_ + size_, (other.size_ - size_) * sizeof(T));
        }

        return sum;
    }

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline nion<T,N,D> nion<T,N,D>::operator-(const nion<T,M,D> &other) const {

        nion<T,N,D> diff; // create a nion to store the difference

        if (size_ >= other.size_) {
            // set the size of the difference
            diff.size_ = size_;

            // subtract the components of the other nion and this nion.
            for (D i = 0; i < other.size_; i++) diff.elem_[i] = elem_[i] - other.elem_[i];

            // copy the remaining values of this nion
            memcpy(diff.elem_ + other.size_, elem_ + other.size_, (size_ - other.size_) * sizeof(T));
        } else {
            // set the size of the difference
            diff.size_ = other.size_;

            // subtract the components of the other nion and this nion.
            for (D i = 0; i < size_; i++) diff.elem_[i] = elem_[i] - other.elem_[i];

            // copy the remaining values of the other nion and negate them
            for (D i = size_; i < other.size_; i++) diff.elem_[i] = -other.elem_[i];
        }

        return diff;
    }

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline nion<T,N,D> nion<T,N,D>::operator+(S scalar) const {
        nion<T,N,D> sum(size_); // create a nion to store the sum

        /// copy the components of this nion to the sum nion.
        memcpy(sum.elem_, elem_, size_ * sizeof(T));

        /// add the scalar to the real component of the sum nion.
        sum.elem_[0] += static_cast<T>(scalar);
        return sum;
    }

    template<typename T, unsigned long int N, typename D, typename S>
    extern constexpr inline nion<T,N,D> operator+(S scalar, const nion<T,N,D> &z) {
        return z + static_cast<T>(scalar);
    }

    template<typename T, unsigned long int N, typename D>
    constexpr inline nion<T,N,D> nion<T,N,D>::operator++() {
        elem_[0]++;
        return *this;
    }

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline nion<T,N,D> nion<T,N,D>::operator-(S scalar) const {
        return *this + (static_cast<T>(-scalar));
    }

    template<typename T, unsigned long int N, typename D, typename S>
    extern constexpr inline nion<T,N,D> operator-(S scalar, const nion<T,N,D> &z) {
        return -z + static_cast<T>(scalar);
    }

    template<typename T, unsigned long int N, typename D>
    constexpr inline nion<T,N,D> nion<T,N,D>::operator--() {
        elem_[0]--;
        return *this;
    }

    /******************************************
    *  ASSIGNMENT AND MULTIPLICATION OPERATORS
    *******************************************/

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline void nion<T,N,D>::operator*=(const nion<T,M,D> &other) {
        *this = *this * other;
    }

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline void nion<T,N,D>::operator*=(S scalar) {
        T product_scalar = static_cast<T>(scalar);
        for (D i = 0; i < size_; i++) elem_[i] *= product_scalar;

    }

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline void nion<T,N,D>::operator/=(const nion<T,M,D> &other) {
        *this = *this / other;
    }

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline void nion<T,N,D>::operator/=(S scalar) {
        T quotient_scalar = static_cast<T>(scalar);
        for (D i = 0; i < size_; i++) elem_[i] /= quotient_scalar;

    }

    /******************************************
    *        MULTIPLICATION OPERATORS
    *******************************************/

    template<typename T, unsigned long int N, typename D>
    constexpr inline nion<T,N,D> nion<T,N,D>::conj() const {
        nion<T,N,D> conjugate = *this; // copy this nion

        // negate all components except the first
        for (D i = 1; i < size_; i++) conjugate.elem_[i] *= -1;

        return conjugate;
    }

    template<typename T, unsigned long int N, typename D>
    constexpr inline nion<T,N,D> nion<T,N,D>::operator-() const {
        nion<T,N,D> negated = *this; // copy this nion

        // negate all components
        for (D i = 0; i < size_; i++) negated.elem_[i] *= -1;

        return negated;
    }

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline nion<T,N,D> nion<T,N,D>::operator*(const nion<T,M,D> &other) const {

        if (size_ == other.size_) {
            nion<T,N,D> product; // create a nion to store the product
            product.size_ = size_;
            
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
        } else { // If the matrix is not square, we need to do the recursive product.

            // if any of the sizes are 1, then the product is just the scalar product.
            if (size_ == 1) return other * elem_[0];
            if (other.size_ == 1) return *this * other.elem_[0];
        }

        /// ********** Recursive Product ********** ///

        nion<T,N-(N>>1),D> a, b;
        nion<T,M-(M>>1),D> c, d; // the four halves of the nions

        // the elements of the first half of the nion
        std::size_t this_half_size = size_ - (size_ >> 1);
        std::size_t other_half_size = other.size_ - (other.size_ >> 1);

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

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline nion<T,N,D> nion<T,N,D>::operator*(S scalar) const {
        nion<T,N,D> product;
        product.size_ = size_;
        T product_scalar = static_cast<T>(scalar);
        
        // compute the product of each element of the nion with the scalar
        for (D i = 0; i < size_; i++) product.elem_[i] = elem_[i] * product_scalar;
        
        return product;
    }

    template<typename T, unsigned long int N, typename D, typename S>
    extern constexpr inline nion<T,N,D> operator*(S scalar, const nion<T,N,D> &z) {
        static_assert(std::is_arithmetic_v<S>, "scalar must be arithmetic");
        return z * static_cast<T>(scalar);
    }

    template<typename T, unsigned long int N, typename D>
    constexpr inline T nion<T,N,D>::abs() const {
        T absVal = 0;
        
        for (D i = 0; i < size_; i++) absVal += elem_[i] * elem_[i];
        
        return absVal;
    }

    template<typename T, unsigned long int N, typename D>
    constexpr inline T nion<T,N,D>::norm() const {
        return std::sqrt(abs());
    }

    template<typename T, unsigned long int N, typename D>
    constexpr inline nion<T,N,D> nion<T,N,D>::inv() const {

        T absolute = abs();
        nion<T,N,D> inverse;
        inverse.size_ = size_;
        T inverse_scalar = static_cast<T>(1) / absolute;
        
        for (D i = 0; i < size_; i++) inverse.elem_[i] = -elem_[i] * inverse_scalar;
        
        inverse.elem_[0] = elem_[0] * inverse_scalar;
        return inverse;
    }

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline nion<T,N,D> nion<T,N,D>::operator/(const nion<T,M,D> &other) const {
        return *this * other.inv();
    }

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline nion<T,N,D> nion<T,N,D>::operator/(S scalar) const {
        return *this * (static_cast<T>(1 / scalar));
    }

    template<typename T, unsigned long int N, typename D, typename S>
    extern constexpr inline nion<T,N,D> operator/(S scalar, const nion<T,N,D> &z) {
        return z.inv() * static_cast<T>(scalar);
    }


    /******************************************
    *            ACCESSOR FUNCTIONS
    *******************************************/

    template<typename T, unsigned long int N, typename D>
    constexpr inline T nion<T,N,D>::operator[](D index) const {
        return this->elem_[index];
    }

    template<typename T, unsigned long int N, typename D>
    constexpr inline T nion<T,N,D>::real() const { return elem_[0]; }

    template<typename T, unsigned long int N, typename D>
    constexpr inline nion<T,N,D> nion<T,N,D>::imag() const {
        nion<T,N,D> imag = *this;
        imag.elem_[0] = 0;
        return imag;
    }

    /******************************************
    *            COMPARISON OPERATORS
    *******************************************/

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline bool nion<T,N,D>::operator==(const nion<T,M,D> &other) const {
        if (this == &other) return true;
        if (size_ != other.size_) return false;
        
        for (D i = 0; i < size_; i++) 
            if (!value_is_similar(elem_[i], other.elem_[i])) return false;
        return true;
    }

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline bool nion<T,N,D>::operator!=(const nion<T,M,D> &other) const {
        if (this == &other) return false;
        if (size_ != other.size_) return true;
        
        for (D i = 0; i < size_; i++)
            if (!value_is_similar(elem_[i], other.elem_[i])) return true;
        return false;
    }

    template<typename T, unsigned long int N, typename D>
    constexpr inline T nion<T,N,D>::rotate_real() const {
        // this yields the shortest rotation of the nion onto the real axis while preserving the norm and the sign of the real component
        // this is useful for comparing nions with different degrees
        return copysign(norm(), real());
    }

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline bool nion<T,N,D>::operator>(const nion<T,M,D> &other) const {
        // nions with degree > 1 are not ordered, but we can arbitrarily order them by their rotation to the real line
        // this is not a good idea, but it's better than nothing (and it's what I'm doing for now)
        return rotate_real() > other.rotate_real();
    }

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline bool nion<T,N,D>::operator<(const nion<T,M,D> &other) const {
        // nions with degree > 1 are not ordered, but we can arbitrarily order them by their rotation to the real line
        // this is not a good idea, but it's better than nothing (and it's what I'm doing for now)
        return rotate_real() < other.rotate_real();
    }

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline bool nion<T,N,D>::operator>=(const nion<T,M,D> &other) const {
        if (rotate_real() > other.rotate_real())
            return true;
        return *this == *other;
    }

    template<typename T, unsigned long int N, typename D>
    template<unsigned long int M>
    constexpr inline bool nion<T,N,D>::operator<=(const nion<T,M,D> &other) const {
        if (rotate_real() < other.rotate_real())
            return true;
        return *this == *other;
    }

    template<typename T, unsigned long int N, typename D, typename S>
    constexpr inline bool value_is_similar(const T a, const S b){
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        return std::fabs(a - static_cast<T>(b)) <= epsilon;
    }

    template<typename T, unsigned long int N, typename D>
    constexpr inline bool nion<T,N,D>::is_real() const {
        for (D i = 1; i < size_; i++)
            if (!value_is_similar(elem_[i], 0)) return false;
        return true;
    }

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline bool nion<T,N,D>::operator==(S scalar) const {
        if (!value_is_similar(real(), static_cast<T>(scalar))) return false;
        return is_real();
    }

    template<typename T, unsigned long int N, typename D, typename S>
    constexpr inline bool operator==(S scalar, const nion<T,N,D> &z) {
        return z == static_cast<T>(scalar);
    }

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline bool nion<T,N,D>::operator!=(S scalar) const {
        if (!value_is_similar(real(), static_cast<T>(scalar))) {
            return true;
        }
        return !is_real();
    }

    template<typename T, unsigned long int N, typename D, typename S>
    constexpr inline bool operator!=(S scalar, const nion<T,N,D> &z) {
        return z != static_cast<T>(scalar);
    }

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline bool nion<T,N,D>::operator>(S scalar) const{
        return rotate_real() > static_cast<T>(scalar);
    }

    template<typename T, unsigned long int N, typename D, typename S>
    constexpr inline bool operator>(S scalar, const nion<T,N,D> &z) {
        return z < static_cast<T>(scalar);
    }

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline bool nion<T,N,D>::operator<(S scalar) const{
        return rotate_real() < static_cast<T>(scalar);
    }

    template<typename T, unsigned long int N, typename D, typename S>
    constexpr inline bool operator<(S scalar, const nion<T,N,D> &z) {
        return z > static_cast<T>(scalar);
    }

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline bool nion<T,N,D>::operator>=(S scalar) const{
        if (*this == nion<T,N,D>(static_cast<T>(scalar), this->size_))
            return true;
        return rotate_real() > static_cast<T>(scalar);
    }

    template<typename T, unsigned long int N, typename D, typename S>
    constexpr inline bool operator>=(S scalar, const nion<T,N,D> &z) {
        return z <= static_cast<T>(scalar);
    }

    template<typename T, unsigned long int N, typename D>
    template<typename S>
    constexpr inline bool nion<T,N,D>::operator<=(S scalar) const{
        if (*this == nion<T,N,D>(static_cast<T>(scalar), this->size_))
            return true;
        return rotate_real() < static_cast<T>(scalar);
    }

    template<typename T, unsigned long int N, typename D, typename S>
    constexpr inline bool operator<=(S scalar, const nion<T,N,D> &z) {
        return z >= static_cast<T>(scalar);
    }

    /******************************************
    *            STREAMING OPERATORS
    *******************************************/

    template<typename T, unsigned long int N, typename D>
    std::ostream &operator<<(std::ostream &os, const nion<T,N,D> &z) {
        T component = z.elem_[0];
        os << "(" << component;

        for (D i = 1; i < z.size_; i++) {
            component = z.elem_[i];
            os << "," << component;
        }
        os << ")";
        return os;
    }

    template<typename T, unsigned long int N, typename D>
    std::istream &operator>>(std::istream &is, nion<T,N,D> &z) {
        for (D i = 0; i < z.size_; i++) {
            is >> z.elem_[i];
        }
        return is;
    }

    /**
    * @brief Converts a nion to a string.
    * @return The string representation of the nion.
    */
    template<typename T, unsigned long int N, typename D>
    inline std::string nion<T,N,D>::to_string() const {
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

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline T real(const nion<T,N,D> &z) {
        return z.real();
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> imag(const nion<T,N,D> &z) {
        return z.imag();
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> conj(const nion<T,N,D> &z) {
        return z.conj();
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline T abs(const nion<T,N,D> &z) {
        return z.abs();
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline T norm(const nion<T,N,D> &z) {
        return z.norm();
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> inv(const nion<T,N,D> &z) {
        return z.inv();
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline T dot(const nion<T,N,D> &lhs, const nion<T,N,D> &rhs) {
        T dotProduct = 0;
        D minDegree = std::min(lhs.size_, rhs.size_);
        for (D i = 0; i < minDegree; i++) {
            dotProduct += lhs.elem_[i] * rhs.elem_[i];
        }
        return dotProduct;
    }

    /*
    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D>
    cross(const nion<T,N,D> &lhs, const nion<T,N,D> &rhs){} //TODO: implement cross product

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D>
    wedge(const nion<T,N,D> &lhs, const nion<T,N,D> &rhs){} //TODO: implement wedge product

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D>
    outer(const nion<T,N,D> &lhs, const nion<T,N,D> &rhs){} //TODO: implement outer product
     */


    /****************************
    *  NION ALGEBRAIC FUNCTIONS *
    *****************************/

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> exp(const nion<T,N,D> &z) {

        // get polar form of nion
        T r = z.real();
        nion<T,N,D> i = z.imag();

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);

        // compute exponential of nion
        T cos_theta;
        T sin_theta;
        T exp_r = std::exp(r);

        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
        if (i_abs < denorm_min)
            return nion<T,N,D>(exp_r, z.size_);

        // compute exponential of nion
        cos_theta = std::cos(i_norm);
        sin_theta = std::sin(i_norm);

        return i*(exp_r * sin_theta / i_norm) + exp_r * cos_theta;

    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> log(const nion<T,N,D> &z) {

        // get polar form of nion
        T r = z.real();
        nion<T,N,D> i = z.imag();

        // compute norms
        T z_abs = z.abs();
        T z_norm = std::sqrt(z_abs);
        T i_abs = z_abs - r * r;
        T i_norm = std::sqrt(z_abs - r * r);
        T theta = std::atan2(i_norm, r);

        // compute natural logarithm of nion
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
        if (i_abs <= denorm_min)
            return std::log(z_norm) + i * theta;
        else
            return i * (theta / i_norm) + std::log(z_norm);
    }

    template<typename T, unsigned long int N, typename D, typename S>
    extern constexpr inline nion<T,N,D> pow(const nion<T,N,D> &base, S power) {
        return exp(power * log(base));
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> pow(const nion<T,N,D> &base, const nion<T,N,D> &power) {
        return exp(power * log(base));
    }

    template<typename T, unsigned long int N, typename D, typename S>
    extern constexpr inline nion<T,N,D> pow(S base, const nion<T,N,D> &power) {
        return exp(power * log(base));
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> sqr(const nion<T,N,D> &base) {
        // get polar form of nion
        T r = base.real();
        nion<T,N,D> i = base.imag();

        // compute norms
        T base_abs = base.abs();
        T base_norm = std::sqrt(base_abs);

        // make unit vector
        T i_abs = base_abs - r * r;
        T i_norm = std::sqrt(i_abs);

        T power_t = static_cast<T>(2);
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();

        T x2 = r * r;
        T y2 = i_abs;
        if (x2 + y2 <= denorm_min)
            return nion<T,N,D>(base_norm * base_norm, base.size_); // if base is zero return zero (0^2 = 0)

        T denom = static_cast<T>(1) / (x2 + y2);
        T cos_2theta = (x2 - y2) * denom;
        T sin_2theta = static_cast<T>(2) * r * i_norm * denom;

        // compute power of nion
        if (i_abs <= denorm_min)
            return std::pow(base_norm, power_t) * (cos_2theta + i * sin_2theta);
        else
            return std::pow(base_norm, power_t) * (cos_2theta + i * (sin_2theta / i_norm));
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> sqrt(const nion<T,N,D> &z) {
        return pow(z, static_cast<T>(0.5));
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> cbrt(const nion<T,N,D> &z) {
        return pow(z, static_cast<T>(1.0) / static_cast<T>(3.0));
    }

    /*******************************************
    *  NION HYPERBOLIC TRIGONOMETRIC FUNCTIONS *
    ********************************************/

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> sinh(const nion<T,N,D> &z) {
        // get polar form of nion
        nion<T,N,D> i = z.imag();

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);

        // calculate scalars
        T e_z = std::exp(z.elem_[0]) / static_cast<T>(2);
        T e_mz = std::exp(-z.elem_[0]) / static_cast<T>(2);

        // compute exponential of nion
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
        if (i_abs <= denorm_min) {
            nion<T,N,D> sin_nion = i * ((e_z + e_mz) * std::sin(i_norm));
            sin_nion += std::cos(i_norm) * (e_z - e_mz);
            return sin_nion;
        } else {
            nion<T,N,D> sin_nion = i * ((e_z + e_mz) * std::sin(i_norm) / i_norm);
            sin_nion += std::cos(i_norm) * (e_z - e_mz);
            return sin_nion;
        }

    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> cosh(const nion<T,N,D> &z) {
        // get polar form of nion
        nion<T,N,D> i = z.imag();

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);

        // calculate scalars
        T e_z = std::exp(z.elem_[0]) / static_cast<T>(2);
        T e_mz = std::exp(-z.elem_[0]) / static_cast<T>(2);

        // compute exponential of nion
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
        if (i_abs <= denorm_min) {
            nion<T,N,D> cos_nion = i * ((e_z - e_mz) * std::sin(i_norm));
            cos_nion += std::cos(i_norm) * (e_z + e_mz);
            return cos_nion;
        } else {
            nion<T,N,D> cos_nion = i * ((e_z - e_mz) * std::sin(i_norm) / i_norm);
            cos_nion += std::cos(i_norm) * (e_z + e_mz);
            return cos_nion;
        }
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> tanh(const nion<T,N,D> &z) {
        return (exp(z*2) - 1) / (exp(z*2) + 1);
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> coth(const nion<T,N,D> &z) {
        return tanh(z).inv();
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> sech(const nion<T,N,D> &z) {
        return cosh(z).inv();
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> csch(const nion<T,N,D> &z) {
        return sinh(z).inv();
    }

    /********************************
    *  NION TRIGONOMETRIC FUNCTIONS *
    *********************************/


    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> sin(const nion<T,N,D> &z) {
        // get the polar form of the nion
        T r = real(z);
        nion<T,N,D> i = imag(z);

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);

        // compute the sine of the nion
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
        if (i_abs <= denorm_min)
            return i * (std::sinh(i_norm) * std::cos(r)) + std::sin(r) * std::cosh(i_norm);
        else
            return i * (std::sinh(i_norm) * std::cos(r) / i_norm) + std::sin(r) * std::cosh(i_norm);
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> cos(const nion<T,N,D> &z) {
        // get the polar form of the nion
        T r = real(z);
        nion<T,N,D> i = imag(z);

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);

        // compute the cosine of the nion
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
        if (i_abs <= denorm_min)
            return -i * (std::sinh(i_norm) * std::sin(r)) + std::cos(r) * std::cosh(i_norm);
        else
            return -i * (std::sinh(i_norm) * std::sin(r) / i_norm) + std::cos(r) * std::cosh(i_norm);
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> tan(const nion<T,N,D> &z) {
        // get the polar form of the nion
        T r = real(z);
        nion<T,N,D> i = imag(z);

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);

        // compute the tangent of the nion
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
        if (i_norm <= denorm_min)
            return (std::tan(r) + i * std::tanh(i_norm)) / (static_cast<T>(1) - i * (std::tan(r) * std::tanh(i_norm)));
        else
            return (std::tan(r) + i * (std::tanh(i_norm) / i_norm)) / (static_cast<T>(1) - i * (std::tan(r) / i_norm * std::tanh(i_norm)));
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> cot(const nion<T,N,D> &z) {
        return tan(z).inv();
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> sec(const nion<T,N,D> &z) {
        return cos(z).inv();
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> csc(const nion<T,N,D> &z) {
        return sin(z).inv();
    }

    /***************************************************
    *  NION INVERSE HYPERBOLIC TRIGONOMETRIC FUNCTIONS *
    ****************************************************/

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> asinh(const nion<T,N,D> &z) {
        return log(z + sqrt(sqr(z) + static_cast<T>(1)));
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> acosh(const nion<T,N,D> &z) {
        // compute the inverse hyperbolic cosine of the nion
        return 2 * log(sqrt((z-1) * static_cast<T>(0.5)) + sqrt((z+1) * static_cast<T>(0.5)));

    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> atanh(const nion<T,N,D> &z) {
        return (log(static_cast<T>(1) + z) - log(static_cast<T>(1) - z)) * static_cast<T>(0.5);
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> acoth(const nion<T,N,D> &z) {
        return (log(static_cast<T>(1) + inv(z)) - log(static_cast<T>(1) - inv(z))) * static_cast<T>(0.5);
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> asech(const nion<T,N,D> &z) {
        return log(sqrt(1/sqr(z) - static_cast<T>(1)) + inv(z));
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> acsch(const nion<T,N,D> &z) {
        return log(sqrt(static_cast<T>(1) + 1/sqr(z)) + inv(z));
    }

    /****************************************
    *  NION INVERSE TRIGONOMETRIC FUNCTIONS *
    *****************************************/

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> asin(const nion<T,N,D> &z) {

        // get the polar form of the nion
        T r = real(z);
        nion<T,N,D> i = z.imag();

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
        if (i_abs > denorm_min)
            i /= i_norm;

        // compute the inv sine of the nion
        return -i * log(sqrt(static_cast<T>(1) - sqr(z)) + (i * r) - i_norm);
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> acos(const nion<T,N,D> &z) {
        return static_cast<T>(M_PI_2) - asin(z);
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> atan(const nion<T,N,D> &z) {
        // get the polar form of the nion
        T r = real(z);
        nion<T,N,D> i = z.imag();

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);
        constexpr T denorm_min = std::numeric_limits<T>::denorm_min();
        if (i_abs > denorm_min)
            i /= i_norm;

        // compute the inv tangent of the nion:
        T one = static_cast<T>(1);
        return static_cast<T>(0.5) * i * (-log((one - i_norm) + (i * r) ) + log((one + i_norm) + (i * -r) ));
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> acot(const nion<T,N,D> &z) {
        return static_cast<T>(M_PI_2) - atan(z);
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> asec(const nion<T,N,D> &z) {
        return acos(inv(z));
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> acsc(const nion<T,N,D> &z) {
        return asin(inv(z));
    }

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> atan2(const nion<T,N,D> &y, const nion<T,N,D> &x) {
        return atan(y / x);
    }

    /***************************
     *   NION GAMMA FUNCTION   *
     **************************/

    template<typename T, unsigned long int N, typename D>
    extern constexpr inline nion<T,N,D> gamma(const nion<T,N,D> &z) {
        // compute the gamma function of the nion
        return exp(-z) * sqrt(inv(z))
               * pow(static_cast<T>(1) / (static_cast<T>(12) * z - inv(static_cast<T>(10) * z)) + z, z)
               * std::sqrt(static_cast<T>(2) * static_cast<T>(M_PI));
    }

