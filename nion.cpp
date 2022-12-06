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

using Nion::nion;
namespace Nion {

    /***************************
    *  NION CONSTRUCTORS
    ***************************/

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(T* vals, D degree) : degree_(degree) {

        /// check if the degree is greater than zero
        ASSERT(degree_ > 0, "The degree of the nion must be greater than zero.");

        /// calculate memory size
        constexpr D alignment_requirement = std::alignment_of<T>::value; // alignment requirement of the type
        bytes_ = degree_ * sizeof(T); // total number of bytes

        /// calculate the number of bytes to add to the memory size to make it the LCM of the alignment requirement
        D min_power_of_two = (bytes_ | (bytes_ - 1)) + 1; // minimum power of two that is greater than or equal to 'bytes'
        alignment_ = std::max(min_power_of_two, alignment_requirement); // minimum alignment size

        /// allocate memory for the nion
        components_ = static_cast<T*>(aligned_alloc(alignment_, bytes_));

        /// copy the values into the nion
        memcpy(components_, vals, bytes_);
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(const std::initializer_list<T> &vals) : degree_(vals.size()) {
        /// check if the degree is greater than zero
        ASSERT(degree_ > 0, "The degree of the nion must be greater than zero.");

        /// calculate memory size
        constexpr D alignment_requirement = std::alignment_of<T>::value; // alignment requirement of the type
        bytes_ = degree_ * sizeof(T); // total number of bytes

        /// calculate the number of bytes to add to the memory size to make it the LCM of the alignment requirement
        D min_power_of_two = (bytes_ | (bytes_ - 1)) + 1; // minimum power of two that is greater than or equal to 'bytes'
        alignment_ = std::max(min_power_of_two, alignment_requirement); // minimum alignment size

        /// allocate memory for the nion
        components_ = static_cast<T*>(aligned_alloc(alignment_, bytes_));

        /// copy the values into the nion
        memcpy(components_, vals.begin(), bytes_);
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(D degree) : degree_(degree) {
        /// check if the degree is greater than zero
        ASSERT(degree_ > 0, "The degree of the nion must be greater than zero.");

        /// calculate memory size
        constexpr D alignment_requirement = std::alignment_of<T>::value; // alignment requirement of the type
        bytes_ = degree_ * sizeof(T); // total number of bytes

        /// calculate the number of bytes to add to the memory size to make it the LCM of the alignment requirement
        D min_power_of_two = (bytes_ | (bytes_ - 1)) + 1; // minimum power of two that is greater than or equal to 'bytes'
        alignment_ = std::max(min_power_of_two, alignment_requirement); // minimum alignment size

        /// allocate memory for the nion
        components_ = static_cast<T*>(aligned_alloc(alignment_, bytes_)); // keep uninitialized
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(const nion<T, D> &other) : degree_(other.degree_), bytes_(other.bytes_), alignment_(other.alignment_) {

        /// allocate memory for the nion
        components_ = static_cast<T*>(aligned_alloc(alignment_, bytes_));

        /// copy the values into the nion
        memcpy(components_, other.components_, bytes_);
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(nion<T, D> &&other) noexcept: degree_(other.degree_), components_(other.components_),
                                                               bytes_(other.bytes_), alignment_(other.alignment_) {
        other.components_ = nullptr;
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline nion<T, D>::nion(S realVal, D degree) : degree_(degree) {
        // check if the degree is greater than zero
        ASSERT(degree_ > 0, "The degree of the nion must be greater than zero.");

        /// calculate memory size
        constexpr D alignment_requirement = std::alignment_of<T>::value; // alignment requirement of the type
        bytes_ = degree_ * sizeof(T); // total number of bytes
        
        /// calculate the number of bytes to add to the memory size to make it the LCM of the alignment requirement
        D min_power_of_two = (bytes_ | (bytes_ - 1)) + 1; // minimum power of two that is greater than or equal to 'bytes'
        alignment_ = std::max(min_power_of_two, alignment_requirement); // minimum alignment size
        
        /// allocate memory for the nion
        components_ = static_cast<T*>(aligned_alloc(alignment_, bytes_));
        memset(components_, 0, bytes_);
        components_[0] = static_cast<T>(realVal);
    }

    template<typename T, typename D>
    constexpr inline nion<T, D>::nion(const nion<T, D> &a, const nion<T, D> &b) : degree_(a.degree_ + b.degree_), bytes_(a.bytes_ + b.bytes_) {
        ASSERT(a.degree_ > 0 && b.degree_ > 0, "The degrees of the nion pair (a, b) must both be greater than zero.");

        /// calculate memory size
        constexpr D alignment_requirement = std::alignment_of<T>::value; // alignment requirement of the type

        /// calculate the number of bytes to add to the memory size to make it the LCM of the alignment requirement
        D min_power_of_two = (bytes_ | (bytes_ - 1)) + 1; // minimum power of two that is greater than or equal to 'bytes'
        alignment_ = std::max(min_power_of_two, alignment_requirement); // minimum alignment size

        /// allocate memory for the nion
        components_ = static_cast<T*>(aligned_alloc(alignment_, bytes_));

        /// copy the values into the nion
        memcpy(components_, a.components_, a.bytes_);
        memcpy(components_ + a.degree_, b.components_, b.bytes_);
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

        /// No need to reallocate memory if this nion's degree is greater than or equal to the other nion's degree
        if (degree_ < other.degree_) {
            // reallocate memory for the nion
            free(components_);
            components_ = static_cast<T *>(aligned_alloc(other.alignment_, other.bytes_));
        }

        /// copy the values into the nion
        degree_ = other.degree_; // set the degree
        bytes_ = other.bytes_; // set the number of bytes
        alignment_ = other.alignment_; // set the alignment

        memcpy(components_, other.components_, bytes_);
        return *this; // return the nion
    }

    template<typename T, typename D>
    constexpr nion<T, D> &Nion::nion<T, D>::operator=(const std::initializer_list<T> &vals) {
        /// check if the degree is the same
        if (degree_ != vals.size()) {
            degree_ = vals.size(); // set the degree
            bytes_ = degree_ * sizeof(T); // set the number of bytes
            alignment_ = std::max((bytes_ | (bytes_ - 1)) + 1, std::alignment_of<T>::value); // set the alignment

            /// check if this nion's degree is greater than or equal to the initializer list's size
            if (degree_ < vals.size()) {
                // reallocate memory for the nion
                free(components_);
                components_ = static_cast<T *>(aligned_alloc(alignment_, bytes_));
            }
        }

        /// copy the values into the nion
        memcpy(components_, vals.begin(), bytes_);
        return *this; // return the nion
    }

    template<typename T, typename D>
    constexpr inline nion<T, D> &nion<T, D>::operator=(nion<T, D> &&other)  noexcept {
        /// check if the nions are the same
        if (&other == this) {
            return *this; // return the nion
        }

        /// free the memory of this nion
        free(components_);

        /// move the values into the nion
        degree_ = other.degree_; // set the degree
        bytes_ = other.bytes_; // set the number of bytes
        alignment_ = other.alignment_; // set the alignment
        components_ = other.components_; // set the components

        other.components_ = nullptr; // set the other nion's components to null
        return *this; // return the nion
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline nion<T, D> &nion<T, D>::operator=(S scalar) {
        /// check if the nion is initialized
        if (degree_ <= 0) {
            degree_ = 1; // set the degree
            bytes_ = sizeof(T); // set the number of bytes
            alignment_ = std::max((bytes_ | (bytes_ - 1)) + 1, std::alignment_of<T>::value); // set the alignment

            /// allocate memory for the nion
            components_ = static_cast<T*>(aligned_alloc(alignment_, bytes_));
        }

        zero(); // set the nion to zero
        components_[0] = static_cast<T>(scalar); // set the real component
        return *this; // return the nion
    }

    /************************************
    *  ASSIGNMENT AND ADDITION OPERATORS
    *************************************/

    template<typename T, typename D>
    constexpr inline void nion<T, D>::resize(D newDegree) {
        ASSERT(newDegree >= degree_, "The new degree must be greater than or equal to the current degree.");
        if (newDegree == degree_) return; // return if the degrees are the same

        /// allocate new memory for the nion
        
        /// calculate memory size
        constexpr D alignment_requirement = std::alignment_of<T>::value; // alignment requirement of the type
        D newbytes_ = newDegree * sizeof(T); // set the number of bytes
        
        /// calculate the number of bytes to add to the memory size to make it the LCM of the alignment requirement
        D min_power_of_two = (newbytes_ | (newbytes_ - 1)) + 1; // minimum power of two that is greater than or equal to 'bytes'
        D newalignment_ = std::max(min_power_of_two, alignment_requirement); // minimum alignment size

        /// allocate memory for the nion
        T* newcomponents_ = static_cast<T*>(aligned_alloc(newalignment_, newbytes_));
        
        /// copy the values into the nion
        memcpy(newcomponents_, components_, bytes_);
        
        /// memset the rest of the memory to zero
        memset(newcomponents_ + degree_, 0, newbytes_ - bytes_);
        
        /// free the old memory
        free(components_); // free the memory of the nion
        
        components_ = newcomponents_; // set the components
        degree_ = newDegree; // set the degree
        bytes_ = newbytes_; // set the number of bytes
        alignment_ = newalignment_; // set the alignment
    }
    
    template<typename T, typename D>
    constexpr inline void nion<T, D>::operator+=(const nion<T, D> &other) {
        // if the degree is less than the other nion, resize this nion.
        if (degree_ < other.degree_)
            resize(other.degree_);

        // add the components of the other nion to this nion.
        #pragma vector aligned always
        for (D i = 0; i < other.degree_; i++) {
            components_[i] += other.components_[i];
        }
    }

    template<typename T, typename D>
    constexpr inline void nion<T, D>::operator-=(const nion<T, D> &other) {
        // if the degree is less than the other nion, resize this nion.
        if (degree_ < other.degree_)
            resize(other.degree_);

        // substract the components of the other nion from this nion.
        #pragma vector aligned always
        for (D i = 0; i < other.degree_; i++) {
            components_[i] -= other.components_[i];
        }
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline void nion<T, D>::operator+=(S scalar) const {
        components_[0] += static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline void nion<T, D>::operator-=(S scalar) const {
        components_[0] -= static_cast<T>(scalar);
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
        sum.components_[0] += static_cast<T>(scalar);
        return sum;
    }

    template<typename T, typename D, typename S>
    static constexpr inline nion<T, D> operator+(S scalar, const nion<T, D> &z) {
        return z + static_cast<T>(scalar);
    }
    
    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator++() {
        components_[0]++;
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
        components_[0]--;
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
        #pragma vector aligned always
        for (D i = 0; i < degree_; i++) {
            components_[i] *= static_cast<T>(scalar);
        }
    }
    
    template<typename T, typename D>
    constexpr inline void nion<T, D>::operator/=(const nion<T, D> &other) {
        *this = *this / other;
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline void nion<T, D>::operator/=(S scalar) {
        #pragma vector aligned always
        for (D i = 0; i < degree_; i++) {
            components_[i] /= static_cast<T>(scalar);
        }
    }

    /******************************************
    *        MULTIPLICATION OPERATORS
    *******************************************/

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::conj() const {
        nion<T, D> conjugate = *this;

        // negate all components except the first
        #pragma vector aligned always
        for (D i = 1; i < degree_; i++) {
            conjugate.components_[i] = -conjugate.components_[i];
        }

        return conjugate;
    }

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator-() const {
        nion<T, D> negated = *this;

        #pragma vector aligned always
        for (D i = 0; i < degree_; i++) {
            negated.components_[i] = -negated.components_[i];
        }

        return negated;
    }
    
    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::operator*(const nion<T, D> &other) const {

        switch (degree_) {
            case 1:
                // if the degree is 1, then the product is just the scalar product.
                if (other.degree_ == 1) {
                    nion<T> product = {components_[0] * other.components_[0]};
                    return product;
                }
                else
                    return other * components_[0];
            case 2:
                if (other.degree_ == 2) { // hard-coded complex product
                    nion<T> product =
                    {
                        components_[0] * other.components_[0] - components_[1] * other.components_[1],
                        components_[1] * other.components_[0] + components_[0] * other.components_[1]
                    };
                    return product;
                }
                break;
            case 4:
                if (other.degree_ == 4) { // hard-coded quaternion product
                    nion<T> product =
                    {
                        components_[0] * other.components_[0] - components_[1] * other.components_[1] - components_[2] * other.components_[2] - components_[3] * other.components_[3],
                        components_[1] * other.components_[0] + components_[0] * other.components_[1] - components_[3] * other.components_[2] + components_[2] * other.components_[3],
                        components_[2] * other.components_[0] + components_[3] * other.components_[1] + components_[0] * other.components_[2] - components_[1] * other.components_[3],
                        components_[3] * other.components_[0] - components_[2] * other.components_[1] + components_[1] * other.components_[2] + components_[0] * other.components_[3]
                    };
                    return product;
                }
                break;
            case 8:
                if (other.degree_ == 8) { // hard-coded octonion product ( I know, it's a bit much )
                    nion<T> product =
                    {
                        components_[0] * other.components_[0] - components_[1] * other.components_[1] - components_[2] * other.components_[2] - components_[3] * other.components_[3] - components_[4] * other.components_[4] - components_[5] * other.components_[5] - components_[6] * other.components_[6] - components_[7] * other.components_[7],
                        components_[1] * other.components_[0] + components_[0] * other.components_[1] - components_[3] * other.components_[2] + components_[2] * other.components_[3] - components_[5] * other.components_[4] + components_[4] * other.components_[5] + components_[7] * other.components_[6] - components_[6] * other.components_[7],
                        components_[2] * other.components_[0] + components_[3] * other.components_[1] + components_[0] * other.components_[2] - components_[1] * other.components_[3] - components_[6] * other.components_[4] - components_[7] * other.components_[5] + components_[4] * other.components_[6] + components_[5] * other.components_[7],
                        components_[3] * other.components_[0] - components_[2] * other.components_[1] + components_[1] * other.components_[2] + components_[0] * other.components_[3] - components_[7] * other.components_[4] + components_[6] * other.components_[5] - components_[5] * other.components_[6] + components_[4] * other.components_[7],
                        components_[4] * other.components_[0] + components_[5] * other.components_[1] + components_[6] * other.components_[2] + components_[7] * other.components_[3] + components_[0] * other.components_[4] - components_[1] * other.components_[5] - components_[2] * other.components_[6] - components_[3] * other.components_[7],
                        components_[5] * other.components_[0] - components_[4] * other.components_[1] + components_[7] * other.components_[2] - components_[6] * other.components_[3] + components_[1] * other.components_[4] + components_[0] * other.components_[5] + components_[3] * other.components_[6] - components_[2] * other.components_[7],
                        components_[6] * other.components_[0] - components_[7] * other.components_[1] - components_[4] * other.components_[2] + components_[5] * other.components_[3] + components_[2] * other.components_[4] - components_[3] * other.components_[5] + components_[0] * other.components_[6] + components_[1] * other.components_[7],
                        components_[7] * other.components_[0] + components_[6] * other.components_[1] - components_[5] * other.components_[2] - components_[4] * other.components_[3] + components_[3] * other.components_[4] + components_[2] * other.components_[5] - components_[1] * other.components_[6] + components_[0] * other.components_[7]
                    };
                    return product;
                }
                break;
            case 16:
                if (other.degree_ == 16) { // hard-coded sedenion product ( I strongly recommend using nowrap if you're looking at this... )
                    nion<T> product =
                    {
                        components_[0] * other.components_[0]  - components_[1] * other.components_[1]  - components_[2] * other.components_[2]  - components_[3] * other.components_[3]  - components_[4] * other.components_[4]  - components_[5] * other.components_[5]  - components_[6] * other.components_[6]  - components_[7] * other.components_[7]  - components_[8] * other.components_[8]  - components_[9] * other.components_[9]  - components_[10] * other.components_[10] - components_[11] * other.components_[11] - components_[12] * other.components_[12] - components_[13] * other.components_[13] - components_[14] * other.components_[14] - components_[15] * other.components_[15],
                        components_[1] * other.components_[0]  + components_[0] * other.components_[1]  - components_[3] * other.components_[2]  + components_[2] * other.components_[3]  - components_[5] * other.components_[4]  + components_[4] * other.components_[5]  - components_[7] * other.components_[6]  + components_[6] * other.components_[7]  - components_[9] * other.components_[8]  + components_[8] * other.components_[9]  - components_[11] * other.components_[10] + components_[10] * other.components_[11] - components_[13] * other.components_[12] + components_[12] * other.components_[13] - components_[15] * other.components_[14] + components_[14] * other.components_[15],
                        components_[0] * other.components_[2]  - components_[1] * other.components_[3]  + components_[2] * other.components_[0]  + components_[3] * other.components_[1]  + components_[4] * other.components_[6]  + components_[5] * other.components_[7]  - components_[6] * other.components_[4]  - components_[7] * other.components_[5]  + components_[8] * other.components_[10] + components_[9] * other.components_[11] - components_[10] * other.components_[8]  - components_[11] * other.components_[9]  - components_[12] * other.components_[14] - components_[13] * other.components_[15] + components_[14] * other.components_[12] + components_[15] * other.components_[13],
                        components_[0] * other.components_[3]  + components_[1] * other.components_[2]  - components_[2] * other.components_[1]  + components_[3] * other.components_[0]  + components_[4] * other.components_[7]  - components_[5] * other.components_[6]  + components_[6] * other.components_[5]  - components_[7] * other.components_[4]  + components_[8] * other.components_[11] - components_[9] * other.components_[10] + components_[10] * other.components_[9]  - components_[11] * other.components_[8]  - components_[12] * other.components_[15] + components_[13] * other.components_[14] - components_[14] * other.components_[13] + components_[15] * other.components_[12],
                        components_[0] * other.components_[4]  - components_[1] * other.components_[5]  - components_[2] * other.components_[6]  - components_[3] * other.components_[7]  + components_[4] * other.components_[0]  + components_[5] * other.components_[1]  + components_[6] * other.components_[2]  + components_[7] * other.components_[3]  + components_[8] * other.components_[12] + components_[9] * other.components_[13] + components_[10] * other.components_[14] + components_[11] * other.components_[15] - components_[12] * other.components_[8]  - components_[13] * other.components_[9]  - components_[14] * other.components_[10] - components_[15] * other.components_[11],
                        components_[0] * other.components_[5]  + components_[1] * other.components_[4]  - components_[2] * other.components_[7]  + components_[3] * other.components_[6]  - components_[4] * other.components_[1]  + components_[5] * other.components_[0]  - components_[6] * other.components_[3]  + components_[7] * other.components_[2]  + components_[8] * other.components_[13] - components_[9] * other.components_[12] + components_[10] * other.components_[15] - components_[11] * other.components_[14] + components_[12] * other.components_[9]  - components_[13] * other.components_[8]  + components_[14] * other.components_[11] - components_[15] * other.components_[10],
                        components_[0] * other.components_[6]  + components_[1] * other.components_[7]  + components_[2] * other.components_[4]  - components_[3] * other.components_[5]  - components_[4] * other.components_[2]  + components_[5] * other.components_[3]  + components_[6] * other.components_[0]  - components_[7] * other.components_[1]  + components_[8] * other.components_[14] - components_[9] * other.components_[15] - components_[10] * other.components_[12] + components_[11] * other.components_[13] + components_[12] * other.components_[10] - components_[13] * other.components_[11] - components_[14] * other.components_[8]  + components_[15] * other.components_[9],
                        components_[0] * other.components_[7]  - components_[1] * other.components_[6]  + components_[2] * other.components_[5]  + components_[3] * other.components_[4]  - components_[4] * other.components_[3]  - components_[5] * other.components_[2]  + components_[6] * other.components_[1]  + components_[7] * other.components_[0]  + components_[8] * other.components_[15] + components_[9] * other.components_[14] - components_[10] * other.components_[13] - components_[11] * other.components_[12] + components_[12] * other.components_[11] + components_[13] * other.components_[10] - components_[14] * other.components_[9]  - components_[15] * other.components_[8],
                        components_[0] * other.components_[8]  - components_[1] * other.components_[9]  - components_[2] * other.components_[10] - components_[3] * other.components_[11] - components_[4] * other.components_[12] - components_[5] * other.components_[13] - components_[6] * other.components_[14] - components_[7] * other.components_[15] + components_[8] * other.components_[0]  + components_[9] * other.components_[1]  + components_[10] * other.components_[2]  + components_[11] * other.components_[3]  + components_[12] * other.components_[4]  + components_[13] * other.components_[5]  + components_[14] * other.components_[6]  + components_[15] * other.components_[7],
                        components_[0] * other.components_[9]  + components_[1] * other.components_[8]  - components_[2] * other.components_[11] + components_[3] * other.components_[10] - components_[4] * other.components_[13] + components_[5] * other.components_[12] + components_[6] * other.components_[15] - components_[7] * other.components_[14] - components_[8] * other.components_[1]  + components_[9] * other.components_[0]  - components_[10] * other.components_[3]  + components_[11] * other.components_[2]  - components_[12] * other.components_[5]  + components_[13] * other.components_[4]  + components_[14] * other.components_[7]  - components_[15] * other.components_[6],
                        components_[0] * other.components_[10] + components_[1] * other.components_[11] + components_[2] * other.components_[8]  - components_[3] * other.components_[9]  - components_[4] * other.components_[14] - components_[5] * other.components_[15] + components_[6] * other.components_[12] + components_[7] * other.components_[13] - components_[8] * other.components_[2]  + components_[9] * other.components_[3]  + components_[10] * other.components_[0]  - components_[11] * other.components_[1]  - components_[12] * other.components_[6]  - components_[13] * other.components_[7]  + components_[14] * other.components_[4]  + components_[15] * other.components_[5],
                        components_[0] * other.components_[11] - components_[1] * other.components_[10] + components_[2] * other.components_[9]  + components_[3] * other.components_[8]  - components_[4] * other.components_[15] + components_[5] * other.components_[14] - components_[6] * other.components_[13] + components_[7] * other.components_[12] - components_[8] * other.components_[3]  - components_[9] * other.components_[2]  + components_[10] * other.components_[1]  + components_[11] * other.components_[0]  - components_[12] * other.components_[7]  + components_[13] * other.components_[6]  - components_[14] * other.components_[5]  + components_[15] * other.components_[4],
                        components_[0] * other.components_[12] + components_[1] * other.components_[13] + components_[2] * other.components_[14] + components_[3] * other.components_[15] + components_[4] * other.components_[8]  - components_[5] * other.components_[9]  - components_[6] * other.components_[10] - components_[7] * other.components_[11] - components_[8] * other.components_[4]  + components_[9] * other.components_[5]  + components_[10] * other.components_[6]  + components_[11] * other.components_[7]  + components_[12] * other.components_[0]  - components_[13] * other.components_[1]  - components_[14] * other.components_[2]  - components_[15] * other.components_[3],
                        components_[0] * other.components_[13] - components_[1] * other.components_[12] + components_[2] * other.components_[15] - components_[3] * other.components_[14] + components_[4] * other.components_[9]  + components_[5] * other.components_[8]  + components_[6] * other.components_[11] - components_[7] * other.components_[10] - components_[8] * other.components_[5]  - components_[9] * other.components_[4]  + components_[10] * other.components_[7]  - components_[11] * other.components_[6]  + components_[12] * other.components_[1]  + components_[13] * other.components_[0]  + components_[14] * other.components_[3]  - components_[15] * other.components_[2],
                        components_[0] * other.components_[14] - components_[1] * other.components_[15] - components_[2] * other.components_[12] + components_[3] * other.components_[13] + components_[4] * other.components_[10] - components_[5] * other.components_[11] + components_[6] * other.components_[8]  + components_[7] * other.components_[9]  - components_[8] * other.components_[6]  - components_[9] * other.components_[7]  - components_[10] * other.components_[4]  + components_[11] * other.components_[5]  + components_[12] * other.components_[2]  - components_[13] * other.components_[3]  + components_[14] * other.components_[0]  + components_[15] * other.components_[1],
                        components_[0] * other.components_[15] + components_[1] * other.components_[14] - components_[2] * other.components_[13] - components_[3] * other.components_[12] + components_[4] * other.components_[11] + components_[5] * other.components_[10] - components_[6] * other.components_[9]  + components_[7] * other.components_[8]  - components_[8] * other.components_[7]  + components_[9] * other.components_[6]  - components_[10] * other.components_[5]  - components_[11] * other.components_[4]  + components_[12] * other.components_[3]  + components_[13] * other.components_[2]  - components_[14] * other.components_[1]  + components_[15] * other.components_[0]
                     };

                    return product;
                }
                break;
            //case 32: TODO: No way. You can't make me. 32x32 products are just too big. I could do it, but it wouldn't be pretty. I'll leave it to you.
        }

        // if the other degree is 1, then the product is just the scalar product.
        if (other.degree_ == 1)
            return *this * other.components_[0];


        ///*** if the degree is not 1, 2, 4, or 8, then we need to do the recursive product ***

        // the degree of the first half of the nion
        std::size_t this_half_degree = degree_ - degree_ / 2;
        std::size_t other_half_degree = other.degree_ - other.degree_ / 2;

        nion<T> a, b, c, d; // the four halves of the nions

        a.degree_ = this_half_degree; // the degree of the first half of the nion
        b.degree_ = degree_ - this_half_degree; // the degree of the second half of the nion (can be one less than the first half)
        c.degree_ = other_half_degree; // same as a for the other nion
        d.degree_ = other.degree_ - other_half_degree; // same as b for the other nion


        /// copy the references to the parts of the nions (this and other) into the parts of the nions (a, b) * (c, d)
        // this is a bit tricky, but it's the only way to do it without copying the data

        a.components_ = components_; // a is the first pair of this
        b.components_ = components_ + this_half_degree; // b is the second pair of this
        c.components_ = other.components_; // c is the first pair of other
        d.components_ = other.components_ + other_half_degree; // d is the second pair of other

        /// compute new bytes and alignments for each of the parts
        a.bytes_ = a.degree_ * sizeof(T);
        b.bytes_ = b.degree_ * sizeof(T);
        c.bytes_ = c.degree_ * sizeof(T);
        d.bytes_ = d.degree_ * sizeof(T);

        constexpr D alignment_requirement = std::alignment_of<T>::value; // alignment requirement of the type T
        a.alignment_ = std::max((a.bytes_ | (a.bytes_ - 1)) + 1, alignment_requirement);
        b.alignment_ = std::max((b.bytes_ | (b.bytes_ - 1)) + 1, alignment_requirement);
        c.alignment_ = std::max((c.bytes_ | (c.bytes_ - 1)) + 1, alignment_requirement);
        d.alignment_ = std::max((d.bytes_ | (d.bytes_ - 1)) + 1, alignment_requirement);

        /// calculate the cayley-dickson product
        nion<T> product(
                (a * c) - (d.conj() * b), // maybe add involution parameter for sign (split hypercomplex numbers)
                (d * a) + (b * c.conj())
        );

        // necessary to prevent the nion pairs from deleting the data
        a.components_ = b.components_ = c.components_ = d.components_ = nullptr;
        return product;

    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline nion<T, D> nion<T, D>::operator*(S scalar) const {
        nion<T, D> product = *this;

        #pragma vector aligned always
        for (D i = 0; i < degree_; i++) {
            product.components_[i] = components_[i] * static_cast<T>(scalar);
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

        #pragma vector aligned always
        for (D i = 0; i < degree_; i++) {
            absVal += components_[i] * components_[i];
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


        T absolute = abs();
        if (absolute < epsilon) {
            // if the absolute value is zero, then use the product definition of the absolute value.
            // zero divisors are possible in nions with degree >= 16, so we need to check for them.
            // this is a bit of a hack, but it works. (like the rest of this library)

            // q* / |q| = q* / sqrt(q * q*)
            nion<T, D> inverse = this->conj() * pow((*this * this->conj()), static_cast<T>(-0.5)).real();
            return inverse;
        }

        nion<T, D> inverse = *this;
        inverse.components_[0] /= absolute;
        #pragma vector aligned always
        for (D i = 1; i < degree_; i++) {
            inverse.components_[i] /= -absolute;
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
        return this->components_[index];
    }
    
    template<typename T, typename D>
    constexpr inline T nion<T, D>::real() const { return components_[0]; }

    template<typename T, typename D>
    constexpr inline nion<T, D> nion<T, D>::imag() const {
        nion<T> imag = *this;
        imag[0] = 0;
        return imag;
    }

    /******************************************
    *            COMPARISON OPERATORS
    *******************************************/

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::operator==(const nion<T, D> &other) const {
        if (degree_ != other.degree_) {
            return false;
        }
        #pragma vector aligned always
        for (D i = 0; i < degree_; i++) {
            if (!value_is_similar(components_[i], other[i])) {
                return false;
            }
        }
        return true;
    }

    template<typename T, typename D>
    constexpr inline bool nion<T, D>::operator!=(const nion<T, D> &other) const {
        if (degree_ != other.degree_) {
            return true;
        }
        #pragma vector aligned always
        for (D i = 0; i < degree_; i++) {
            if (!value_is_similar(components_[i], other[i])) {
                return true;
            }
        }
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
        #pragma vector aligned always
        for (D i = 1; i < degree_; i++) {
            if (!value_is_similar(components_[i], 0)) {
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
        return rotate_real() > static_cast<T>(scalar);
    }

    template<typename T, typename D, typename S>
    static constexpr inline bool operator>(S scalar, const nion<T, D> &z) {
        return z < static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline bool nion<T,D>::operator<(S scalar) const{
        return rotate_real() < static_cast<T>(scalar);
    }

    template<typename T, typename D, typename S>
    static constexpr inline bool operator<(S scalar, const nion<T, D> &z) {
        return z > static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline bool nion<T,D>::operator>=(S scalar) const{
        if (*this == nion<T, D>(static_cast<T>(scalar), this->degree_))
            return true;
        return rotate_real() > static_cast<T>(scalar);
    }

    template<typename T, typename D, typename S>
    static constexpr inline bool operator>=(S scalar, const nion<T, D> &z) {
        return z <= static_cast<T>(scalar);
    }

    template<typename T, typename D>
    template<typename S>
    constexpr inline bool nion<T,D>::operator<=(S scalar) const{
        if (*this == nion<T, D>(static_cast<T>(scalar), this->degree_))
            return true;
        return rotate_real() < static_cast<T>(scalar);
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
        T component = z.components_[0];
        os << "(" << component;

        for (D i = 1; i < z.degree_; i++) {
            component = z.components_[i];
            os << "," << component;
        }
        os << ")";
        return os;
    }

    template<typename T, typename D>
    static std::istream &operator>>(std::istream &is, nion<T, D> &z) {
        for (D i = 0; i < z.degree_; i++) {
            is >> z.components_[i];
        }
        return is;
    }

    /**
    * @brief Converts a nion to a string.
    * @return The string representation of the nion.
    */
    template<typename T, typename D>
    inline std::string nion<T, D>::to_string() const {
        std::string nion_string = "(" + std::to_string(components_[0]);
        for (D i = 1; i < degree_; i++) {
            nion_string += "," + std::to_string(components_[i]);
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
        D minDegree = std::min(lhs.degree_, rhs.degree_);
        #pragma vector aligned always
        for (D i = 0; i < minDegree; i++) {
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
    static constexpr inline nion<T, D> exp(const nion<T, D> &z) {

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
            return nion<T, D>(exp_r, z.degree_);

        // compute exponential of nion
        cos_theta = std::cos(i_norm);
        sin_theta = std::sin(i_norm);

        return i*(exp_r * sin_theta / i_norm) + exp_r * cos_theta;

    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> log(const nion<T, D> &z) {

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
    static constexpr inline nion<T, D> pow(const nion<T, D> &base, S power) {

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
                return nion<T, D>(std::pow(base_norm, static_cast<T>(power)), base.degree_); // if base is positive
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
                    return nion<T, D>(1, z.degree_);
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
    static constexpr inline nion<T, D> pow(const nion<T, D> &base, const nion<T, D> &power) {
        return exp(power * log(base));
    }

    template<typename T, typename D, typename S>
    static constexpr inline nion<T, D> pow(S base, const nion<T, D> &power) {
        return pow(nion<T, D>(static_cast<T>(base), power.degree_), power);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> sqr(const nion<T, D> &base) {
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
            return nion<T>(base_norm * base_norm, base.degree_); // if base is zero return zero (0^2 = 0)

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
    static constexpr inline nion<T, D> sqrt(const nion<T, D> &z) {
        return pow(z, static_cast<T>(0.5));
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> cbrt(const nion<T, D> &z) {
        return pow(z, static_cast<T>(1.0) / static_cast<T>(3.0));
    }

    /*******************************************
    *  NION HYPERBOLIC TRIGONOMETRIC FUNCTIONS *
    ********************************************/

    template<typename T, typename D>
    static constexpr inline nion<T, D> sinh(const nion<T, D> &z) {
        // get polar form of nion
        nion<T, D> i = z.imag();

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);

        // calculate scalars
        T e_z = std::exp(z[0]) / static_cast<T>(2);
        T e_mz = std::exp(-z[0]) / static_cast<T>(2);

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
    static constexpr inline nion<T, D> cosh(const nion<T, D> &z) {
        // get polar form of nion
        nion<T, D> i = z.imag();

        // make unit vector
        T i_abs = i.abs();
        T i_norm = std::sqrt(i_abs);

        // calculate scalars
        T e_z = std::exp(z[0]) / static_cast<T>(2);
        T e_mz = std::exp(-z[0]) / static_cast<T>(2);

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
    static constexpr inline nion<T, D> tanh(const nion<T, D> &z) {
        return (exp(z*2) - 1) / (exp(z*2) + 1);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> coth(const nion<T, D> &z) {
        return tanh(z).inv();
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> sech(const nion<T, D> &z) {
        return cosh(z).inv();
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> csch(const nion<T, D> &z) {
        return sinh(z).inv();
    }

    /********************************
    *  NION TRIGONOMETRIC FUNCTIONS *
    *********************************/


    template<typename T, typename D>
    static constexpr inline nion<T, D> sin(const nion<T, D> &z) {
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
    static constexpr inline nion<T, D> cos(const nion<T, D> &z) {
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
    static constexpr inline nion<T, D> tan(const nion<T, D> &z) {
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
    static constexpr inline nion<T, D> cot(const nion<T, D> &z) {
        return tan(z).inv();
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> sec(const nion<T, D> &z) {
        return cos(z).inv();
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> csc(const nion<T, D> &z) {
        return sin(z).inv();
    }

    /***************************************************
    *  NION INVERSE HYPERBOLIC TRIGONOMETRIC FUNCTIONS *
    ****************************************************/

    template<typename T, typename D>
    static constexpr inline nion<T, D> asinh(const nion<T, D> &z) {
        return log(z + sqrt(sqr(z) + static_cast<T>(1)));
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> acosh(const nion<T, D> &z) {
        // compute the inverse hyperbolic cosine of the nion
        return 2 * log(sqrt((z-1) * static_cast<T>(0.5)) + sqrt((z+1) * static_cast<T>(0.5)));

    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> atanh(const nion<T, D> &z) {
        return (log(static_cast<T>(1) + z) - log(static_cast<T>(1) - z)) * static_cast<T>(0.5);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> acoth(const nion<T, D> &z) {
        return (log(static_cast<T>(1) + inv(z)) - log(static_cast<T>(1) - inv(z))) * static_cast<T>(0.5);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> asech(const nion<T, D> &z) {
        return log(sqrt(1/sqr(z) - static_cast<T>(1)) + inv(z));
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> acsch(const nion<T, D> &z) {
        return log(sqrt(static_cast<T>(1) + 1/sqr(z)) + inv(z));
    }

    /****************************************
    *  NION INVERSE TRIGONOMETRIC FUNCTIONS *
    *****************************************/

    template<typename T, typename D>
    static constexpr inline nion<T, D> asin(const nion<T, D> &z) {

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
    static constexpr inline nion<T, D> acos(const nion<T, D> &z) {
        return static_cast<T>(M_PI_2l) - asin(z);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> atan(const nion<T, D> &z) {
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
    static constexpr inline nion<T, D> acot(const nion<T, D> &z) {
        return static_cast<T>(M_PI_2l) - atan(z);
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> asec(const nion<T, D> &z) {
        return acos(inv(z));
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> acsc(const nion<T, D> &z) {
        return asin(inv(z));
    }

    template<typename T, typename D>
    static constexpr inline nion<T, D> atan2(const nion<T, D> &y, const nion<T, D> &x) {
        return atan(y / x);
    }

    /***************************
     *   NION GAMMA FUNCTION   *
     **************************/

    template<typename T, typename D>
    static constexpr inline nion<T, D> gamma(const nion<T, D> &z) {
        // compute the gamma function of the nion
        return exp(-z) * sqrt(inv(z))
               * pow(static_cast<T>(1) / (static_cast<T>(12) * z - inv(static_cast<T>(10) * z)) + z, z)
               * std::sqrt(static_cast<T>(2) * static_cast<T>(M_PIl));
    }
}
