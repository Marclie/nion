/*
 Copyright 2023 Marcus Dante Liebenthal

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

#ifndef NIOLLOC_HPP
#define NIOLLOC_HPP

#ifndef ASSERT
#define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << message << ": file=" << __FILE__ \
                      << ", line=" << __LINE__ << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#endif

#include <cmath>
#include <limits>
#include <type_traits>
#include <iostream>
#include <stdexcept>
#include <cstring>
#include <concepts>
#include <vector>
#include <array>

// macro to determine if heap should be used for memory (zero enables heap)
#define NION_USE_HEAP 0

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
 * @file niolloc.hpp
 * @brief Templated class that manages memory for Cayley-Dickson construction
 * @author Marcus Dante Liebenthal
 * @version 1.0
 * @date 2022-10-08
 */
template<arith_ops T, // type of the coefficients (required to have arithmetic operations)
        std::size_t N = NION_USE_HEAP> // type of the size
// (default is zero, where heap is used for memory instead of stack)
struct niolloc; // forward declaration

// create a constexpr function that counts the number of bits in an integer
template <typename T>
constexpr std::size_t count_size_bits(T n) {
    std::size_t bits = 0; // initialize bits to zero
    while (n > 0) { // while `n` still has bits
        n >>= 1; // remove the least significant bit
        ++bits; // increment the number of bits
    } return bits; // return the number of bits
}

template<arith_ops T, std::size_t N>
struct niolloc {

    // determine the minimum number of bits required to represent N at compile time
    static constexpr std::size_t N_BITS = count_size_bits(N);
    static constexpr bool on_heap = N == NION_USE_HEAP; // determine if heap should be used for memory
    static constexpr bool is_trivial = std::is_trivial_v<T>; // determine if T is a trivial type

    // set the internal integer type to the smallest width that can hold N
    using D = std::conditional_t<on_heap, std::size_t, // default to std::size_t if N is 0 for heap allocation
            std::conditional_t<N_BITS <= 8, uint8_t,     // max size of 256
                    std::conditional_t<N_BITS <= 16, uint16_t,   // max size of 65536
                            std::conditional_t<N_BITS <= 32, uint32_t,   // max size of 4294967296
                                    std::size_t>>>>; // max size of 18446744073709551616

    /// coefficients
    // if N is NION_USE_HEAP, then use heap for memory, else use stack
    using elem_type = std::conditional_t<on_heap, T*, T[N]>; // container for coefficients

    elem_type vals_{}; // declare an array of coefficients
    D size_; // number of coefficients

    /**
     * @brief free the memory on the heap
     * @note this function is a no-op if the memory is on the stack
     */
    constexpr inline void destroy() {
        if constexpr (on_heap) {
            if constexpr (is_trivial) {
                free(vals_); // free memory on heap (if trivially constructible)
            } else {
                delete[] vals_; // free memory on heap (if not trivially constructible)
            }
            vals_ = nullptr; // set pointer to null
        }
    }

    /**
     * @brief allocate memory on the heap
     * @param size number of elements to allocate
     * @note this function is a no-op if the memory is on the stack
     */
    template<integral_type Z>
    constexpr inline void create(Z size) {
        if constexpr (on_heap) {
            if constexpr (is_trivial) {
                if (vals_ != nullptr)
                    throw std::runtime_error("Cannot create memory for nion on heap if it already exists");
                vals_ = (T *) std::malloc(size * sizeof(T)); // allocate memory on heap
            } else {
                vals_ = new T[size]; // allocate memory on heap
            }
        }
    }

    /**
     * @brief realloc memory on the heap to larger size
     * @param new_size new number of elements to allocate
     * @note this function is a no-op if the memory is on the stack
     */
    template<integral_type Z>
    constexpr inline void expand(Z size) {
        if (size_ >= size) {size_ = size; return;} // if the new size is smaller than the old size, only change size
        if constexpr (on_heap) {
            if (vals_ == nullptr) {
                create(size); return; // if the memory is not allocated, allocate it and return
            }
            if constexpr (is_trivial) { // allocate memory on heap (if trivially constructible)

                auto vals_new = (T*) std::realloc(vals_, size * sizeof(T)); // reallocate memory on heap
                if (vals_new == nullptr)
                    throw std::runtime_error("realloc for nion failed"); // throw error if realloc failed

                vals_ = vals_new; // set pointer to new memory

            } else { // create a new array and copy old values
                destroy(); // free memory on heap
                create(size); // allocate new memory on heap
            }
        }
        size_ = size; // set size to new size
    }

    /**
     * @brief copy memory from source values to this niolloc
     * @param src source values
     * @param count number of elements to copy
     * @param this_off offset for values in this niolloc (default 0)
     * @param src_off offset in the source niolloc (default 0)
     */
    template<arith_ops S,
            integral_type Z1 = int, integral_type Z2 = int, integral_type Z3 = int> // types of the offsets
    requires (std::is_convertible_v<S,T>)
    constexpr inline void copy(const S* src, Z1 count, Z2 this_off = 0, Z3 src_off = 0) {
        if (count <= 0)
            throw std::runtime_error("count must be greater than zero"); // throw error if count is less than zero
        // just use memcpy if the element type is trivially copyable, and the types are the same
        if constexpr (is_trivial && std::is_same_v<S,T>) {
            std::memcpy(vals_ + this_off, src + src_off, count * sizeof(S));
        } else { // else, copy each value one by one
            for (D i = 0; i < count; ++i) {
                vals_[this_off + i] = static_cast<T>(src[src_off + i]);
            }
        }
    }

    /**
     * @brief copy memory from one niolloc to another
     * @param dest destination niolloc
     * @param src source niolloc
     * @param count number of elements to copy
     * @param dest_off offset in the destination niolloc (default 0)
     * @param src_off offset in the source niolloc (default 0)
     * @note combinations of types and sizes are checked at compile time
     */
    template<arith_ops S, arith_ops U, // types of the destination and source niollocs
            std::size_t M, std::size_t P, // sizes of the destination and source niollocs
            integral_type Z1 = int, integral_type Z2 = int, integral_type Z3 = int> // types of the offsets
    requires (std::is_convertible_v<S,U>)
    constexpr inline static void nioncpy(niolloc<S,M> &dest, const niolloc<U,P> &src, Z1 count,
                                      Z2 dest_off = 0, Z3 src_off = 0) {
//        ASSERT(dest_off + count <= dest.size_, "niolloc::nioncpy: destination offset out of bounds");
//        ASSERT(src_off + count <= src.size_, "niolloc::nioncpy: source offset out of bounds");
        dest.copy(src.vals_, count, dest_off, src_off); // call the copy function for the destination niolloc
    }

    /**
     * @brief copy memory from source niolloc to this niolloc
     * @param src source values
     * @param count number of elements to copy
     * @param this_off offset for values in this niolloc (default 0)
     * @param src_off offset in the source niolloc (default 0)
     */
    template<arith_ops S, std::size_t M,
            integral_type Z1 = int, integral_type Z2 = int, integral_type Z3 = int> // types of the offsets
    requires (std::is_convertible_v<S,T>)
    constexpr inline void copy(const niolloc<S,M> &src, Z1 count, Z2 this_off = 0, Z3 src_off = 0) {
        nioncpy(*this, src, count, this_off, src_off); // call the static copy function
    }

    /**
     * @brief move memory from source values to this niolloc
     * @param src source values
     * @note src is assumed to have the same size as this niolloc's size
     */
    template<arith_ops S, std::size_t M> // types of the offsets and count
    requires (std::is_convertible_v<S,T>)
    constexpr inline void steal(niolloc<S,M> &src) {
        // if the user wants to use the heap, the source is on the heap, and the types are the same:
        if constexpr (on_heap && niolloc<S,M>::on_heap && std::is_same_v<S,T>) {
            if (src.vals_ != vals_) {
                if (vals_ != nullptr) destroy(); // free memory on heap
                vals_ = src.vals_; // steal the pointer
                src.vals_ = nullptr; // set the other pointer to null
            }
        } else { // else, the user cannot steal the pointer
            copy(src.vals_, size_); // copy the values into the nion
        }
    }

    // default constructor
    constexpr inline niolloc<T,N>() : size_(1) {
        create(size_); // allocate memory on heap
    }

    // destructor
    constexpr inline ~niolloc() {
        destroy(); // free memory on heap
    }

    /**
     * @brief fill the nion with a given value
     * @tparam S a type that is convertible or equal to T
     * @param val value to fill the nion with
     * @param this_off offset for values in this niolloc (default 0)
     */
    template <arith_ops S> requires (std::is_convertible_v<S, T>)
    constexpr inline void fill(const S& value, D start = 0, D end = 0) {
        if (end == 0) end = size_; // if 'end' is zero, set it to the size of the nion
        if constexpr (is_trivial && std::is_same_v<S,T>) {
            std::memset(vals_ + start, value, (end - start) * sizeof(T)); // fill the nion with the value
        } else {
            for (D i = start; i < end; ++i) {
                vals_[i] = static_cast<T>(value); // fill the nion with the value
            }
        }
    }

    /**
     * @brief construct an empty niolloc object with a given size
     * @tparam Z integral type
     * @param size size of the nion
     */
    template<integral_type Z>
    constexpr inline explicit niolloc<T,N>(Z size) : size_(size) {
        // check if the degree is greater than zero
        ASSERT(size > 0, "The size of the nion must be greater than zero.");
        ASSERT(size_ <= N || on_heap,
               "The size of the nion is too large. "
               "consider increasing the template parameter, N, or using the default value for heap allocation."
        );

        // allocate memory
        create(size_);
    }

    /**
     * @brief construct single element in the niolloc
     * @tparam S type of the element to construct (cannot be integral to avoid ambiguity)
     */
    template<arith_ops S> requires (!std::is_integral_v<S> && std::is_convertible_v<S,T>)
    constexpr inline explicit niolloc<T,N>(const S &value) : size_(1) {
        create(1); // allocate memory for the element
        size_ = 1; // set the size to one (just in case)
        vals_[0] = static_cast<T>(value); // set the value
    }

    /**
     * @brief assign single element in the niolloc
     * @tparam S type of the element to construct (cannot be integral to avoid ambiguity)
     */
    template<arith_ops S> requires (!std::is_integral_v<S> && std::is_convertible_v<S,T>)
    constexpr inline niolloc<T,N> &operator=(const S& value) {
        size_ = 1; // set the size to one
        vals_[0] = static_cast<T>(value); // set the value
        return *this;
    }

    /**
     * @brief construct a niolloc object from an array of values and a size
     * @tparam S a type that is convertible or equal to T
     * @tparam Z integral type
     * @param vals array of values
     * @param size size of the nion
     */
    template <arith_ops S, integral_type Z> requires (std::is_convertible_v<S, T>)
    constexpr inline niolloc<T,N>(const S* vals, Z size) : size_(size) {

        /// check if the degree is greater than zero and less than the maximum size
        ASSERT(size > 0, "The size of the nion must be greater than zero.");
        ASSERT(size_ <= N || on_heap,
               "The size of the nion is too large. "
               "consider increasing the template parameter, N, or using the default value for heap allocation."
        );

        create(size_); // allocate memory on heap
        copy(vals, size); // copy the values into the nion
    }

    /**
     * @brief construct a niolloc object from a stack array of values and a size
     * @tparam S a type that is convertible or equal to T
     * @tparam Z integral type
     * @param vals array of values
     * @param size size of the nion
     */
    template <arith_ops S, integral_type Z> requires (std::is_convertible_v<S, T>)
    constexpr inline niolloc<T,N>(S vals[], Z size) : size_(size) {

        /// check if the degree is greater than zero and less than the maximum size
        ASSERT(size > 0, "The size of the nion must be greater than zero.");
        ASSERT(size_ <= N || on_heap,
               "The size of the nion is too large. "
               "consider increasing the template parameter, N, or using the default value for heap allocation."
        );

        create(size_); // allocate memory on heap
        copy(vals, size); // copy the values into the nion
    }

    /**
     * @brief Construct a new niolloc object from an initializer list
     * @param vals initializer list of values to copy into the nion
     */
    template<arith_ops S> requires (std::is_same_v<T,S> || std::is_convertible_v<S,T>)
    constexpr inline niolloc<T,N>(const std::initializer_list<S> &vals) : size_(vals.size()) {
        /// check if the degree is greater than zero
        ASSERT(size_ > 0, "The size of the nion must be greater than zero.");
        ASSERT(size_ <= N || on_heap,
               "The size of the nion is too large. "
               "consider increasing the template parameter, N, or using the default value for heap allocation."
        );

        create(size_); // allocate memory on heap
        copy(vals.begin(), size_); // copy the values into the nion
    }

    /**
     * @brief Construct a new niolloc object from a vector
     * @param vals vector of values to copy into the nion
     */
    template<arith_ops S> requires (std::is_same_v<T,S> || std::is_convertible_v<S,T>)
    constexpr inline explicit niolloc<T,N>(const std::vector<S> &vals) : size_(vals.size()) {
        /// check if the degree is greater than zero
        ASSERT(size_ > 0, "The size of the nion must be greater than zero.");
        ASSERT(size_ <= N || on_heap,
               "The size of the nion is too large. "
               "consider increasing the template parameter, N, or using the default value for heap allocation."
        );

        create(size_); // allocate memory on heap
        copy(vals.data(), size_); // copy the values into the nion
    }

    /**
     * @brief Copy constructor of a niolloc object
     * @param other niolloc object to copy
     */
    constexpr inline niolloc<T,N>(const niolloc<T,N> &other) : size_(other.size_) {
        ASSERT(size_ <= N || on_heap,
               "The size of the nion is too large. "
               "consider increasing the template parameter, N, or using the default value for heap allocation."
        );

        create(size_); // allocate memory on heap
        copy(other.vals_, size_); // copy the values into the nion
    }

    /**
     * Copy assignment operator of a niolloc object of the same type
     * @param other niolloc object to copy
     * @return reference to the current niolloc object
     */
    constexpr inline niolloc<T,N> &operator=(const niolloc<T,N> &other) {

        if (this == &other) // if the objects are the same
            return *this; // return the current object

        ASSERT(other.size_ <= N || on_heap,
               "The size of the nion is too large. "
               "consider increasing the template parameter, N, or using the default value for heap allocation.");

        expand(other.size_); // reallocate memory if necessary
        size_ = other.size_; // set the new size

        // copy the values into the nion
        copy(other.vals_, size_);

        return *this;
    }

    /**
     * @brief Copy constructor of a new niolloc object from a niolloc object of a different type
     * @tparam S type of the other niolloc object
     * @tparam M memory type and size of the other niolloc object
     * @param other niolloc object to copy
     */
    template<arith_ops S, std::size_t M> requires (
        std::is_convertible_v<S,T> // the types are convertible
        && (!std::is_same_v<T,S> || N != M) // either the types are different or the sizes are different
    )
    constexpr inline explicit niolloc<T,N>(const niolloc<S,M> &other) : size_(other.size_) {
        ASSERT(size_ <= N || on_heap,
               "The size of the other nion is too large. "
               "consider increasing the template parameter, N, or using the default value for heap allocation."
        );

        create(size_); // allocate memory on heap
        copy(other.vals_, size_); // copy the values into the nion
    }

    /**
     * Copy assignment operator of a niolloc object of a different type
     * @tparam S type of the other niolloc object
     * @tparam M memory type and size of the other niolloc object
     * @param other niolloc object to copy
     * @return reference to the current niolloc object
     */
    template<arith_ops S, std::size_t M> requires (
        std::is_convertible_v<S,T> // the types are convertible
        && (!std::is_same_v<T,S> || N != M) // either the types are different or the sizes are different
    )
    constexpr inline niolloc<T,N> &operator=(const niolloc<S,M> &other) {

        if (this == &other) // if the objects are the same
            return *this; // return the current object

        ASSERT(other.size_ <= N || on_heap,
               "The size of the nion is too large. "
               "consider increasing the template parameter, N, or using the default value for heap allocation.");

        expand(other.size_); // reallocate memory if necessary
        size_ = other.size_; // set the new size

        // copy the values into the nion
        copy(other.vals_, size_);

        return *this;
    }

    /**
     * @brief Move constructor of a new niolloc object from a niolloc object of the same type
     * @param other niolloc object to move
     */
    constexpr inline niolloc<T,N>(niolloc<T,N> &&other)  noexcept : size_(other.size_) {
        ASSERT(size_ <= N || on_heap,
               "The size of the nion is too large. "
               "consider increasing the template parameter, N, or using the default value for heap allocation."
        );

        steal(other); // steal the other object's memory
    }

    /**
     * Move assignment operator of a niolloc object of the same type
     * @param other niolloc object to move
     * @return reference to the current niolloc object
     */
    constexpr inline niolloc<T,N> &operator=(niolloc<T,N> &&other) noexcept {

        if (this == &other) // return if the objects are the same
            return *this;

        ASSERT(other.size_ <= N || on_heap,
               "The size of the nion is too large. "
               "consider increasing the template parameter, N, or using the default value for heap allocation.");

        size_ = other.size_; // set the new size
        steal(other); // steal the other object's memory

        return *this;
    }

    /**
     * @brief Move constructor of a new niolloc object from a niolloc object of a different type
     * @tparam M memory type and size of the other niolloc object
     * @param other niolloc object to move
     * @note a shallow copy cannot be performed if the arith_ops types are different
     */
    template<std::size_t M> requires (N != M)
    constexpr inline explicit niolloc<T,N>(niolloc<T,M> &&other) noexcept {
        ASSERT(other.size_ <= N || on_heap,
               "The size of the nion is too large. "
               "consider increasing the template parameter, N, or using the default value for heap allocation.");

        size_ = other.size_; // set the new size
        steal(other); // steal the other object's memory
    }

    /**
     * Move assignment operator of a niolloc object of a different type
     * @tparam M memory type and size of the other niolloc object
     * @param other niolloc object to move
     * @return reference to the current niolloc object
     * @note a shallow copy cannot be performed if the arith_ops types are different
     */
    template<std::size_t M> requires (N != M)
    constexpr inline niolloc<T,N> &operator=(niolloc<T,M> &&other) noexcept {

        ASSERT(other.size_ <= N || on_heap,
               "The size of the nion is too large. "
               "consider increasing the template parameter, N, or using the default value for heap allocation.");

        size_ = other.size_; // set the new size
        steal(other); // steal the other object's memory

        return *this;
    }

    /**
     * @brief Cast operator to a niolloc object of a different type
     * @tparam M memory type and size of the other niolloc object
     * @return a niolloc object of a different type
     * @note this is a deep copy
     */
    template<arith_ops S, std::size_t M>
    constexpr inline explicit operator niolloc<S,M>() const {
        if constexpr (std::is_same_v<T,S> && N == M) { // if the types are the same and the sizes are the same
            return *this; // return the current object
        }

        return niolloc<S,M>(*this); // return the new niolloc object
    }

    /**
     * @brief overload of the subscript operator
     * @param i index of the element to access
     * @tparam Z type of the index
     * @return reference to the element at the given index
     */
    template<integral_type Z>
    constexpr inline T &operator[](Z i) { return vals_[i]; } // we do not check the bounds

    template<integral_type Z>
    constexpr inline const T &operator[](Z i) const { return vals_[i]; } // const overload

    /**
     * @brief define the begin iterator
     * @return pointer to the first element
     */
    constexpr inline T* begin() { return vals_; }
    constexpr inline const T* begin() const { return vals_; } // const overload

    /**
     * @brief define the end iterator
     * @return pointer to one past the last element
     */
    constexpr inline T* end() { return vals_ + size_; }
    constexpr inline const T* end() const { return vals_ + size_; } // const overload

    /**
     * @brief define the next iterator
     * @param current current iterator
     * @return pointer to the next element
     */
    constexpr inline T* next(T* current) {
        T* end_pos = end();
        if (current < end_pos) return ++current;
        else return end_pos;
    }
    constexpr inline const T* next(const T* current) const { // const overload
        T* end_pos = end();
        if (current < end_pos) return ++current;
        else return end_pos;
    }

};


#endif //NIOLLOC_HPP