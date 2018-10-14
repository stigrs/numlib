// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_ARRAY0_H
#define NUMLIB_ARRAY0_H

namespace Numlib {

// The type Array<T, 0> is not really an array. It stores a single scalar
// of type T and can only be converted to a reference to that type.
//
template <class T>
class Array<T, 0> {
public:
    static constexpr std::size_t rank = 0;

    using value_type = T;

    Array() = default;

    // Initialize the pseudo-array.
    Array(const T& x) : elem(x) {}

    // Assignment.
    Array& operator=(const T& value)
    {
        elem = value;
        return *this;
    }

    T& operator()() { return elem; }
    const T& operator()() const { return elem; }

    operator T&() { return elem; }
    operator const T&() { return elem; }

private:
    T elem;
};

} // namespace Numlib

#endif // NUMLIB_ARRAY0_H

