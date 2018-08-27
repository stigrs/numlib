// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_SLICE_H
#define NUMLIB_MATRIX_SLICE_H

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4245) // signed/unsigned error caused by size_t(-1)
#endif

#include <cstddef>
#include <iostream>

namespace Numlib {

// A slice describes a sequence of elements in some dimension (or row) of a
// matrix. It is a triple comprised of a starting index, a number of elements,
// and the stride between subsequent elements.
//
struct slice {
    slice() : start(-1), length(-1), stride(1) {}

    explicit slice(std::size_t s) : start(s), length(-1), stride(1) {}

    slice(std::size_t s, std::size_t l, std::size_t n = 1)
        : start(s), length(l), stride(n)
    {
    }

    std::size_t operator()(std::size_t i) const { return start + i * stride; }

    std::size_t start;
    std::size_t length;
    std::size_t stride;
};

inline std::ostream& operator<<(std::ostream& to, const slice& s)
{
    to << '(' << s.start << ", " << s.length << ", " << s.stride << ')';
    return to;
}

} // namespace Numlib

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // NUMLIB_MATRIX_SLICE_H
