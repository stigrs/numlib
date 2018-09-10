// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_SLICE_H
#define NUMLIB_SLICE_H

#include <cstddef>
#include <iostream>

namespace Numlib {

// A slice describes a sequence of elements in some dimension (or row) of a
// matrix. It is a triple comprised of a starting index, a number of elements,
// and the stride between subsequent elements.
//
struct slice {
    slice() : start(-1), length(-1), stride(1) {}

    explicit slice(std::ptrdiff_t s) : start(s), length(-1), stride(1) {}

    slice(std::ptrdiff_t s, std::ptrdiff_t l, std::ptrdiff_t n = 1)
        : start(s), length(l), stride(n)
    {
    }

    std::ptrdiff_t operator()(std::ptrdiff_t i) const
    {
        return start + i * stride;
    }

    std::ptrdiff_t start;
    std::ptrdiff_t length;
    std::ptrdiff_t stride;
};

inline std::ostream& operator<<(std::ostream& to, const slice& s)
{
    to << '(' << s.start << ", " << s.length << ", " << s.stride << ')';
    return to;
}

} // namespace Numlib

#endif // NUMLIB_SLICE_H
