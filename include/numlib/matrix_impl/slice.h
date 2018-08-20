// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_SLICE_H
#define NUMLIB_MATRIX_SLICE_H

#include <cstddef>
#include <iostream>

namespace Numlib {

// A slice describes a sequence of elements in some dimension (or row) of a
// matrix. It is a triple comprised of a starting index, a number of elements,
// and the stride between subsequent elements.
//
// The special member Slice::all represents the selection of all elements
// in a particular dimension.
struct Slice {
    Slice() : start(0), length(0), stride(1) {}

    explicit Slice(std::size_t s) : start(s), length(0), stride(1) {}

    Slice(std::size_t s, std::size_t l, std::size_t n)
        : start(s), length(l), stride(n)
    {
    }

	static Slice all;

	std::size_t start;
	std::size_t length;
	std::size_t stride;
};

std::ostream& operator<<(std::ostream& to, const Slice& s)
{
    to << '(' << s.start << ", " << s.length << ", " << s.stride << ')';
	return to;
}

} // namespace Numlib 

#endif // NUMLIB_MATRIX_SLICE_H
