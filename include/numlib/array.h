// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_ARRAY_H
#define NUMLIB_ARRAY_H

#include <numlib/traits.h>

namespace Numlib {

// N-dimensional dense array class.
//
// Features:
// - Create 0D, 1D, 2D, 3D, and 4D arrays.
// - Continuous storage of elements.
// - Column-major storage order (base 0).
// - Value semantics.
// - Range-checked element access unless NDEBUG is defined.
// - Sub-array views (slicing).
//
// Note:
// - The general Array template exists only to allow specializations.
// - Array indexing uses signed integers.
// - Use e.g. Intel MKL for improved numerical performance.
//
template <class T, std::size_t N>
class Array {
private:
    Array();
};

} // namespace Numlib

#include <numlib/array_impl/array0.h>
#include <numlib/array_impl/array1.h>
#include <numlib/array_impl/array2.h>
#include <numlib/array_impl/array3.h>
#include <numlib/array_impl/array4.h>
#include <numlib/array_impl/operations.h>
#include <numlib/array_impl/type_alias.h>

#endif // NUMLIB_ARRAY_H

