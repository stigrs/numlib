// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_H
#define NUMLIB_MATRIX_H

#include <cstddef>
#include <numlib/traits/traits.h>
#include <numlib/matrix_impl/traits.h>

namespace Numlib {

// Forward declarations:

template <std::size_t N>
struct Matrix_slice;

template <typename T, std::size_t N>
class Matrix;

// Matrix initializer:
template <typename T, std::size_t N>
using Matrix_initializer = typename Matrix_impl::Matrix_init<T, N>::type;

} // namespace Numlib

// Support classes and algorithms:
#include <numlib/matrix_impl/support.h>
#include <numlib/matrix_impl/matrix_slice.h>
#include <numlib/matrix_impl/iterator.h>

// Matrix classes:
#include <numlib/matrix_impl/matrix_ref.h>
#include <numlib/matrix_impl/matrix.h>

// Arithmetic operations:
#include <numlib/matrix_impl/operations.h>

#endif // NUMLIB_MATRIX_H
