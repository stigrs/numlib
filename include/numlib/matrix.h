// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_H
#define NUMLIB_MATRIX_H

#include <numlib/traits.h>
#include <numlib/matrix_impl/traits.h>

namespace Numlib {

//------------------------------------------------------------------------------
//
// Slice:

struct slice;

//------------------------------------------------------------------------------
//
// Dense N-dimensional matrix class:

template <std::size_t N>
struct Matrix_slice;

template <typename T, std::size_t N>
class Matrix_ref;

template <typename T, std::size_t N>
class Matrix;

// Matrix initializer:
template <typename T, std::size_t N>
using Matrix_initializer = typename Matrix_impl::Matrix_init<T, N>::type;

//------------------------------------------------------------------------------
//
// Band matrix class:

template <typename T>
class Band_matrix;

//------------------------------------------------------------------------------
//
// Packed matrix class:

// Enumeration of triangular storage schemes.
enum Uplo_scheme { upper_triang, lower_triang };

template <typename T, Uplo_scheme Uplo = lower_triang>
class Packed_matrix;

//------------------------------------------------------------------------------
//
// Sparse vector class:

template <typename T>
class Sparse_vector;

} // namespace Numlib

// Support classes and algorithms:
#include <numlib/matrix_impl/slice.h>
#include <numlib/matrix_impl/support.h>
#include <numlib/matrix_impl/matrix_slice.h>
#include <numlib/matrix_impl/iterator.h>

// Matrix classes:
#include <numlib/matrix_impl/matrix_ref.h>
#include <numlib/matrix_impl/matrix.h>
#include <numlib/band_impl/band_matrix.h>
#include <numlib/packed_impl/packed_matrix.h>
#include <numlib/sparse_impl/sparse_vector.h>

// Arithmetic operations:
#include <numlib/matrix_impl/operations.h>
#include <numlib/band_impl/operations.h>
#include <numlib/packed_impl/operations.h>

// Type aliases:
#include <numlib/matrix_impl/type_alias.h>

#endif // NUMLIB_MATRIX_H
