// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_PACKED_MATRIX_TYPE_ALIAS_H
#define NUMLIB_PACKED_MATRIX_TYPE_ALIAS_H

namespace Numlib {

// Provides convenient type aliases:

// Symmetric matrix held in packed storage format.
template <typename T, Uplo_scheme Uplo>
using Symm_mat = Packed_matrix<T, Uplo>;

// Packed matrix held in upper triangular format.
template <typename T>
using Upper_triang_mat = Packed_matrix<T, up>;

// Packed matrix held in lower triangular format.
template <typename T>
using Lower_triang_mat = Packed_matrix<T, lo>;

} // namespace Numlib

#endif // NUMLIB_PACKED_MATRIX_TYPE_ALIAS
