// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_PACKED_MATRIX_TYPE_ALIAS_H
#define NUMLIB_PACKED_MATRIX_TYPE_ALIAS_H

namespace Numlib {

// Provides convenient type aliases:

template <typename T>
using Upper_triang_mat = Packed_matrix<T, upper_triang>;

template <typename T>
using Lower_triang_mat = Packed_matrix<T, lower_triang>;

} // namespace Numlib

#endif // NUMLIB_PACKED_MATRIX_TYPE_ALIAS
