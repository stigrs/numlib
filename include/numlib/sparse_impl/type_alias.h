// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_SPARSE_TYPE_ALIAS_H
#define NUMLIB_SPARSE_TYPE_ALIAS_H

namespace Numlib {

// Provides convenient type aliases:

template <typename T>
using Sp_vec = Sparse_vector<T>;

template <typename T>
using Sp_mat = Sparse_matrix<T>;

} // namespace Numlib

#endif // NUMLIB_SPARSE_TYPE_ALIAS
