// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_DENSE_MATRIX_TYPE_ALIAS_H
#define NUMLIB_DENSE_MATRIX_TYPE_ALIAS_H

namespace Numlib {

// Provides convenient type aliases:

template <typename T>
using Vec = Matrix<T, 1>;

template <typename T>
using Mat = Matrix<T, 2>;

template <typename T>
using Cube = Matrix<T, 3>;

template <typename T>
using Hypercube = Matrix<T, 4>;

} // namespace Numlib

#endif // NUMLIB_DENSE_MATRIX_TYPE_ALIAS
