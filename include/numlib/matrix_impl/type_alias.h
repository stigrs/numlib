// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_TYPE_ALIAS_H
#define NUMLIB_MATRIX_TYPE_ALIAS_H

#include <complex>

namespace num {

// Provides convenient type aliases:

using vec = Matrix<double, 1>;
using mat = Matrix<double, 2>;
using cube = Matrix<double, 3>;

using dvec = Matrix<double, 1>;
using dmat = Matrix<double, 2>;
using dcube = Matrix<double, 3>;

using zvec = Matrix<std::complex<double>, 1>;
using zmat = Matrix<std::complex<double>, 2>;
using zcube = Matrix<std::complex<double>, 3>;

using ivec = Matrix<int, 1>;
using imat = Matrix<int, 2>;
using icube = Matrix<int, 3>;

using uvec = Matrix<unsigned, 1>;
using umat = Matrix<unsigned, 2>;
using ucube = Matrix<unsigned, 3>;

} // namespace num

#endif // NUMLIB_MATRIX_TYPE_ALIAS
