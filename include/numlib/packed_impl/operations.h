// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_PACKED_MATRIX_OPERATIONS_H
#define NUMLIB_PACKED_MATRIX_OPERATIONS_H

#include <iomanip>
#include <iostream>

namespace Numlib {

//------------------------------------------------------------------------------
//
// Output to stream:

template <typename T, Uplo_scheme Uplo>
std::ostream& operator<<(std::ostream& to, const Packed_matrix<T, Uplo>& ap)
{
    to << ap.rows() << " x " << ap.cols() << "\n[";
    for (std::ptrdiff_t i = 0; i < ap.rows(); ++i) {
        for (std::ptrdiff_t j = 0; j < ap.cols(); ++j) {
            to << std::setw(9) << ap(i, j) << " ";
        }
        if (i != ap.rows() - 1) {
            to << "\n ";
        }
    }
    to << "]\n";
    return to;
}

} // namespace Numlib

#endif // NUMLIB_PACKED_MATRIX_OPERATIONS_H
