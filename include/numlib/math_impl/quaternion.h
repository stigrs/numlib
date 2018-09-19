// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_QUATERNION_H
#define NUMLIB_MATH_QUATERNION_H

#include <numlib/matrix.h>

namespace Numlib {

// Return rotation matrix for the qiven quaternions.
//
// Notes:
//   Uses the XYZ convention
//
// Algorithm:
//   https://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
//
// Arguments:
//   Vector with quaternions (q0, q1, q2, q3)
//
// Returns:
//   Rotation matrix
//
Mat<double> quat2rotm(const Vec<double>& quat);

} // namespace Numlib

#endif // NUMLIB_MATH_QUATERNION_H
