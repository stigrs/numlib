// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_EULER_H
#define NUMLIB_MATH_EULER_H

#include <numlib/matrix.h>

namespace Numlib {

// Return rotation matrix around x, y, and z axes.
//
// Convention:
//   Rm = Rm(Z(phi), Y(psi), X(theta)) = ZYX
//   Axis rotation sequence: 3, 2, 1
//
// Algorithm:
//   http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf
//
// Arguments:
//   z: rotation angle in degrees around z-axis (yaw/phi)
//   y: rotation angle in degrees around y-axis (roll/psi)
//   x: rotation angle in radians around x-axis (pitch/theta)
//
// Returns:
//   3x3 rotation matrix giving same rotation as for given angles
//
Mat<double> eul2rotm(double z = 0.0, double y = 0.0, double x = 0.0);

// Return Euler angles from rotation matrix.
//
// Convention:
//   Rm = Rm(Z(phi), Y(psi), X(theta)) = ZYX
//   Axis rotation scheme: 3, 2, 1
//
// Algorithm:
//   http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf
//
// Arguments:
//   rotation matrix
//
// Returns:
//   vector with z, y, x
//   z: rotation angle in degrees around z-axis (yaw/phi)
//   y: rotation angle in degrees around y-axis (roll/psi)
//   x: rotation angle in degrees around x-axis (pitch/theta)
//
Vec<double> rotm2eul(const Mat<double>& rotm);

// Return quaternions from Euler angles.
//
// Convention:
//   Rm = Rm(Z(phi), Y(psi), X(theta)) = ZYX
//   Axis rotation sequence: 3, 2, 1
//
// Algorithm:
//   http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf
//
// Arguments:
//   z: rotation angle in degrees around z-axis (yaw/phi)
//   y: rotation angle in degrees around y-axis (roll/psi)
//   x: rotation angle in degrees around x-axis (pitch/theta)
//
// Returns:
//   array of quaternions [q0, q1, q2, q3]
//
Vec<double> eul2quat(double z = 0.0, double y = 0.0, double x = 0.0);

} // namespace Numlib

#endif // NUMLIB_MATH_EULER_H
