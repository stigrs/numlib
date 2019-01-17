// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_GEOMETRY_H
#define NUMLIB_MATH_GEOMETRY_H

#include <numlib/constants.h>
#include <numlib/matrix.h>
#include <numlib/math.h>
#include <cmath>

namespace Numlib {

// Convert radians to degrees.
inline double radtodeg(double rad) { return rad * 180.0 / Constants::pi; }

// Convert degrees to radians.
inline double degtorad(double deg) { return deg * Constants::pi / 180.0; }

// Find hypotenuse without overflow or destructive underflow.
double hypot(const double a, const double b);

// Spherical to cartesian.
inline void
sph2cart(double azi, double zen, double r, double& x, double& y, double& z)
{
    x = r * std::cos(azi) * std::sin(zen);
    y = r * std::sin(azi) * std::sin(zen);
    z = r * std::cos(zen);
}

// Cartesian to spherical.
inline void
cart2sph(double x, double y, double z, double& azi, double& zen, double& r)
{
    r = hypot(hypot(x, y), z);
    azi = std::atan2(y, x);
    zen = std::acos(z / r);
}

// Cylindrical to cartesian.
inline void
pol2cart(double th, double r, double zin, double& x, double& y, double& z)
{
    x = r * std::cos(th);
    y = r * std::sin(th);
    z = zin;
}

// Cartesian to cylindrical.
inline void
cart2pol(double x, double y, double z, double& th, double& r, double& zout)
{
    th = std::atan2(y, x);
    r = hypot(x, y);
    zout = z;
}

// Polar to cartesian.
inline void pol2cart(double th, double r, double& x, double& y)
{
    x = r * std::cos(th);
    y = r * std::sin(th);
}

// Cartesian to polar.
inline void cart2pol(double x, double y, double& th, double& r)
{
    th = std::atan2(y, x);
    r = hypot(x, y);
}

// Compute distance between two points.
inline double distance(const Vec<double>& a, const Vec<double>& b)
{
    return norm(b - a);
}

// Compute angle in degrees between three points.
inline double
angle(const Vec<double>& a, const Vec<double>& b, const Vec<double>& c)
{
    const auto& ab = normalize(a - b);
    const auto& bc = normalize(c - b);
    return radtodeg(std::acos(dot(ab, bc)));
}

// Compute dihedral angle in degrees given four points.
double dihedral(const Vec<double>& a,
                const Vec<double>& b,
                const Vec<double>& c,
                const Vec<double>& d);

// Compute the pair-wise distances between observations in n-dim. space.
void pdist_matrix(Mat<double>& dm, const Mat<double>& mat);

// Compute centroid of a XYZ coordinate matrix.
inline Vec<double> centroid(const Mat<double>& xyz)
{
    return sum(xyz, 1) / narrow_cast<double>(xyz.rows());
}

// Perform translation.
void translate(Mat<double>& xyz, double dx, double dy, double dz);

// Perform rotation given a rotation matrix.
void rotate(Mat<double>& xyz, const Mat<double>& rotm);

// Create two-dmensional grid based on coordinates in the vectors x and y.
template <typename T>
void meshgrid(const Vec<T>& x, const Vec<T>& y, Mat<T>& xx, Mat<T>& yy)
{
    Index nc = x.size();
    Index nr = y.size();

    xx.resize(nr, nc);
    yy.resize(nr, nc);

    // xx is a matrix where each row is a copy of x:

    for (Index i = 0; i < nr; ++i) {
        xx.row(i) = x;
    }

    // yy is a matrix where each column is a copy of y:

    for (Index j = 0; j < nc; ++j) {
        yy.column(j) = y;
    }
}

} // namespace Numlib

#endif // NUMLIB_MATH_GEOMETRY_H
