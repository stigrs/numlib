// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <numlib/constants.h>
#include <limits>

Numlib::Mat<double> Numlib::eul2rotm(double z, double y, double x)
{
    using namespace Numlib;

    Mat<double> result = identity(3);
    if ((z != 0.0) || (y != 0.0) || (x != 0.0)) {
        const double tol = 2.0 * std::numeric_limits<double>::epsilon();

        z *= Constants::pi / 180.0;
        y *= Constants::pi / 180.0;
        x *= Constants::pi / 180.0;

        double c1 = std::cos(z);
        double c2 = std::cos(y);
        double c3 = std::cos(x);
        double s1 = std::sin(z);
        double s2 = std::sin(y);
        double s3 = std::sin(x);

        double m11 = c1 * c2;
        if (std::abs(m11) < tol) {
            m11 = 0.0;
        }
        double m12 = c1 * s2 * s3 - s1 * c3;
        if (std::abs(m12) < tol) {
            m12 = 0.0;
        }
        double m13 = c1 * s2 * c3 + s1 * s3;
        if (std::abs(m13) < tol) {
            m13 = 0.0;
        }
        double m21 = s1 * c2;
        if (std::abs(m21) < tol) {
            m21 = 0.0;
        }
        double m22 = s1 * s2 * s3 + c1 * c3;
        if (std::abs(m22) < tol) {
            m22 = 0.0;
        }
        double m23 = s1 * s2 * c3 - c1 * s3;
        if (std::abs(m23) < tol) {
            m23 = 0.0;
        }
        double m31 = -s2;
        if (std::abs(m31) < tol) {
            m31 = 0.0;
        }
        double m32 = c2 * s3;
        if (std::abs(m32) < tol) {
            m32 = 0.0;
        }
        double m33 = c2 * c3;
        if (std::abs(m33) < tol) {
            m33 = 0.0;
        }
        result = {{m11, m12, m13}, {m21, m22, m23}, {m31, m32, m33}};
    }
    return result;
}

Numlib::Vec<double> Numlib::rotm2eul(const Numlib::Mat<double>& rotm)
{
    using namespace Numlib;

    const double tol = 2.0 * std::numeric_limits<double>::epsilon();

    double m11 = rotm(0, 0);
    double m21 = rotm(1, 0);
    double m31 = rotm(2, 0);
    double m32 = rotm(2, 1);
    double m33 = rotm(2, 2);

    double z = 0.0;
    double y = 0.0;
    double x = 0.0;

    if ((std::abs(m11) < tol) && (std::abs(m21) < tol)) {
        z = 0.0;
    }
    else if ((std::abs(m11) < tol) && (std::abs(m21) > tol)) {
        z = sign(0.5, m21) * Constants::pi;
    }
    else {
        z = std::atan(m21 / m11);
    }
    double div = std::sqrt(1.0 - m31 * m31);
    if ((std::abs(div) < tol) && (std::abs(m31) < tol)) {
        y = 0.0;
    }
    else if ((std::abs(div) < tol) && (std::abs(m31) > tol)) {
        y = sign(0.5, -m31) * Constants::pi;
    }
    else {
        y = std::atan(-m31 / div);
    }
    if ((std::abs(m33) < tol) && (std::abs(m32) < tol)) {
        x = 0.0;
    }
    else if ((std::abs(m33) < tol) && (std::abs(m32) > tol)) {
        x = sign(0.5, m32) * Constants::pi;
    }
    else {
        x = std::atan(m32 / m33);
    }
    z *= 180.0 / Constants::pi;
    y *= 180.0 / Constants::pi;
    x *= 180.0 / Constants::pi;

    return {z, y, x};
}

Numlib::Vec<double> Numlib::eul2quat(double z, double y, double x)
{
    using namespace Numlib;

    Vec<double> result = {1.0, 0.0, 0.0, 0.0};
    if ((z != 0.0) || (y != 0.0) || (x != 0.0)) {
        const double tol = 2.0 * std::numeric_limits<double>::epsilon();

        z *= 0.5 * Constants::pi / 180.0;
        y *= 0.5 * Constants::pi / 180.0;
        x *= 0.5 * Constants::pi / 180.0;

        double c1 = std::cos(z);
        double c2 = std::cos(y);
        double c3 = std::cos(x);
        double s1 = std::sin(z);
        double s2 = std::sin(y);
        double s3 = std::sin(x);

        double q0 = s1 * s2 * s3 + c1 * c2 * c3;
        if (std::abs(q0) < tol) {
            q0 = 0.0;
        }
        double q1 = -s1 * s2 * c3 + s3 * c1 * c2;
        if (std::abs(q1) < tol) {
            q1 = 0.0;
        }
        double q2 = s1 * s3 * c2 + s2 * c1 * c3;
        if (std::abs(q2) < tol) {
            q2 = 0.0;
        }
        double q3 = s1 * c2 * c3 - s2 * s3 * c1;
        if (std::abs(q3) < tol) {
            q3 = 0.0;
        }
        result = {q0, q1, q2, q3};
    }
    return result;
}
