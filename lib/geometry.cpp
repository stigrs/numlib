// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>

double Numlib::hypot(const double a, const double b)
{
    double aa = std::abs(a);
    double ab = std::abs(b);

    if (aa > ab) {
        return aa * std::sqrt(1.0 + std::pow(ab / aa, 2.0));
    }
    return ab == 0.0 ? 0.0 : ab * std::sqrt(1.0 + std::pow(aa / ab, 2.0));
}

double Numlib::dihedral(const Numlib::Vec<double>& a,
                        const Numlib::Vec<double>& b,
                        const Numlib::Vec<double>& c,
                        const Numlib::Vec<double>& d)
{
    auto ab = Numlib::normalize(b - a);
    auto bc = Numlib::normalize(c - b);
    auto cd = Numlib::normalize(d - c);
    auto n1 = Numlib::cross(ab, bc);
    auto n2 = Numlib::cross(bc, cd);
    auto m = Numlib::cross(n1, bc);
    double x = Numlib::dot(n1, n2);
    double y = Numlib::dot(m, n2);

    double tau = Numlib::radtodeg(std::atan2(y, x));
    if (std::abs(tau) < 1.0e-8) { // avoid very small angles close to zero
        tau = 0.0;
    }
    return tau;
}

void Numlib::pdist_matrix(Numlib::Mat<double>& dm,
                          const Numlib::Mat<double>& mat)
{
    dm.resize(mat.rows(), mat.rows());

    for (Index i = 0; i < dm.rows(); ++i) {
        for (Index j = i; j < dm.cols(); ++j) {
            dm(i, j) = 0.0;
            if (i != j) {
                dm(i, j) = Numlib::norm(mat.row(i) - mat.row(j));
                dm(j, i) = dm(i, j);
            }
        }
    }
}

void Numlib::translate(Numlib::Mat<double>& xyz,
                       double dx,
                       double dy,
                       double dz)
{
    assert(xyz.cols() == 3);

    for (Index i = 0; i < xyz.rows(); ++i) {
        xyz(i, 0) += dx;
        xyz(i, 1) += dy;
        xyz(i, 2) += dz;
    }
}

void Numlib::rotate(Numlib::Mat<double>& xyz, const Numlib::Mat<double>& rotm)
{
    assert(rotm.rows() == 3 && rotm.cols() == 3);

    for (Index i = 0; i < xyz.rows(); ++i) {
        xyz.row(i) = rotm * xyz.row(i);
    }
}
