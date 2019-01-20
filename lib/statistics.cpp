// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <algorithm>
#include <cmath>

double Numlib::harmmean(const Numlib::Vec<double>& x)
{
    double sumi = 0.0;
    for (Index i = 0; i < x.size(); ++i) {
        assert(x(i) != 0.0);
        sumi += 1.0 / x(i);
    }
    assert(sumi != 0.0);
    return x.size() / sumi;
}

double Numlib::median(Numlib::Vec<double>& x)
{
    auto first = x.begin();
    auto last = x.end();
    auto mid = first + (last - first) / 2;

    std::nth_element(first, mid, last);
    double med = *mid;

    if ((x.size() % 2) == 0) { // size is even
        mid = first + (last - first) / 2 - 1;
        std::nth_element(first, mid, last);
        med = (med + *mid) / 2.0;
    }
    return med;
}

double Numlib::var(const Numlib::Vec<double>& x)
{
    // Two-pass algorithm:
    double n = static_cast<double>(x.size());
    double xmean = mean(x);
    double sum2 = 0.0;

    for (Index i = 0; i < x.size(); ++i) {
        sum2 += std::pow(x(i) - xmean, 2.0);
    }
    return sum2 / (n - 1.0);
}

double Numlib::mad(const Numlib::Vec<double>& x)
{
    double xmean = mean(x);
    double sumdev = 0.0;

    for (Index i = 0; i < x.size(); ++i) {
        sumdev = std::abs(x(i) - xmean);
    }
    return sumdev / x.size();
}

double Numlib::rms(const Numlib::Vec<double>& x)
{
    double sum2 = 0.0;

    for (Index i = 0; i < x.size(); ++i) {
        sum2 += x(i) * x(i);
    }
    return std::sqrt(sum2 / x.size());
}

double Numlib::rmsd(const Numlib::Mat<double>& a, const Numlib::Mat<double>& b)
{
    assert(a.rows() == b.rows());
    assert(a.cols() == b.cols());

    double sum2 = 0.0;
    for (Index i = 0; i < a.rows(); ++i) {
        double dist = 0.0;
        for (Index j = 0; j < a.cols(); ++j) {
            dist += std::pow(a(i, j) - b(i, j), 2.0);
        }
        sum2 += dist;
    }
    return std::sqrt(sum2 / a.rows());
}

double Numlib::cov(const Numlib::Vec<double>& x, const Numlib::Vec<double>& y)
{
    assert(x.size() == y.size() && !x.empty());

    double xmean = mean(x);
    double ymean = mean(y);
    double cov = 0.0;

    for (Index i = 0; i < x.size(); ++i) {
        double a = x(i) - xmean;
        double b = y(i) - ymean;
        cov += a * b / (x.size() - 1.0);
    }
    return cov;
}

double Numlib::kabsch_rmsd(const Numlib::Mat<double>& p,
                           const Numlib::Mat<double>& q)
{
    assert(same_extents(p, q));

    // Centroid of P:
    Vec<double> p0 = centroid(p);

    // Centroid of Q:
    Vec<double> q0 = centroid(q);

    // Translate P to center the origin:
    Mat<double> pc(p);
    translate(pc, -p0(0), -p0(1), -p0(2));

    // Translate Q to center the origin:
    Mat<double> qc(q);
    translate(qc, -q0(0), -q0(1), -q0(2));

    // Cross-covariance matrix:
    Mat<double> h = transpose(pc) * qc;

    // Singular value decomposition:
    Vec<double> s;
    Mat<double> u;
    Mat<double> vt;

    svd(h, s, u, vt);

    Mat<double> v = transpose(vt);
    Mat<double> ut = transpose(u);

    // Ensure right-handedness:
    double d = det(v * ut);
    Mat<double> eye = identity(3);
    eye(2, 2) = std::copysign(1.0, d);

    // Optimal rotation matrix:
    Mat<double> rotm = v * eye * ut;

    // Rotate coordinates:
    rotate(pc, rotm);

    // Compute RMSD:
    return rmsd(pc, qc);
}
