// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <algorithm>

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
