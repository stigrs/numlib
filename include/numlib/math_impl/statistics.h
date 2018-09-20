// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_STATISTICS_H
#define NUMLIB_MATH_STATISTICS_H

#include <numlib/matrix.h>
#include <cmath>

namespace Numlib {

// Arithmetic mean.
inline double mean(const Vec<double>& x) { return sum(x) / x.size(); }

// Geometric mean.
inline double geomean(const Vec<double>& x)
{
    return std::pow(prod(x), static_cast<double>(1.0 / x.size()));
}

// Harmonic mean.
double harmmean(const Vec<double>& x);

// Median.
double median(Vec<double>& x);

// Variance.
double var(const Vec<double>& x);

// Standard deviation.
inline double stddev(const Vec<double>& x) { return std::sqrt(var(x)); }

// Covariance.
double cov(const Vec<double>& x, const Vec<double>& y);

// Mean absolute deviation.
double mad(const Vec<double>& x);

// Root-mean-square deviation.
double rms(const Vec<double>& x);

// Root-mean-square displacement.
double rmsd(const Mat<double>& a, const Mat<double>& b);

// Minimum root-mean-square displacement between two paired set of points.
//
// Algorithm:
//   Kabsch algorithm (https://en.wikipedia.org/wiki/Kabsch_algorithm)
//
double kabsch_rmsd(const Mat<double>& p, const Mat<double>& q);

} // namespace Numlib

#endif // NUMLIB_MATH_STATISTICS_H
