// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <cassert>

double Numlib::Math::trapz(double xlo, double xup, const Vec<double>& y)
{
    assert(!y.empty());

    const double step = std::abs(xup - xlo) / static_cast<double>(y.size() - 1);
    double ans = 0.0;

    for (std::size_t i = 1; i < y.size(); ++i) {
        ans += 0.5 * (y(i) + y(i - 1));
    }
    return ans *= step;
}
