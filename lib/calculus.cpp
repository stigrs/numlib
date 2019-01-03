// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

double Numlib::trapz(double xlo, double xup, const Vec<double>& y)
{
    assert(!y.empty());

    const double step = std::abs(xup - xlo) / (y.size() - 1);
    double ans = 0.0;

#pragma omp parallel for shared(y) reduction(+ : ans)
    for (Index i = 1; i < y.size(); ++i) {
        ans += 0.5 * (y(i) + y(i - 1));
    }
    return ans *= step;
}

void Numlib::rk4(
    std::function<Numlib::Vec<double>(double t, const Numlib::Vec<double>&)> f,
    Numlib::Vec<double>& y,
    double& t0,
    double t1,
    double dt)
{
    int nsteps = Numlib::round<int>((t1 - t0) / dt);

    for (int i = 0; i < nsteps; ++i) {
        auto k1 = dt * f(t0, y);
        auto k2 = dt * f(t0 + dt / 2.0, y + k1 / 2.0);
        auto k3 = dt * f(t0 + dt / 2.0, y + k2 / 2.0);
        auto k4 = dt * f(t0 + dt, y + k3);
        y += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        t0 += dt;
    }
}

void Numlib::rk4(
    std::function<Numlib::Vec<double>(double t, const Numlib::Vec<double>&)> f,
    Numlib::Vec<double>& y,
    double& t0,
    double t1,
    int nsteps)
{
    double dt = (t1 - t0) / nsteps;

    for (int i = 0; i < nsteps; ++i) {
        auto k1 = dt * f(t0, y);
        auto k2 = dt * f(t0 + dt / 2.0, y + k1 / 2.0);
        auto k3 = dt * f(t0 + dt / 2.0, y + k2 / 2.0);
        auto k4 = dt * f(t0 + dt, y + k3);
        y += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        t0 += dt;
    }
}

