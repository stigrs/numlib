// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <algorithm>
#include <limits>

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
    std::function<
        void(double t, const Numlib::Vec<double>&, Numlib::Vec<double>&)> f,
    Numlib::Vec<double>& y,
    double& t0,
    double t1,
    double dt)
{
    int nsteps = Numlib::round<int>((t1 - t0) / dt);

    Numlib::Vec<double> k1(y.size());
    Numlib::Vec<double> k2(y.size());
    Numlib::Vec<double> k3(y.size());
    Numlib::Vec<double> k4(y.size());
    Numlib::Vec<double> yn(y.size());

    for (int it = 0; it < nsteps; ++it) {
        f(t0, y, k1);
#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            yn(i) = y(i) + 0.5 * dt * k1(i);
        }

        f(t0 + dt / 2.0, yn, k2);
#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            yn(i) = y(i) + 0.5 * dt * k2(i);
        }

        f(t0 + dt / 2.0, yn, k3);
#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            yn(i) = y(i) + dt * k3(i);
        }

        f(t0 + dt, yn, k4);
#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            y(i) += dt * (k1(i) + 2.0 * k2(i) + 2.0 * k3(i) + k4(i)) / 6.0;
        }
        t0 += dt;
    }
}

void Numlib::rk4(
    std::function<
        void(double t, const Numlib::Vec<double>&, Numlib::Vec<double>&)> f,
    Numlib::Vec<double>& y,
    double& t0,
    double t1,
    int nsteps)
{
    double dt = (t1 - t0) / nsteps;

    Numlib::Vec<double> k1(y.size());
    Numlib::Vec<double> k2(y.size());
    Numlib::Vec<double> k3(y.size());
    Numlib::Vec<double> k4(y.size());
    Numlib::Vec<double> yn(y.size());

    for (int it = 0; it < nsteps; ++it) {
        f(t0, y, k1);
#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            yn(i) = y(i) + 0.5 * dt * k1(i);
        }

        f(t0 + dt / 2.0, yn, k2);
#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            yn(i) = y(i) + 0.5 * dt * k2(i);
        }

        f(t0 + dt / 2.0, yn, k3);
#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            yn(i) = y(i) + dt * k3(i);
        }

        f(t0 + dt, yn, k4);
#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            y(i) += dt * (k1(i) + 2.0 * k2(i) + 2.0 * k3(i) + k4(i)) / 6.0;
        }
        t0 += dt;
    }
}

void Numlib::dopri5(
    std::function<
        void(double t, const Numlib::Vec<double>&, Numlib::Vec<double>&)> f,
    Numlib::Vec<double>& y,
    double& t0,
    double t1,
    double atol,
    double rtol,
    int maxstep)
{
    const double a21 = 1.0 / 5.0;
    const double a31 = 3.0 / 40.0;
    const double a32 = 9.0 / 40.0;
    const double a41 = 44.0 / 45.0;
    const double a42 = -56.0 / 15.0;
    const double a43 = 32.0 / 9.0;
    const double a51 = 19372.0 / 6561.0;
    const double a52 = -25360.0 / 2187.0;
    const double a53 = 64448.0 / 6561.0;
    const double a54 = -212.0 / 729.0;
    const double a61 = 9017.0 / 3168.0;
    const double a62 = -355.0 / 33.0;
    const double a63 = 46732.0 / 5247.0;
    const double a64 = 49.0 / 176.0;
    const double a65 = -5103.0 / 18656.0;
    const double a71 = 35.0 / 384.0;
    const double a72 = 0.0;
    const double a73 = 500.0 / 1113.0;
    const double a74 = 125.0 / 192.0;
    const double a75 = -2187.0 / 6784.0;
    const double a76 = 11.0 / 84.0;

    const double b1p = 35.0 / 384.0;
    const double b2p = 0.0;
    const double b3p = 500.0 / 1113.0;
    const double b4p = 125.0 / 192.0;
    const double b5p = -2187.0 / 6784.0;
    const double b6p = 11.0 / 84.0;
    const double b7p = 0.0;

    const double b1 = 5179.0 / 57600.0;
    const double b2 = 0.0;
    const double b3 = 7571.0 / 16695.0;
    const double b4 = 393.0 / 640.0;
    const double b5 = -92097.0 / 339200.0;
    const double b6 = 187.0 / 2100.0;
    const double b7 = 1.0 / 40.0;

    const double c1 = 0.0;
    const double c2 = 1.0 / 5.0;
    const double c3 = 3.0 / 10.0;
    const double c4 = 4.0 / 5.0;
    const double c5 = 8.0 / 9.0;
    const double c6 = 1.0;
    const double c7 = 1.0;

    const double eps = std::numeric_limits<double>::epsilon();
    const double hmin = 2.0 * eps;
    const double hmax = t1 - t0;

    double h = std::max(hmin, hmax);

    Numlib::Vec<double> k1(y.size());
    Numlib::Vec<double> k2(y.size());
    Numlib::Vec<double> k3(y.size());
    Numlib::Vec<double> k4(y.size());
    Numlib::Vec<double> k5(y.size());
    Numlib::Vec<double> k6(y.size());
    Numlib::Vec<double> k7(y.size());
    Numlib::Vec<double> yn(y.size());

    int istep = 0;
    while (istep < maxstep) {
        // Compute function values:
        f(t0 + c1 * h, y, k1);
#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            yn(i) = y(i) + h * (a21 * k1(i));
        }
        f(t0 + c2 * h, yn, k2);
#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            yn(i) = y(i) + h * (a31 * k2(i) + a32 * k2(i));
        }
        f(t0 + c3 * h, yn, k3);
#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            yn(i) = y(i) + h * (a41 * k3(i) + a42 * k3(i) + a43 * k3(i));
        }
        f(t0 + c4 * h, yn, k4);
#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            yn(i) = y(i) +
                    h * (a51 * k4(i) + a52 * k4(i) + a53 * k4(i) + a54 * k4(i));
        }
        f(t0 + c5 * h, yn, k5);
#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            yn(i) = y(i) + h * (a61 * k5(i) + a62 * k5(i) + a63 * k5(i) +
                                a64 * k5(i) + a65 * k5(i));
        }
        f(t0 + c6 * h, yn, k6);
#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            yn(i) = y(i) + h * (a71 * k6(i) + a72 * k6(i) + a73 * k6(i) +
                                a74 * k6(i) + a75 * k6(i) + a76 * k6(i));
        }
        f(t0 + c7 * h, yn, k7);

        // Compute solution:

#pragma omp parallel for
        for (Index i = 0; i < y.size(); ++i) {
            yn(i) =
                y(i) + h * (b1 * k1(i) + b2 * k2(i) + b3 * k3(i) + b4 * k4(i) +
                            b5 * k5(i) + b6 * k6(i) + b7 * k7(i));
        }
        // Compute error:

        double maxerr = 0.0;
        double maxtol = 0.0;
        for (Index i = 0; i < y.size(); ++i) {
            double ei = std::abs((b1 - b1p) * k1(i) + (b2 - b2p) * k2(i) +
                                 (b3 - b3p) * k3(i) + (b4 - b4p) * k4(i) +
                                 (b5 - b5p) * k5(i) + (b6 - b6p) * k6(i) +
                                 (b7 - b7p) * k7(i));
            double di = std::max(rtol * std::abs(yn(i)), atol);
            if (ei > maxerr) {
                maxerr = ei;
            }
            if (di > maxtol) {
                maxtol = di;
            }
        }
        if (maxerr <= maxtol) { // accept solution
            y = yn;
            t0 += h;
            istep = 0;
        }

        // Adjust stepsize (https://en.wikipedia.org/wiki/Adaptive_stepsize):

        double s = 0.9 * std::pow(maxtol / maxerr, 0.2);

        h *= std::min(std::max(s, 0.3), 2.0);

        if (t0 + h > t1) {
            h = t1 - t0;
        }
        if (h > hmax) {
            h = hmax;
        }
        else if (h < hmin) {
            h = hmin;
        }
        if (t0 >= t1) { // integration succeeded
            break;
        }
        ++istep;
    }
    if (t0 < t1) {
        throw Math_error("integration failed to converge");
    }
}

