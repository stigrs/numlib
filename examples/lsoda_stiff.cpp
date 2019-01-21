#include <numlib/math.h>
#include <iostream>
#include <iomanip>

void fsys(double /* t */, double* y, double* ydot, void* /* data */)
{
    ydot[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
    ydot[2] = 3.0e7 * y[1] * y[1];
    ydot[1] = -ydot[0] - ydot[2];
}

int main()
{
    Numlib::Vec<double> rtol = {1.0e-4, 1.0e-8, 1.0e-4};
    Numlib::Vec<double> atol = {1.0e-6, 1.0e-10, 1.0e-6};

    Numlib::Odeint ode(fsys, rtol, atol);

    Numlib::Vec<double> y = {1.0, 0.0, 0.0};

    double t0 = 0.0;
    double t1 = 0.4;

    for (int i = 0; i < 12; ++i) {
        ode.integrate(y, t0, t1);
        std::cout << std::scientific << std::setprecision(6) << "At t = " << t0
                  << ", y = " << y(0) << " " << y(1) << " " << y(2) << '\n';
        t1 *= 10.0;
    }
}
