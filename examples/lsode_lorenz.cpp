#include <numlib/math.h>
#include <iostream>

void lorenz(int* /* neq */, double* /* t */, double* y, double* ydot)
{
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0 / 3.0;

    ydot[0] = sigma * (y[1] - y[0]);
    ydot[1] = R * y[0] - y[1] - y[0] * y[2];
    ydot[2] = -b * y[2] + y[0] * y[1];
}

int main()
{
#ifdef ENABLE_ODEPACK
    Numlib::Lsode ode(lorenz);

    Numlib::Vec<double> y = {10.0, 1.0, 1.0};

    double t0 = 0.0;
    double t1 = 0.1;

    for (int i = 0; i < 5; ++i) {
        ode.integrate(t0, t1, y);
        std::cout << "At t = " << t0 << ", y = " << y(0) << " " << y(1) << " "
                  << y(2) << '\n';
        t1 += 0.1;
    }
#else
    std::cout << "Fortran compiler is needed\n";
#endif
}
