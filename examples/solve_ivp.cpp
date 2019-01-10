#include <numlib/matrix.h>
#include <numlib/math.h>
#include <iostream>

void lorenz(const Numlib::Vec<double>& x,
            Numlib::Vec<double>& dxdt,
            const double /* t */)
{
    const double sigma(10.0);
    const double R(28.0);
    const double b(8.0 / 3.0);

    dxdt(0) = sigma * (x(1) - x(0));
    dxdt(1) = R * x(0) - x(1) - x(0) * x(2);
    dxdt(2) = -b * x(2) + x(0) * x(1);
}

int main()
{
    Numlib::Vec<double> x = {10.0, 1.0, 1.0};

    double t0 = 0.0;
    double t1 = 0.1;

    for (int i = 0; i < 5; ++i) {
        Numlib::solve_ivp(lorenz, x, t0, t1);
        t1 += 0.1;
        std::cout << "At t = " << t0 << ", y = " << x(0) << " " << x(1) << " "
                  << x(2) << '\n';
    }
}

