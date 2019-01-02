#include <iostream>
#include <cmath>
#include <numlib/matrix.h>
#include <numlib/constants.h>
#include <numlib/math.h>

int main()
{
    using namespace Numlib;

    double xlo = 2.1;
    double xup = 3.6;

    Vec<double> y = {3.2, 2.7, 2.9, 3.5, 4.1, 5.2};

    std::cout << "trapezoidal: " << trapz(xlo, xup, y) << '\n';

    double a = 0.0;
    double b = Constants::pi;

    std::cout << "5-point gaussian quadrature: "
              << quad<5>([](double x) { return std::sin(x); }, a, b) << '\n';

#ifdef ENABLE_QUADPACK
    std::cout << "qags: " << qags([](double& x) { return std::sin(x); }, a, b)
              << '\n';
#endif
}
