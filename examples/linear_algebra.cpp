#include <iostream>
#include <numlib/matrix.h>
#include <numlib/math.h>

int main()
{
    using namespace Numlib;

    Vec<double> x = linspace(1.0, 10.0, 10);
    std::cout << "x = \n" << x << "\n\n";

    Vec<double> y = 2.0 * ones<Vec<double>>(10);

    std::cout << "min(x) =    " << min(x) << '\n'
              << "max(x) =    " << max(x) << '\n'
              << "sum(x) =    " << sum(x) << '\n'
              << "prod(x) =   " << prod(x) << '\n'
              << "norm(x) =   " << norm(x) << '\n'
              << "dot(x, y) = " << dot(x, y) << "\n\n";

    Vec<double> a = {3.0, -3.0, 1.0};
    Vec<double> b = {4.0, 9.0, 2.0};
    std::cout << "cross(a, b) = \n" << cross(a, b) << '\n';
}
