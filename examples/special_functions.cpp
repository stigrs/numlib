#include <iostream>
#include <numlib/math.h>
#include <complex>

int main()
{
    using namespace Faddeeva;

    std::complex<double> z(-1.0, 1.0);
    std::cout << "Faddeeva function: " << w(z) << '\n';

    z = {1.0, 2.0};
    std::cout << "Error function: " << erf(z) << '\n';

    z = {2.0, 1.0};
    std::cout << "Dawson function: " << Dawson(z) << '\n';
}
