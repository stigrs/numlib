#include <iostream>
#include <numlib/constants.h>
#include <numlib/matrix.h>
#include <numlib/math.h>

double moda(int k, int l, int m)
{
    return narrow_cast<double>(narrow_cast<long long>(k * l) % m);
}

int main()
{
#ifdef USE_MKL
    using namespace Numlib;

    const double twopi = 2.0 * Constants::pi;
    const int n = 6;
    int h = -1;

    Vec<double> x(n);

    double factor = (2 * (n - h) % n == 0) ? 1.0 : 2.0;
    for (int i = 0; i < n; ++i) {
        double phase = moda(i, h, n) / n;
        x(i) = factor * std::cos(twopi * phase) / n;
    }
    auto y = fft(x);
    std::cout << "x = " << x << "\n\n"
              << "y = fft(x):\n"
              << y << "\n\n"
              << "x = ifft(y):\n"
              << ifft(y) << '\n';

#else
    std::cout << "Intel MKL is required\n";
#endif
}
