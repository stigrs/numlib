#include <iostream>
#include <numlib/matrix.h>
#include <numlib/math.h>

int main()
{
    using namespace Numlib;

    Mat<double> a = {{0.12, -8.19, 7.69, -2.26, -4.71},
                     {-6.91, 2.22, -5.12, -9.08, 9.96},
                     {-3.33, -8.94, -6.72, -4.40, -9.98},
                     {3.97, 3.33, -2.74, -7.92, -3.20}};

    Mat<double> b = {{7.30, 0.47, -6.28},
                     {1.33, 6.58, -3.42},
                     {2.68, -1.71, 3.46},
                     {-9.62, -0.79, 0.41},
                     {0.00, 0.00, 0.00}};

    lstsq(a, b);
    std::cout << b << '\n';
}
