#include <iostream>
#include <numlib/matrix.h>
#include <numlib/math.h>

int main()
{
#ifdef USE_MKL
    using namespace Numlib;

    Mat<double> A = {{0.0, 2.0, 0.0, 1.0, 0.0},
                     {4.0, -1.0, -1.0, 0.0, 0.0},
                     {0.0, 0.0, 0.0, 3.0, -6.0},
                     {-2.0, 0.0, 0.0, 0.0, 2.0},
                     {0.0, 0.0, 4.0, 2.0, 0.0}};

    Sp_mat<double> SA = gather(A);

    Mat<double> B = {{8.0}, {-1.0}, {-18.0}, {8.0}, {20.0}};
    Mat<double> x;

    linsolve(SA, B, x);
    std::cout << x << '\n';
#else
    std::cout << "Intel MKL is required\n";
#endif
}
