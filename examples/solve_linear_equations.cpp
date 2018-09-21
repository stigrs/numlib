#include <iostream>
#include <numlib/matrix.h>
#include <numlib/math.h>

int main()
{
    using namespace Numlib;

    Mat<double> A = {{1.0, 2.0, 3.0}, {2.0, 3.0, 4.0}, {3.0, 4.0, 1.0}};
    Mat<double> B = {{14.0}, {20.0}, {14.0}};

    linsolve(A, B);
    std::cout << B << '\n';
}
