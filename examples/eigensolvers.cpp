#include <iostream>
#include <complex>
#include <numlib/matrix.h>
#include <numlib/math.h>

int main()
{
    using namespace Numlib;

    Mat<double> sa = hilbert(5);
    Vec<double> eval;

    eigs(sa, eval);
    std::cout << "Eigenvalues for symmetric matrix:\n"
              << eval << "\n\n"
              << "Eigenvectors for symmetric matrix:\n"
              << sa << '\n';

    Mat<double> ga = {{1.0, 5.0, 4.0, 2.0},
                     {-2.0, 3.0, 6.0, 4.0},
                     {5.0, 1.0, 0.0, -1.0},
                     {2.0, 3.0, -4.0, 0.0}};

    Vec<std::complex<double>> geval;
    Mat<std::complex<double>> gevec;

    eig(ga, gevec, geval);
    std::cout << "Eigenvalues for general matrix:\n"
              << geval << "\n\n"
              << "Right eigenvectors for general matrix:\n"
              << gevec << '\n';
}
