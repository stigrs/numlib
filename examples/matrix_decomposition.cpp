#include <iostream>
#include <numlib/matrix.h>
#include <numlib/math.h>

int main()
{
    using namespace Numlib;

    Mat<double> a = {
        {12.0, -51.0, 4.0}, {6.0, 167.0, -68.0}, {-4.0, 24.0, -41.0}};

    Mat<double> q;
    Mat<double> r;

    qr(a, q, r);

    Mat<double> qr = q * r;
    std::cout << "QR decomposition (A = QR):\n" << qr << '\n';

    Mat<double> m = {
        {8.79, 9.93, 9.83, 5.45, 3.16},   {6.11, 6.91, 5.04, -0.27, 7.98},
        {-9.15, -7.93, 4.86, 4.85, 3.01}, {9.57, 1.64, 8.83, 0.74, 5.8},
        {-3.49, 4.02, 9.80, 10.00, 4.27}, {9.84, 0.15, -8.99, -6.02, -5.31}};

    Vec<double> s;
    Mat<double> u;
    Mat<double> vt;

    svd(m, s, u, vt);

    Mat<double> sigma = zeros<Mat<double>>(m.rows(), m.cols());
    sigma.diag() = s;

    std::cout << "SVD decomposition (M = USVT):\n" << u * sigma * vt << '\n';
}
