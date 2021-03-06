#include <iostream>
#include <numlib/matrix.h>
#include <numlib/math.h>

int main()
{
    using namespace Numlib;

    Vec<double> a = {3.0,  13.0, 7.0,  5.0,  21.0, 23.0, 39.0, 23.0,
                     40.0, 23.0, 14.0, 12.0, 56.0, 23.0, 29.0};

    std::cout << "mean(a) =   " << mean(a) << '\n'
              << "median(a) = " << median(a) << '\n'
              << "stddev(a) = " << stddev(a) << '\n'
              << "rms(a) =    " << rms(a) << "\n\n";

    Vec<double> b = {3.0,  13.0, 7.0,  5.0,  21.0, 23.0, 23.0,
                     40.0, 23.0, 14.0, 12.0, 56.0, 23.0, 29.0};
    Vec<double> c = {3.0,  13.0, 7.0,  5.0,  21.0, 23.0, 39.0,
                     23.0, 40.0, 23.0, 14.0, 12.0, 56.0, 23.0};

    std::cout << "cov(b, c) = " << cov(b, c) << "\n\n";

    Mat<double> p = {{0.8147, 0.9058, 0.1270}, {0.9134, 0.6324, 0.0975},
                     {0.2785, 0.5469, 0.9575}, {0.9649, 0.1576, 0.9706},
                     {0.9572, 0.4854, 0.8003}, {0.1419, 0.4218, 0.9157}};

    Mat<double> q = {{0.7922, 0.9595, 0.6557}, {0.0357, 0.8491, 0.9340},
                     {0.6787, 0.7577, 0.7431}, {0.3922, 0.6555, 0.1712},
                     {0.7060, 0.0318, 0.2769}, {0.0462, 0.0971, 0.8235}};

    std::cout << "rmsd(p, q) = " << kabsch_rmsd(p, q) << '\n';
}
