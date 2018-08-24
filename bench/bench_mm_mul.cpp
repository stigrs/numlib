// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <armadillo>
#include <chrono>
#include <iostream>

using Timer = std::chrono::duration<double, std::milli>;

void print(int n, int m, const Timer& t_arma, const Timer& t_mm_mul)
{
    std::cout << "Matrix-matrix multiplication:\n"
              << "-----------------------------\n"
              << "size =        " << n << " x " << m << '\n'
              << "mm_mul/arma = " << t_mm_mul.count() / t_arma.count()
              << "\n\n";
}

void benchmark(int n, int m)
{
    using namespace Numlib;

    arma::mat a1 = arma::ones<arma::mat>(n, m);
    arma::mat a2 = arma::ones<arma::mat>(m, n);
    auto t1 = std::chrono::high_resolution_clock::now();
    arma::mat a3 = a1 * a2;
    auto t2 = std::chrono::high_resolution_clock::now();
    Timer t_arma = t2 - t1;

    Mat<double> b1 = ones<Mat<double>>(n, m);
    Mat<double> b2 = ones<Mat<double>>(m, n);
    t1 = std::chrono::high_resolution_clock::now();
    auto b3 = b1 * b2;
    t2 = std::chrono::high_resolution_clock::now();
    Timer t_mm_mul = t2 - t1;

    print(n, m, t_arma, t_mm_mul);
}

int main()
{
    int n = 10;
    int m = 5;
    benchmark(n, m);

    n = 100;
    m = 50;
    benchmark(n, m);

    n = 1000;
    m = 500;
    benchmark(n, m);
}
