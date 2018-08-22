// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <armadillo>
#include <chrono>
#include <iostream>
#include <valarray>

typedef std::chrono::duration<double, std::milli> Timer;

void print(int n, const Timer& t_arma, const Timer& t_numlib,
           const Timer& t_val)
{
    std::cout << "Vector addition:\n"
              << "----------------\n"
              << "size =        " << n << '\n'
              << "numlib/arma = " << t_numlib.count() / t_arma.count() << "\n"
              << "numlib/val =  " << t_numlib.count() / t_val.count() << "\n\n";
}

void benchmark(int n)
{
    arma::vec aa(n);
    arma::vec ab(n);
    aa.fill(1.0);
    ab.fill(1.0);
    auto t1 = std::chrono::high_resolution_clock::now();
    ab = 2.0 * aa + ab;
    auto t2 = std::chrono::high_resolution_clock::now();
    Timer t_arma = t2 - t1;

    num::Matrix<double, 1> va(n);
    num::Matrix<double, 1> vb(n);

    va = 1.0;
    vb = 1.0;

    t1 = std::chrono::high_resolution_clock::now();
    vb = 2.0 * va + vb;
    t2 = std::chrono::high_resolution_clock::now();
    Timer t_numlib = t2 - t1;

    std::valarray<double> wa(1.0, n);
    std::valarray<double> wb(1.0, n);
    t1 = std::chrono::high_resolution_clock::now();
    wb = 2.0 * wa + wb;
    t2 = std::chrono::high_resolution_clock::now();
    Timer t_val = t2 - t1;

    print(n, t_arma, t_numlib, t_val);

    for (int i = 0; i < n; ++i) {
        if (ab(i) != vb(i)) {
            std::cout << "Different\n";
        }
    }
}

int main()
{
    int n = 10;
    benchmark(n);

    n = 100;
    benchmark(n);

    n = 1000;
    benchmark(n);

    n = 10000;
    benchmark(n);

    n = 100000;
    benchmark(n);
}
