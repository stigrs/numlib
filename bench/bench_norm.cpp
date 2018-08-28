// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <numlib/math.h>
#include <armadillo>
#include <chrono>
#include <iostream>

typedef std::chrono::duration<double, std::milli> Timer;

void print(int n, const Timer& t_arma, const Timer& t_numlib)
{
    std::cout << "Vector norm:\n"
              << "---------------\n"
              << "size =        " << n << '\n'
              << "numlib/arma = " << t_numlib.count() / t_arma.count()
              << "\n\n";
}

void benchmark(int n)
{
    arma::vec aa(n);
    aa.fill(1.0);
    double nrm;
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int it = 0; it < 10; ++it) {
        nrm = arma::norm(aa);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    Timer t_arma = t2 - t1;

    Numlib::Vec<double> va(n);
    va = 1.0;

    t1 = std::chrono::high_resolution_clock::now();
    for (int it = 0; it < 10; ++it) {
        nrm = norm(va);
    }
    t2 = std::chrono::high_resolution_clock::now();
    Timer t_numlib = t2 - t1;

    print(n, t_arma, t_numlib);
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
