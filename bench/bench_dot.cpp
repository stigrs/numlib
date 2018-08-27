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
#include <valarray>
#include <numeric>

using Timer = std::chrono::duration<double, std::micro>;

void print(int n, const Timer& t_arma, const Timer& t_num, const Timer& t_val)
{
    std::cout << "Dot product:\n"
              << "------------\n"
              << "size =            " << n << '\n'
              << "numlib/arma =     " << t_num.count() / t_arma.count() << "\n"
              << "numlib/valarray = " << t_num.count() / t_val.count()
              << "\n\n";
}

void benchmark(int n)
{
    arma::vec aa(n);
    arma::vec ab(n);
    aa.fill(1.0);
    ab.fill(2.0);
    auto t1 = std::chrono::high_resolution_clock::now();
    arma::dot(aa, ab);
    auto t2 = std::chrono::high_resolution_clock::now();
    Timer t_arma = t2 - t1;

    Numlib::Vec<double> na(n);
    Numlib::Vec<double> nb(n);
    na = 1.0;
    nb = 2.0;
    t1 = std::chrono::high_resolution_clock::now();
    Numlib::dot(na, nb);
    t2 = std::chrono::high_resolution_clock::now();
    Timer t_num = t2 - t1;

    std::valarray<double> va(1.0, n);
    std::valarray<double> vb(2.0, n);

    t1 = std::chrono::high_resolution_clock::now();
    std::inner_product(std::begin(va), std::end(va), std::begin(vb), 0.0);
    t2 = std::chrono::high_resolution_clock::now();
    Timer t_val = t2 - t1;

    print(n, t_arma, t_num, t_val);
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
