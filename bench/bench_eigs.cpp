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

using Timer = std::chrono::duration<double, std::milli>;

void print(int n, const Timer& t_arma, const Timer& t_eigs)
{
    std::cout << "Eigenvalues for symmetric matrix:\n"
              << "---------------------------------\n"
              << "size =      " << n << " x " << n << '\n'
              << "eigs/arma = " << t_eigs.count() / t_arma.count() << "\n\n";
}

void benchmark(int n)
{
    // arma::mat a1 = arma::randu<arma::mat>(n, n);
    // arma::mat a2 = a1.t() * a1;
    // arma::mat eigvec(n, n);
    // arma::vec eigval(n);
    auto t1 = std::chrono::high_resolution_clock::now();
    // arma::eig_sym(eigval, eigvec, a2);
    auto t2 = std::chrono::high_resolution_clock::now();
    Timer t_arma = t2 - t1;

    Numlib::Mat<double> b1 = Numlib::randu<Numlib::Mat<double>>(n, n);
    Numlib::Mat<double> b2 = Numlib::transpose(b1) * b1;
    Numlib::Vec<double> wr(n);
    t1 = std::chrono::high_resolution_clock::now();
    Numlib::eigs(b2, wr);
    t2 = std::chrono::high_resolution_clock::now();
    Timer t_eigs = t2 - t1;

    print(n, t_arma, t_eigs);
}

int main()
{
    int n = 10;
    benchmark(n);

    n = 100;
    benchmark(n);

    n = 500;
    benchmark(n);
}
