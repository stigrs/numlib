
#include <numlib/matrix.h>
#include <numlib/math.h>
#include <armadillo>
#include <chrono>
#include <iostream>

using Timer = std::chrono::duration<double, std::milli>;

void print(int n, int m, const Timer& t_arma, const Timer& t_num)
{
    std::cout << "Matrix transpose:\n"
              << "-----------------\n"
              << "size =     " << n << " x " << m << '\n'
              << "numlib/arma = " << t_num.count() / t_arma.count() << "\n\n";
}

void benchmark(int n, int m)
{
    arma::mat m1 = arma::ones<arma::mat>(n, m);
    auto t1 = std::chrono::high_resolution_clock::now();
    m1.t();
    auto t2 = std::chrono::high_resolution_clock::now();
    Timer t_arma = t2 - t1;

    Numlib::Mat<double> m2(n, m);
    m2 = 1.0;
    t1 = std::chrono::high_resolution_clock::now();
    auto mt = transpose(m2);
    t2 = std::chrono::high_resolution_clock::now();
    Timer t_num = t2 - t1;

    print(n, m, t_arma, t_num);
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
