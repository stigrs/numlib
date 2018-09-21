#include <iostream>
#include <numlib/matrix.h>
#include <numlib/math.h>

int main()
{
    using namespace Numlib;

    Mat<double> a = {{1.0, 5.0, 2.0, 0.0},
                     {5.0, 2.0, 5.0, 2.0},
                     {2.0, 5.0, 3.0, 5.0},
                     {0.0, 2.0, 5.0, 4.0}};

    Band_matrix<double> ab(/* kl = */ 3, /* ku = */ 3, a);
    Mat<double> evec;
    Vec<double> eval;

    eigs(ab, evec, eval);
    std::cout << eval << '\n';
}
