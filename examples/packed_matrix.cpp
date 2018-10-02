#include <iostream>
#include <numlib/matrix.h>
#include <numlib/math.h>

int main()
{
#ifdef USE_MKL
    using namespace Numlib;

    auto a = hilbert<>(5);
    Symm_mat<double, lower_triang> ap(a);

    Mat<double> evec;
    Vec<double> eval(5); // correct size must be allocated

    eigs(ap, evec, eval);
    std::cout << eval << '\n';
#else
    std::cout << "Intel MKL is required\n";
#endif
}
