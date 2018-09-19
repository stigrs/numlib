#include <iostream>
#include <numlib/math.h>

double f(double x) { return x * x; }

int main() { std::cout << Numlib::dfdx(f, 2.0) << '\n'; }
