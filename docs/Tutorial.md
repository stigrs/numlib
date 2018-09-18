## Tutorial

### Basic Matrix Uses

Example program:

	#include <iostream>
	#include <numlib/matrix.h>

	int main()
	{
	    using namespace Numlib;

	    Mat<int> m(3, 3);
	    std::cout << "rank =      " << m.rank() << '\n'
	              << "size =      " << m.size() << '\n'
	              << "rows =      " << m.rows() << '\n'
	              << "cols =      " << m.cols() << '\n'
	              << "extent(0) = " << m.extent(0) << '\n'
	              << "extent(1) = " << m.extent(1) << '\n';
	}

	rank =      2
	size =      9
	rows =      3
	cols =      3
	extent(0) = 3
	extent(1) = 3
