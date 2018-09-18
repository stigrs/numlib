# Tutorial

## Basic Matrix Uses

Example program:

	#include <iostream>
	#include <numlib/matrix.h>

	int main()
	{
		using namespace Numlib;

		Mat<int> m1(3, 3);
		std::cout << "rank =      " << m1.rank() << '\n'
				  << "size =      " << m1.size() << '\n'
				  << "rows =      " << m1.rows() << '\n'
				  << "cols =      " << m1.cols() << '\n'
				  << "extent(0) = " << m1.extent(0) << '\n'
				  << "extent(1) = " << m1.extent(1) << "\n\n";

		m1.resize(5, 5);
		std::cout << "rank =      " << m1.rank() << '\n'
				  << "size =      " << m1.size() << '\n'
				  << "rows =      " << m1.rows() << '\n'
				  << "cols =      " << m1.cols() << '\n'
				  << "extent(0) = " << m1.extent(0) << '\n'
				  << "extent(1) = " << m1.extent(1) << "\n\n";

		Mat<int> m2 = {{1, 2, 3}, {4, 5, 6}};
		std::cout << m2 << '\n';
	}

Generated output:

	rank =      2
	size =      9
	rows =      3
	cols =      3
	extent(0) = 3
	extent(1) = 3

	rank =      2
	size =      25
	rows =      5
	cols =      5
	extent(0) = 5
	extent(1) = 5

	2 x 3
	[        1         2         3
			 4         5         6 ]

