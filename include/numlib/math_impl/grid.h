//
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_GRID_H
#define NUMLIB_GRID_H

#include <numlib/traits.h>
#include <numlib/matrix.h>
#include <iostream>
#include <cassert>

namespace Numlib {

// Create two-dmensional grid based on coordinates in the vectors x and y.
template <typename T>
void meshgrid(const Vec<T>& x, const Vec<T>& y, Mat<T>& xx, Mat<T>& yy)
{
    Index nr = y.size();
    Index nc = x.size();

    xx.resize(nr, nc);
    yy.resize(nr, nc);

    for (Index i = 0; i < nr; ++i) {
        xx.row(i) = x;
    }
    for (Index j = 0; j < nc; ++j) {
        yy.column(j) = y;
    }
}

// Create three-dmensional grid based on coordinates in the vectors x, y and z.
template <typename T>
void meshgrid(const Vec<T>& x,
              const Vec<T>& y,
              const Vec<T>& z,
              Cube<T>& xx,
              Cube<T>& yy,
              Cube<T>& zz)
{
    Index n1 = y.size();
    Index n2 = x.size();
    Index n3 = z.size();

    xx.resize(n1, n2, n3);
    yy.resize(n1, n2, n3);
    zz.resize(n1, n2, n3);

    for (Index i = 0; i < n1; ++i) {
        for (Index j = 0; j < n2; ++j) {
            for (Index k = 0; k < n3; ++k) {
                xx(i, j, k) = x(j);
            }
        }
    }
    for (Index i = 0; i < n1; ++i) {
        for (Index j = 0; j < n2; ++j) {
            for (Index k = 0; k < n3; ++k) {
                yy(i, j, k) = y(i);
            }
        }
    }
    for (Index i = 0; i < n1; ++i) {
        for (Index j = 0; j < n2; ++j) {
            for (Index k = 0; k < n3; ++k) {
                zz(i, j, k) = z(k);
            }
        }
    }
}

// Class for handling grids with arithmetic progression.
//
class Grid {
public:
    using size_type = Index;

    Grid() : a0{0}, an{0}, d{0}, n{0} {}

    // Generate [0, 1, ..., amax].
    Grid(double amax) : a0{0}, an{amax}, d{1}
    {
        n = 1 + static_cast<size_type>(amax);
    }

    // Generate [amin, amin + 1, .... amax].
    Grid(double amin, double amax) : a0{amin}, an{amax}, d{1}
    {
        n = 1 + narrow_cast<size_type>(amax - amin);
    }

    // Generate the list [amin, amin + dd, ..., amax].
    Grid(double amin, double amax, double dd) : a0{amin}, an{amax}, d{dd}
    {
        n = 1 + narrow_cast<size_type>((amax - amin) / dd);
    }

    // Initialize data by reading an input file.
    Grid(std::istream& from, const std::string& key) : a0{0}, an{0}, d{0}, n{0}
    {
        set(from, key);
    }

    // Set new grid data by reading an input file.
    void set(std::istream& from, const std::string& key);

    // Set new grid data.
    void set(double amin, double amax, double dd);

    size_type size() const { return n; }
    double start() const { return a0; }
    double max() const { return an; }
    double step() const { return d; }
    bool empty() const { return n == 0; }

    double operator[](size_type i) const;
    double operator()(size_type i) const;

private:
    double a0;   // start value
    double an;   // maximum value
    double d;    // step size (common difference)
    size_type n; // number of steps (size of grid)
};

inline void Grid::set(double amin, double amax, double dd)
{
    a0 = amin;
    an = amax;
    d = dd;
    n = 1 + static_cast<size_type>((amax - amin) / dd);
}

inline double Grid::operator[](size_type i) const
{
    assert(i < n);
    return a0 + i * d;
}

inline double Grid::operator()(size_type i) const
{
    assert(i < n);
    return a0 + i * d;
}

// Non-member operators:

std::ostream& operator<<(std::ostream& to, const Grid& g);

} // namespace Numlib

#endif // NUMLIB_GRID_H
