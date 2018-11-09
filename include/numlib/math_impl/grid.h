//
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_GRID_H
#define NUMLIB_GRID_H

#include <numlib/traits.h>
#include <iostream>
#include <cassert>

namespace Numlib {

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

