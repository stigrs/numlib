// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math_impl/grid.h>
#include <string>

void Numlib::Grid::set(std::istream& from, const std::string& key)
{
    a0 = 0.0;
    d = 1.0;
    an = 100.0;

    from.clear();
    from.seekg(0, std::ios_base::beg); // search from the beginning

    std::string token;
    while (from >> token) {
        if (token == key) {
            while (from >> token) {
                if (token == "End") {
                    break;
                }
                if (token == "min") {
                    from >> a0;
                }
                if (token == "step") {
                    from >> d;
                }
                if (token == "max") {
                    from >> an;
                }
            }
        }
    }
    assert(an >= a0);
    assert(d > 0.0);

    n = 1 + narrow_cast<size_type>((an - a0) / d);
}

std::ostream& Numlib::operator<<(std::ostream& to, const Numlib::Grid& g)
{
    to << "Min value:\t" << g.start() << '\n'
       << "Max value:\t" << g.max() << '\n'
       << "Step size:\t" << g.step() << '\n';
    return to;
}

