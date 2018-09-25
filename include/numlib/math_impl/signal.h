// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_SIGNAL_H
#define NUMLIB_MATH_SIGNAL_H

#include <numlib/matrix.h>

#ifdef USE_MKL
#include <complex>
#include <mkl.h>
#endif

namespace Numlib {

// Vector convolution.
template <typename M>
Enable_if<Matrix_type<M>(), Vec<typename M::value_type>> conv(const M& a,
                                                              const M& b)
{
    static_assert(M::order == 1, "bad rank for vector convolution");
    using value_type = typename M::value_type;

    const auto na = a.size();
    const auto nb = b.size();
    const auto nc = na + nb - 1;

    Vec<value_type> res(nc);

    for (Index i = 0; i < nc; ++i) {
        res(i) = value_type{0};
        auto jmin = (i >= (nb - 1)) ? (i - (nb - 1)) : 0;
        auto jmax = (i < (na - 1)) ? i : (na - 1);
        for (Index j = jmin; j <= jmax; ++j) {
            res(i) += a(j) * b(i - j);
        }
    }
    return res;
}

#ifdef USE_MKL
// Forward real-to-complex 1D discrete Fourier transform.
Vec<std::complex<double>> fft(const Vec<double>& x);
#endif

#ifdef USE_MKL
// Inverse complex-to-real 1D discrete Fourier transform.
Vec<double> ifft(const Vec<std::complex<double>>& y);
#endif

} // namespace Numlib

#endif // NUMLIB_MATH_SIGNAL_H
