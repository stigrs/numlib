// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>

#ifdef USE_MKL
Numlib::Vec<std::complex<double>> Numlib::fft(const Numlib::Vec<double>& x)
{
    Numlib::Vec<std::complex<double>> y(x.size());

    MKL_LONG status = 0;
    MKL_LONG n = narrow_cast<MKL_LONG>(x.size());

    // Create DFTI descriptor:
    DFTI_DESCRIPTOR_HANDLE hand = 0;
    status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_REAL, 1, n);
    if (status != DFTI_NO_ERROR) {
        DftiFreeDescriptor(&hand);
        throw Math_error("DftiCreateDescriptor failed");
    }

    // Set configuration out-of-place:
    status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    if (status != DFTI_NO_ERROR) {
        DftiFreeDescriptor(&hand);
        throw Math_error("DftiSetValue out-of-place failed");
    }

    // Set configuration CCE storage:
    status =
        DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    if (status != DFTI_NO_ERROR) {
        DftiFreeDescriptor(&hand);
        throw Math_error("DftiSetValue CCE storage failed");
    }

    // Commit the descriptor:
    status = DftiCommitDescriptor(hand);
    if (status != DFTI_NO_ERROR) {
        DftiFreeDescriptor(&hand);
        throw Math_error("DftiCommitDescriptor failed");
    }

    // Compute forward transform:
    status = DftiComputeForward(hand, (void*) x.data(), y.data());
    if (status != DFTI_NO_ERROR) {
        DftiFreeDescriptor(&hand);
        throw Math_error("DftiComputeForward failed");
    }

    // Unpack CCE format:
    for (MKL_LONG i = n / 2 + 1; i < n; ++i) {
        y(i).real(y(n - i).real());
        y(i).imag((-1) * y(n - i).imag());
    }

    // Clean-up and return:
    DftiFreeDescriptor(&hand);
    return y;
}
#endif

#ifdef USE_MKL
Numlib::Vec<double> Numlib::ifft(const Numlib::Vec<std::complex<double>>& y)
{
    Numlib::Vec<double> x(y.size());

    MKL_LONG status = 0;
    MKL_LONG n = narrow_cast<MKL_LONG>(y.size());

    // Create DFTI descriptor:
    DFTI_DESCRIPTOR_HANDLE hand = 0;
    status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_REAL, 1, n);
    if (status != DFTI_NO_ERROR) {
        DftiFreeDescriptor(&hand);
        throw Math_error("DftiCreateDescriptor failed");
    }

    // Set configuration out-of-place:
    status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    if (status != DFTI_NO_ERROR) {
        DftiFreeDescriptor(&hand);
        throw Math_error("DftiSetValue out-of-place failed");
    }

    // Set configuration CCE storage:
    status =
        DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    if (status != DFTI_NO_ERROR) {
        DftiFreeDescriptor(&hand);
        throw Math_error("DftiSetValue CCE storage failed");
    }

    // Set backward scale factor so that backward transform yields inverse
    // of forward transform:
    double bscale = 1.0 / n;
    status = DftiSetValue(hand, DFTI_BACKWARD_SCALE, bscale);
    if (status != DFTI_NO_ERROR) {
        DftiFreeDescriptor(&hand);
        throw Math_error("DftiSetValue backward scale factor failed");
    }

    // Commit the descriptor:
    status = DftiCommitDescriptor(hand);
    if (status != DFTI_NO_ERROR) {
        DftiFreeDescriptor(&hand);
        throw Math_error("DftiCommitDescriptor failed");
    }

    // Compute backward transform:
    status = DftiComputeBackward(hand, (void*) y.data(), x.data());
    if (status != DFTI_NO_ERROR) {
        DftiFreeDescriptor(&hand);
        throw Math_error("DftiComputeBackward failed");
    }

    // Clean-up and return:
    DftiFreeDescriptor(&hand);
    return x;
}
#endif
