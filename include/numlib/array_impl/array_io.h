////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2017 Stig Rune Sellevag. All rights reserved.
//
// This code is licensed under the MIT License (MIT).
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SRS_ARRAY_IO_H
#define SRS_ARRAY_IO_H

#include <iomanip>
#include <iostream>


namespace srs {

// Non-member I/O operators for Array<T, N>:

template <class T>
std::ostream& operator<<(std::ostream& to, const Array<T, 1>& a)
{
    using size_type = typename Array<T, 1>::size_type;

    to << a.size() << '\n' << "[ ";
    for (size_type i = 0; i < a.size(); ++i) {
        to << std::setw(9) << a(i) << " ";
        if (!((i + 1) % 7) && (i != (a.size() - 1))) {
            to << "\n  ";
        }
    }
    to << ']';
    return to;
}

template <class T>
std::ostream& operator<<(std::ostream& to, const Array_ref<T, 1>& a)
{
    using size_type = typename Array_ref<T, 1>::size_type;

    to << a.size() << '\n' << "[ ";
    for (size_type i = 0; i < a.size(); ++i) {
        to << std::setw(9) << a(i) << " ";
        if (!((i + 1) % 7) && (i != (a.size() - 1))) {
            to << "\n  ";
        }
    }
    to << ']';
    return to;
}

template <class T>
std::istream& operator>>(std::istream& from, Array<T, 1>& a)
{
    using size_type = typename Array<T, 1>::size_type;

    size_type n;
    from >> n;
    if (n < 1) {
        throw Array_error("Array<T, 1>::operator>> : bad size");
    }
    a.resize(n);

    char ch;
    from >> ch;
    if (ch != '[') {
        throw Array_error("Array<T, 1>::operator>> : '[' missing");
    }
    for (size_type i = 0; i < n; ++i) {
        from >> a(i);
    }
    from >> ch;
    if (ch != ']') {
        throw Array_error("Array<T, 1>::operator>> : ']' missing");
    }
    return from;
}

template <class T>
std::ostream& operator<<(std::ostream& to, const Array<T, 2>& a)
{
    using size_type = typename Array<T, 2>::size_type;

    to << a.rows() << " x " << a.cols() << "\n[";
    for (size_type i = 0; i < a.rows(); ++i) {
        for (size_type j = 0; j < a.cols(); ++j) {
            to << std::setw(9) << a(i, j) << " ";
        }
        if (i != (a.rows() - 1)) {
            to << "\n ";
        }
    }
    to << "]\n";
    return to;
}

template <class T>
std::ostream& operator<<(std::ostream& to, const Array_ref<T, 2>& a)
{
    using size_type = typename Array_ref<T, 2>::size_type;

    to << a.rows() << " x " << a.cols() << "\n[";
    for (size_type i = 0; i < a.rows(); ++i) {
        for (size_type j = 0; j < a.cols(); ++j) {
            to << std::setw(9) << a(i, j) << " ";
        }
        if (i != (a.rows() - 1)) {
            to << "\n ";
        }
    }
    to << "]\n";
    return to;
}

template <class T>
std::istream& operator>>(std::istream& from, Array<T, 2>& a)
{
    using size_type = typename Array<T, 2>::size_type;

    size_type m;
    size_type n;
    char ch;

    from >> m >> ch >> n;
    if (ch != 'x') {
        throw Array_error("Array<T, 2>::operator>> : 'x' missing");
    }
    if ((m < 1) || (n < 1)) {
        throw Array_error("Array<T, 2>::operator>> : bad size");
    }
    a.resize(m, n);

    from >> ch;
    if (ch != '[') {
        throw Array_error("Array<T, 2>::operator>> : '[' missing");
    }
    for (size_type i = 0; i < m; ++i) {
        for (size_type j = 0; j < n; ++j) {
            from >> a(i, j);
        }
    }
    from >> ch;
    if (ch != ']') {
        throw Array_error("Array<T, 2>::operator>> : ']' missing");
    }
    return from;
}

template <class T>
std::ostream& operator<<(std::ostream& to, const Array<T, 3>& a)
{
    using size_type = typename Array<T, 3>::size_type;

    size_type n1 = a.rows();
    size_type n2 = a.cols();
    size_type n3 = a.depths();
    to << n1 << " x " << n2 << " x " << n3 << '\n' << "[ ";

    for (size_type k = 0; k < n3; ++k) {
        for (size_type i = 0; i < n1; ++i) {
            for (size_type j = 0; j < n2; ++j) {
                to << std::setw(9) << a(i, j, k) << " ";
            }
            if (i != (n1 - 1)) {
                to << "\n  ";
            }
        }
        if (k != (n3 - 1)) {
            to << "\n\n  ";
        }
    }
    to << ']';
    return to;
}

template <class T>
std::ostream& operator<<(std::ostream& to, const Array_ref<T, 3>& a)
{
    using size_type = typename Array_ref<T, 3>::size_type;

    size_type n1 = a.rows();
    size_type n2 = a.cols();
    size_type n3 = a.depths();
    to << n1 << " x " << n2 << " x " << n3 << '\n' << "[ ";

    for (size_type k = 0; k < n3; ++k) {
        for (size_type i = 0; i < n1; ++i) {
            for (size_type j = 0; j < n2; ++j) {
                to << std::setw(9) << a(i, j, k) << " ";
            }
            if (i != (n1 - 1)) {
                to << "\n  ";
            }
        }
        if (k != (n3 - 1)) {
            to << "\n\n  ";
        }
    }
    to << ']';
    return to;
}

template <class T>
std::istream& operator>>(std::istream& from, Array<T, 3>& a)
{
    using size_type = typename Array<T, 3>::size_type;

    size_type n1;
    size_type n2;
    size_type n3;
    char ch1;
    char ch2;

    from >> n1 >> ch1 >> n2 >> ch2 >> n3;
    if ((ch1 != 'x') || (ch2 != 'x')) {
        throw Array_error("Array<T, 3>::operator>> : 'x' missing");
    }
    if ((n1 < 1) || (n2 < 1) || (n3 < 1)) {
        throw Array_error("Array<T, 3>::operator>> : bad size");
    }
    a.resize(n1, n2, n3);

    from >> ch1;
    if (ch1 != '[') {
        throw Array_error("Array<T, 3>::operator>> : '[' missing");
    }
    for (size_type k = 0; k < n3; ++k) {
        for (size_type i = 0; i < n1; ++i) {
            for (size_type j = 0; j < n2; ++j) {
                from >> a(i, j, k);
            }
        }
    }
    from >> ch1;
    if (ch1 != ']') {
        throw Array_error("Array<T, 3>::operator>> : ']' missing");
    }
    return from;
}

}  // namespace srs

#endif  // SRS_ARRAY_IO_H
