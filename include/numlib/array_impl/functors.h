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

#ifndef SRS_FUNCTORS_H
#define SRS_FUNCTORS_H


namespace srs {

template <class T>
struct Assign {
    void operator()(T& a, const T& c) { a = c; }
};

template <class T>
struct Add_assign {
    void operator()(T& a, const T& c) { a += c; }
};

template <class T>
struct Mul_assign {
    void operator()(T& a, const T& c) { a *= c; }
};

template <class T>
struct Minus_assign {
    void operator()(T& a, const T& c) { a -= c; }
};

template <class T>
struct Div_assign {
    void operator()(T& a, const T& c) { a /= c; }
};

template <class T>
struct Mod_assign {
    void operator()(T& a, const T& c) { a %= c; }
};

template <class T>
struct Or_assign {
    void operator()(T& a, const T& c) { a |= c; }
};

template <class T>
struct Xor_assign {
    void operator()(T& a, const T& c) { a ^= c; }
};

template <class T>
struct And_assign {
    void operator()(T& a, const T& c) { a &= c; }
};

template <class T>
struct Not_assign {
    void operator()(T& a) { a = !a; }
};

template <class T>
struct Not {
    void operator()(T& a) { !a; }
};

template <class T>
struct Unary_minus {
    void operator()(T& a) { -a; }
};

template <class T>
struct Complement {
    void operator()(T& a) { ~a; }
};

}  // namespace srs

#endif  // SRS_FUNCTORS_H
