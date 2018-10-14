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

#ifndef SRS_ARRAY_REF_H
#define SRS_ARRAY_REF_H

namespace srs {

//
// N-dimensional dense array reference class.
//
// Note:
// - An Array_ref is a reference to memory in an Array specified by a slice.
//   An Array_ref does not owns its elements. As such, the created Array_ref is
//   not alias safe and does not take into account that the underlying array
//   memory could be freed (e.g. due to any operation involving a resize of the
//   array or that the array goes out of scope).
// - The general Array_ref template exists only to allow specializations.
//
template <class T, int N>
class Array_ref {
private:
    Array_ref();
};

}  // namespace srs

#include <srs/array_impl/array_ref1.h>
#include <srs/array_impl/array_ref2.h>
#include <srs/array_impl/array_ref3.h>

#endif  // SRS_ARRAY_REF_H
