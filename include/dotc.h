//
//  dotc.h
//
//  Purpose
//  =======
//
//  Computes the dot product of the conjugate of a vector and another vector
//
//  Returns
//  =======
//
//  the dot product:  x^H y
//
//  Arguments
//  =========
//
//  n       length of vectors x and y
//
//  sum     initial value of accumulation
//
//  x       vector of length n
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//
//  y       vector of length n
//
//  incy    stride of vector y; if negative, y is stored in reverse order
//

#ifndef __dotc__
#define __dotc__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    complex<T> dotc(size_t n, complex<T> sum, complex<T> *x, ptrdiff_t incx, complex<T> *y, ptrdiff_t incy)
    {
        if((incx==1)&&(incy==1))
        {
            for(size_t i=0;i<n;i++)
                sum+=conj(x[i])*y[i];
        }
        else
        {
            size_t ix=incx>0?0:(1-n)*incx;
            size_t iy=incy>0?0:(1-n)*incy;
            for(size_t i=0;i<n;i++)
            {
                sum+=conj(x[ix])*y[iy];
                ix+=incx;
                iy+=incy;
            }
        }
        return sum;
    }
}
#endif
