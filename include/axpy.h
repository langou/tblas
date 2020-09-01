//
//  axpy.h
//
//  Purpose
//  =======
//
//  Computes the sum of a constant multiple of a vector and another vector.
//
//      y <- alpha * x + y
//
//  Arguments
//  =========
//
//  n       length of vectors x and y
//
//  alpha   constant multiple of x
//
//  x       vector of length n
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//
//  y       vector of length n
//
//  incy    stride of vector y; if negative, y is stored in reverse order
//

#ifndef __axpy__
#define __axpy__

#include <cstddef>

using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void axpy(size_t n, T alpha, T *x, ptrdiff_t incx, T *y, ptrdiff_t incy)
    {
        if((incx==1)&&(incy==1))
        {
            for(size_t i=0;i<n;i++)
                y[i]+=alpha*x[i];
        }
        else
        {
            size_t ix=(incx>0)?0:(1-n)*incx;
            size_t iy=(incy>0)?0:(1-n)*incy;
            for(size_t i=0;i<n;i++)
            {
                y[iy]+=alpha*x[ix];
                ix+=incx;
                iy+=incy;
            }
        }
    }
}
#endif
