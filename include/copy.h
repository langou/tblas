//
//  copy.h
//
//  Purpose
//  =======
//
//  Copies vector x to vector y.
//
//      y <- x
//
//  Arguments
//  =========
//
//  n       length of vectors x and y
//
//  x       vector of length n
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//
//  y       vector of length n
//
//  incy    stride of vector y; if negative, y is stored in reverse order
//

#ifndef __copy__
#define __copy__

#include <cstddef>

using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void copy(size_t n, T *x, ptrdiff_t incx, T *y, ptrdiff_t incy)
    {
        if((incx==1)&&(incy==1))
        {
            for(size_t i=0;i<n;i++)
                y[i]=x[i];
        }
        else
        {
            size_t ix=incx>0?0:(1-n)*incx;
            size_t iy=incy>0?0:(1-n)*incy;
            for(size_t i=0;i<n;i++)
            {
                y[iy]=x[ix];
                ix+=incx;
                iy+=incy;
            }
        }
    }
}
#endif
