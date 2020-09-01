//
//  swap.h
//
//  Purpose
//  =======
//
//  Interchanges the elements of vectors x and y.
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

#ifndef __swap__
#define __swap__

#include <complex>
#include <cstddef>
#include <utility>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void swap(size_t n, T *x, ptrdiff_t incx, T *y, ptrdiff_t incy)
    {
        using std::swap;
        if((incx==1)&&(incy==1))
        {
            for(size_t i=0;i<n;i++)
                swap(x[i],y[i]);
        }
        else
        {
            size_t ix=incx>0?0:(1-n)*incx;
            size_t iy=incy>0?0:(1-n)*incy;
            for(size_t i=0;i<n;i++)
            {
                swap(x[ix],y[iy]);
                ix+=incx;
                iy+=incy;
            }
        }
    }
}
#endif
