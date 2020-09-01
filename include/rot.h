//
//  rot.h
//
//  Purpose
//  =======
//
//  Applies real plane rotations to a series of points (x,y).
//
//      (x,y) <- (c * x + s * y, -s * x + c * y)
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
//  c       scalar representing cos(theta) for rotation angle theta
//
//  s       scalar representing sin(theta) for rotation angle theta
//

#ifndef __rot__
#define __rot__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void rot(size_t n, T *x, ptrdiff_t incx, T *y, ptrdiff_t incy, T c, T s)
    {
        if((incx==1)&&(incy==1))
        {
            for(size_t i=0;i<n;i++)
            {
                T temp=c*x[i]+s*y[i];
                y[i]=c*y[i]-s*x[i];
                x[i]=temp;
            }
        }
        else
        {
            size_t ix=incx>0?0:(1-n)*incx;
            size_t iy=incy>0?0:(1-n)*incy;
            for(size_t i=0;i<n;i++)
            {
                T temp=c*x[ix]+s*y[iy];
                y[iy]=c*y[iy]-s*x[ix];
                x[ix]=temp;
                ix+=incx;
                iy+=incy;
            }
        }
    }
    
    template <typename T>
    void rot(size_t n, complex<T> *x, ptrdiff_t incx, complex<T> *y, ptrdiff_t incy, T c, T s)
    {
        if((incx==1)&&(incy==1))
        {
            for(size_t i=0;i<n;i++)
            {
                complex<T> temp=c*x[i]+s*y[i];
                y[i]=c*y[i]-s*x[i];
                x[i]=temp;
            }
        }
        else
        {
            size_t ix=incx>0?0:(1-n)*incx;
            size_t iy=incy>0?0:(1-n)*incy;
            for(size_t i=0;i<n;i++)
            {
                complex<T> temp=c*x[ix]+s*y[iy];
                y[iy]=c*y[iy]-s*x[ix];
                x[ix]=temp;
                ix+=incx;
                iy+=incy;
            }
        }
    }
}
#endif
