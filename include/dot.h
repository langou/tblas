//
//  dot.h
//
//  Purpose
//  =======
//
//  Computes the dot product of two vectors in mixed precision.
//
//  Returns
//  =======
//
//  the dot product:  x^T y
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

#ifndef __dot__
#define __dot__

#include <cstddef>

using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T1,typename T2>
    T2 dot(size_t n, T2 sum, T1 *x, ptrdiff_t incx, T1 *y, ptrdiff_t incy)
    {
        if((incx==1)&&(incy==1))
        {
            for(size_t i=0;i<n;i++)
                sum+=static_cast<T2>(x[i])*static_cast<T2>(y[i]);
        }
        else
        {
            size_t ix=(incx>0)?0:(1-n)*(incx);
            size_t iy=(incy>0)?0:(1-n)*(incy);
            for(size_t i=0;i<n;i++)
            {
                sum+=static_cast<T2>(x[ix])*static_cast<T2>(y[iy]);
                ix+=incx;
                iy+=incy;
            }
        }
        return sum;
    }
}
#endif
