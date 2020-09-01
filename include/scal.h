//
//  scal.h
//
//  Purpose
//  =======
//
//  Multiplies a vector x by a scalar alpha:
//
//      x <- alpha * x
//
//  Arguments
//  =========
//
//  n       length of vector x
//
//  alpha   scalar factor of x
//
//  x       vector of length n
//
//  incx    stride of vector x
//

#ifndef __scal__
#define __scal__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;

namespace tblas
{
    template <typename T>
    void scal(size_t n, T alpha, T *x, size_t incx=1)
    {
        for(size_t i=0;i<n;i++)
            x[i*incx]*=alpha;
    }
    
    template <typename T>
    void scal(size_t n, T alpha, complex<T> *x, size_t incx=1)
    {
        for(size_t i=0;i<n;i++)
            x[i*incx]*=alpha;
    }
}
#endif
