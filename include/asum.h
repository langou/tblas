//
//  asum.h
//
//  Purpose
//  =======
//
//  Sums absolute values (of real and imaginary parts) of a vector.
//
//  Returns
//  =======
//
//  abs(x[0]) + ... + abs(x[n-1]))
//
//      or
//
//  abs(Re(x[0])) + abs(Im(x[0])) + ... + abs(Re(x[n-1])) + abs(Im(x[n-1]))
//
//  Arguments
//  =========
//
//  n       length of vector x
//
//  x       vector of length n
//
//  incx    stride of vector x
//

#ifndef __asum__
#define __asum__

#include <complex>
#include <cstddef>
#include <cmath>

using std::complex;
using std::size_t;

namespace tblas
{
    template <typename T>
    T asum(size_t n, T *x, size_t incx)
    {
        using std::abs;
        T sum(0.0);
        for(size_t i=0;i<n*incx;i+=incx)
            sum+=abs(x[i]);
        return sum;
    }

    template <typename T>
    T asum(size_t n, complex<T> *x, size_t incx)
    {
        using std::abs;
        T sum(0.0);
        for(size_t i=0;i<n*incx;i+=incx)
            sum+=abs(real(x[i]))+abs(imag(x[i]));
        return sum;
    }
}
#endif
