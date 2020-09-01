//
//  nrm2.h
//
//  Purpose
//  =======
//
//  Compute the Euclidean norm of a vector.
//
//  Returns
//  =======
//
//  ||x||, the 2-norm of the vector x
//
//  Arguments
//  =========
//
//  n       length of the vector x
//
//  x       vector of length n
//
//  incx    stride of the vector x
//

#ifndef __nrm2__
#define __nrm2__

#include <complex>
#include <cmath>
#include <cstddef>

using std::complex;
using std::size_t;

namespace tblas
{
    template <typename T>
    T nrm2(size_t n, T *x, size_t incx=1)
    {
        using std::abs;
        using std::sqrt;
        const T one(1.0);
        const T zero(0.0);
        T scale=zero;
        T ssq=one;
        for(size_t i=0;i<n*incx;i+=incx)
        {
            if(x[i]!=zero)
            {
                T a=abs(x[i]);
                if(scale<a)
                {
                    T b=scale/a;
                    ssq=one+ssq*b*b;
                    scale=a;
                }
                else
                {
                    T b=a/scale;
                    ssq+=b*b;
                }
            }
        }
        return scale*sqrt(ssq);
    }
    
    template <typename T>
    T nrm2(size_t n, complex<T> *x, size_t incx=1)
    {
        using std::abs;
        using std::sqrt;
        const T one(1.0);
        const T zero(0.0);
        T scale=zero;
        T ssq=one;
        for(size_t i=0;i<n*incx;i+=incx)
        {
            if(real(x[i])!=zero)
            {
                T a=abs(real(x[i]));
                if(scale<a)
                {
                    T b=scale/a;
                    ssq=one+ssq*b*b;
                    scale=a;
                }
                else
                {
                    T b=a/scale;
                    ssq+=b*b;
                }
            }
            if(imag(x[i])!=zero)
            {
                T a=abs(imag(x[i]));
                if(scale<a)
                {
                    T b=scale/a;
                    ssq=one+ssq*b*b;
                    scale=a;
                }
                else
                {
                    T b=a/scale;
                    ssq+=b*b;
                }
            }
        }
        return scale*sqrt(ssq);
    }
}
#endif
