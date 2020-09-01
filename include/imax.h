//
//  imax.h
//
//  Purpose
//  =======
//
//  Find the index of the maximal element of a vector.
//
//  Returns
//  =======
//
//  index i such that |Re(x[i])| + |Im(x[i])| is maximal
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

#ifndef __imax__
#define __imax__

#include <cmath>
#include <cstddef>
#include <complex>
using std::complex;
using std::size_t;

namespace tblas
{
    template <typename T>
    size_t imax(size_t n, T *x, size_t incx=1)
    {
        using std::abs;
        size_t index=0;
        if(n>1)
        {
            T maxval=abs(x[0]);
            for(size_t i=1;i<n;i++)
            {
                if(abs(x[i*incx])>maxval)
                {
                    index=i;
                    maxval=abs(x[i*incx]);
                }
            }
        }
        return index;
    }
    
    template <typename T>
    size_t imax(size_t n, complex<T> *x, size_t incx=1)
    {
        using std::abs;
        size_t index=0;
        if(n>1)
        {
            T maxval=abs(real(x[0]))+abs(imag(x[0]));
            for(size_t i=1;i<n;i++)
            {
                T absval=abs(real(x[i*incx]))+abs(imag(x[i*incx]));
                if(absval>maxval)
                {
                    index=i;
                    maxval=absval;
                }
            }
        }
        return index;
    }
}
#endif
