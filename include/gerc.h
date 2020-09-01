//
//  gerc.h
//
//  Purpose
//  =======
//
//  Performs a vector outer product of complex vectors x and y,
//
//      A <- alpha * x * y^H + A
//
//  where A is a matrix, x and y are vectors, and alpha is a scalar.
//
//  Arguments
//  =========
//
//  m       specifies the number of rows of the matrix A
//
//  n       specifies the number of columns of the matrix A
//
//  alpha   scalar multiple of the vector outer project
//
//  x       vector of length m
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//
//  y       vector of length m if trans='N', or length n otherwise
//
//  incy    stride of vector y; if negative, y is stored in reverse order
//
//  A       matrix of size m-by-n matrix
//
//  ldA     specifies the column length of A, must be at least m
//


#ifndef __gerc__
#define __gerc__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void gerc(size_t m, size_t n, complex<T> alpha, complex<T> *x, ptrdiff_t incx, complex<T> *y, ptrdiff_t incy, complex<T> *A, size_t ldA)
    {
        const complex<T> zero(0.0);
        
        if((m==0)||(n==0)||(alpha==zero))
            return;
        
        size_t jy=(incy>0)?0:(1-n)*incy;
        size_t kx=(incx>0)?0:(1-m)*incx;
        
        if(incx==1)
        {
            for(size_t j=0;j<n;j++)
            {
                for(size_t i=0;i<m;i++)
                    A[i]+=x[i]*alpha*conj(y[jy]);
                jy+=incy;
                A+=ldA;
            }
        }
        else
        {
            for(size_t j=0;j<n;j++)
            {
                size_t ix=kx;
                for(size_t i=0;i<m;i++)
                {
                    A[i]+=x[ix]*alpha*conj(y[jy]);
                    ix+=incx;
                }
                A+=ldA;
                jy+=incy;
            }
        }
    }
}
#endif
