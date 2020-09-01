//
//  spr.h
//
//  Purpose
//  =======
//
//  Performs the vector outer product
//
//      A <- alpha * x * x^T + A
//
//  where A is a symmetric matrix in packed storage format,
//  x is a vector, and alpha is a scalar.
//
//  Arguments
//  =========
//
//  uplo    specifies whether A is stored as upper ('U') or lower ('L') triangular
//
//  n       specifies the order of the symmetric matrix A
//
//  alpha   scalar multiple of the vector outer project
//
//  x       vector of length n
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//
//  A       triangular matrix of order n stored in packed format in a
//          vector of length n(n+1)/2
//

#ifndef __spr__
#define __spr__

#include <cstddef>

using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void spr(char uplo, size_t n, T alpha, T *x, ptrdiff_t incx, T *A)
    {
        if(incx==1)
        {
            if(uplo=='U')
            {
                for(size_t j=0;j<n;j++)
                {
                    T t=alpha*x[j];
                    for(size_t i=0;i<=j;i++)
                        A[i]+=x[i]*t;
                    A+=j+1;
                }
            }
            else if(uplo=='L')
            {
                for(size_t j=0;j<n;j++)
                {
                    T t=alpha*x[j];
                    for(size_t i=j;i<n;i++)
                        A[i-j]+=x[i]*t;
                    A+=n-j;
                }
            }
        }
        else
        {
            size_t kx=(incx>0)?0:(1-n)*incx;
            if(uplo=='U')
            {
                size_t jx=kx;
                for(size_t j=0;j<n;j++)
                {
                    T t=alpha*x[jx];
                    size_t ix=kx;
                    for(size_t i=0;i<=j;i++)
                    {
                        A[i]+=x[ix]*t;
                        ix+=incx;
                    }
                    A+=j+1;
                    jx+=incx;
                }
            }
            else if(uplo=='L')
            {
                size_t jx=kx;
                for(size_t j=0;j<n;j++)
                {
                    T t=alpha*x[jx];
                    size_t ix=jx;
                    for(size_t i=j;i<n;i++)
                    {
                        A[i-j]+=x[ix]*t;
                        ix+=incx;
                    }
                    jx+=incx;
                    A+=n-j;
                }
            }
        }
    }
}
#endif
