//
//  spr2.h
//
//  Purpose
//  =======
//
//  Performs the vector outer product of vectors x and y,
//
//      A <- alpha * x * y^T + alpha * y * x^T + A
//
//  where A is a symmetric matrix in packed format,
//  x and y are vectors, and alpha is a scalar.
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
//  y       vector of length n
//
//  incy    stride of vector y; if negative, y is stored in reverse order
//
//  A       triangular matrix of order n stored in packed format in a
//          vector of length n(n+1)/2
//

#ifndef __spr2__
#define __spr2__

#include <cstddef>

using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void spr2(char uplo, size_t n, T alpha, T *x, ptrdiff_t incx, T *y, ptrdiff_t incy, T *A)
    {
        if((incx==1)&&(incy==1))
        {
            if(uplo=='U')
            {
                for(size_t j=0;j<n;j++)
                {
                    T tx=alpha*x[j];
                    T ty=alpha*y[j];
                    for(size_t i=0;i<=j;i++)
                        A[i]+=x[i]*ty+y[i]*tx;
                    A+=j+1;
                }
            }
            else if(uplo=='L')
            {
                for(size_t j=0;j<n;j++)
                {
                    T tx=alpha*x[j];
                    T ty=alpha*y[j];
                    for(size_t i=j;i<n;i++)
                        A[i-j]+=x[i]*ty+y[i]*tx;
                    A+=n-j;
                }
            }
        }
        else
        {
            size_t kx=(incx>0)?0:(1-n)*incx;
            size_t ky=(incy>0)?0:(1-n)*incy;
            if(uplo=='U')
            {
                size_t jx=kx;
                size_t jy=ky;
                for(size_t j=0;j<n;j++)
                {
                    T tx=alpha*x[jx];
                    T ty=alpha*y[jy];
                    size_t ix=kx;
                    size_t iy=ky;
                    for(size_t i=0;i<=j;i++)
                    {
                        A[i]+=x[ix]*ty+y[iy]*tx;
                        ix+=incx;
                        iy+=incy;
                    }
                    A+=j+1;
                    jx+=incx;
                    jy+=incy;
                }
            }
            else if(uplo=='L')
            {
                size_t jx=kx;
                size_t jy=ky;
                for(size_t j=0;j<n;j++)
                {
                    T tx=alpha*x[jx];
                    T ty=alpha*y[jy];
                    size_t ix=jx;
                    size_t iy=jy;
                    for(size_t i=j;i<n;i++)
                    {
                        A[i-j]+=x[ix]*ty+y[iy]*tx;
                        ix+=incx;
                        iy+=incy;
                    }
                    jx+=incx;
                    jy+=incy;
                    A+=n-j;
                }
            }
        }
    }
}
#endif
