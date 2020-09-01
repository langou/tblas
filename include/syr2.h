//
//  syr2.h
//
//  Purpose
//  =======
//
//  Performs the symmetric vector outer product of vectors x and y,
//
//      A <- alpha * x * y^H + alpha * y * x^T + A
//
//  where A is a symmetric matrix, x and y are vectors, and alpha is a scalar.
//
//  Arguments
//  =========
//
//  uplo    specifies whether A is accessed as upper ('U') or lower ('L') triangular
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
//  A       symmetric matrix of order n
//
//  ldA     specifies the column length of A, must be at least n
//

#ifndef __syr2__
#define __syr2__

#include <cstddef>

using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void syr2(char uplo, size_t n, T alpha, T *x, ptrdiff_t incx, T *y, ptrdiff_t incy, T *A, size_t ldA)
    {
        const T zero(0.0);
        
        if((n==0)||(alpha==zero))
            return;
        
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
                    A+=ldA;
                }
            }
            else if(uplo=='L')
            {
                for(size_t j=0;j<n;j++)
                {
                    T tx=alpha*x[j];
                    T ty=alpha*y[j];
                    for(size_t i=j;i<n;i++)
                        A[i]+=x[i]*ty+y[i]*tx;
                    A+=ldA;
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
                    A+=ldA;
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
                        A[i]+=x[ix]*ty+y[iy]*tx;
                        ix+=incx;
                        iy+=incy;
                    }
                    jx+=incx;
                    jy+=incy;
                    A+=ldA;
                }
            }
        }
    }
}
#endif
