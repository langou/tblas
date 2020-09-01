//
//  hpr2.h
//
//  Purpose
//  =======
//
//  Performs the vector outer product of vectors x and y,
//
//      A <- alpha * x * y^H + conj(alpha) * y * x^H + A
//
//  where A is a Hermitian matrix in packed format,
//  x and y are vectors, and alpha is a scalar.
//
//  Arguments
//  =========
//
//  uplo    specifies whether A is stored as upper ('U') or lower ('L') triangular
//
//  n       specifies the order of the Hermitian matrix A
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
//  A       Hermitian matrix of order n stored in packed format in a
//          vector of length n(n+1)/2
//

#ifndef __hpr2__
#define __hpr2__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void hpr2(char uplo, size_t n, complex<T> alpha, complex<T> *x, ptrdiff_t incx, complex<T> *y, ptrdiff_t incy, complex<T> *A)
    {
        const complex<T> zero(0.0);
        
        if((n==0)||(alpha==zero))
            return;
        
        if((incx==1)&&(incy==1))
        {
            if(uplo=='U')
            {
                for(size_t j=0;j<n;j++)
                {
                    complex<T> tx=conj(alpha*x[j]);
                    complex<T> ty=alpha*conj(y[j]);
                    for(size_t i=0;i<j;i++)
                        A[i]+=x[i]*ty+y[i]*tx;
                    A[j]=real(A[j])+real(x[j]*ty+y[j]*tx);
                    A+=j+1;
                }
            }
            else if(uplo=='L')
            {
                for(size_t j=0;j<n;j++)
                {
                    complex<T> tx=conj(alpha*x[j]);
                    complex<T> ty=alpha*conj(y[j]);
                    A[0]=real(A[0])+real(x[j]*ty+y[j]*tx);
                    for(size_t i=j+1;i<n;i++)
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
                    complex<T> tx=conj(alpha*x[jx]);
                    complex<T> ty=alpha*conj(y[jy]);
                    size_t ix=kx;
                    size_t iy=ky;
                    for(size_t i=0;i<j;i++)
                    {
                        A[i]+=x[ix]*ty+y[iy]*tx;
                        ix+=incx;
                        iy+=incy;
                    }
                    A[j]=real(A[j])+real(x[jx]*ty+y[jy]*tx);
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
                    complex<T> tx=conj(alpha*x[jx]);
                    complex<T> ty=alpha*conj(y[jy]);
                    size_t ix=jx;
                    size_t iy=jy;
                    A[0]=real(A[0])+real(x[jx]*ty+y[jy]*tx);
                    for(size_t i=j+1;i<n;i++)
                    {
                        ix+=incx;
                        iy+=incy;
                        A[i-j]+=x[ix]*ty+y[iy]*tx;
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
